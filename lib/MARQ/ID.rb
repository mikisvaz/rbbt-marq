require 'MARQ'
require 'rbbt/sources/organism'
require 'DBcache'

module ID

  #DEFAULT_FORMAT_ALL = 'Entrez Gene Id'
  DEFAULT_FORMAT_ALL = nil
  DEFAULT_FORMATS = {}

  def self.id_position(org, id)
    @@supported[org] ||= Organism.supported_ids(org)
    Organism.id_position(@@supported[org], id)
  end


  def self.DB_load(org, to = 0)
    identifier_file = File.join(Rbbt.datadir,'organisms',org,'identifiers')
    tablename = "ID_#{org.to_s.strip}_#{to.to_s.strip}"
    
    if DBcache.has_table?(tablename)
      DBcache.drop(tablename)
    end

    DBcache.create(tablename, 'CHAR(50)', ['CHAR(50)'])

    total_fields = Organism.supported_ids(org).length
    
    total_fields.times{|field|
      File.open(identifier_file).each{|l|
        next if l =~ /^#/
        codes = l.chomp.split(/\t/)

        native = codes[to]
        next if native.nil? || native == "" 

        other  = codes[field]
        next if other.nil? || other == ""
        

        #codes.collect{|c| c.split("|")}.flatten.compact.select{|c| c != ""}.uniq.each{|code|
        other.split("|").each{|code|
          begin
            DBcache.fast_add(tablename, code.downcase, [native])
          rescue
            puts $!.message
          end
        }
      }
    }
  end

  def self.translate_DB(org, genes, options = {})
    to = options[:to] || DEFAULT_FORMATS[org] || DEFAULT_FORMAT_ALL

    if to
      to = id_position(org, to)
    else 
      to = 0
    end

    tablename = "ID_#{org.to_s.strip}_#{to.to_s.strip}"
    DB_load(org, to) unless DBcache.has_table?(tablename)
    genes = genes.collect{|gene| gene.strip.downcase}
    DBcache.load(tablename, genes).values_at(*genes).collect{|gene| gene.first if gene}
  end
  
  def self.translate_index(org, genes, options = {})
    genes = genes.collect{|gene| gene.strip if gene}
    to = options[:to] || DEFAULT_FORMATS[org] || DEFAULT_FORMAT_ALL
    from = options[:from]
    @indexes ||= {}
    if  @indexes[org.to_s + to.to_s + from.to_s].nil?
      puts "Loading #{ org }"
      options = {}
      options[:other] = [from] if from
      options[:native] = to if to
      options[:case_sensitive] = false
      
      @indexes[org.to_s + to.to_s + from.to_s] = Organism.id_index(org, options)
    end

    index =  @indexes[org.to_s + to.to_s + from.to_s]
    
    return genes.collect{|code| 
      code.nil? || code == "" ? nil : index[code]
    }
  end

  @@supported = {}
  def self.translate_grep(org, genes, options = {})
    to = options[:to] || DEFAULT_FORMATS[org] || DEFAULT_FORMAT_ALL
    from = options[:from]

    options = {}
    options[:extra]  = [id_position(org, from)] if from
    options[:native] = id_position(org,to) if to
    options[:case_sensitive] ||= false

    genes = genes.collect{|gene| gene.strip if gene}
    genes = genes.collect{|gene| gene.downcase} unless options[:case_sensitive] 

    #genes_re = '\(^\|[' + "\t" + '|]\)' + genes.join('\|') + '\($\|[' + "\t" + '|]\)'
    #cmd = "cat #{File.join(Rbbt.datadir,'organisms',org,'identifiers')}|grep -i '#{genes_re}'"
    genes_re = '(?:^|[\t\|])(?:' + genes.join('|') + ')(?:$|[\t\|])'
    cmd = "cat #{File.join(Rbbt.datadir,'organisms',org,'identifiers')}|ruby -lne 'puts $_ if $_ =~ /#{genes_re}/i'"

    index = Index.index(IO::popen(cmd), options)

    index.values_at(*genes)

  end


  def self.AILUN_translate(platform, genes)
    index = Open.to_hash("ftp://ailun.stanford.edu/ailun/annotation/geo/#{platform}.annot.gz", :fix => proc{|l| l.match(/^(.*?)\t(.*?)\t.*/); $1.downcase + "\t" + $2 })
    index.values_at(*genes.collect{|code| code.strip.downcase}).collect{|v| v.nil? ? nil : v.first.first}
  end

  class << self
    alias_method :translate, :translate_DB
  end

end

if __FILE__ == $0
  require 'benchmark'


  Organism.all.each{|org|
    ID.DB_load(org)
  }


  #num = 4000

  #genes = File.open(Rbbt.datadir + '/organisms/human/identifiers').collect{|l| l.split("\t")[2] }.select{|c| c != ""}[1000,num]

  #trans = nil
  #puts Benchmark.measure{
  #  trans = ID.translate('human', genes)
  #}
  #p trans
  #
end
