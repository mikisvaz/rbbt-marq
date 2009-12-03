require 'MARQ'
require 'rbbt/util/open'
require 'rbbt/sources/organism'

module GEO

  CACHE_DIR = File.join(MARQ.cachedir,'GEO')
  FileUtils.mkdir_p CACHE_DIR unless File.exists? CACHE_DIR

  GEO_SOFT="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&view=full&form=text&acc="
  def self.get_soft(item)
    item = item.strip
    cache_file = File.join(CACHE_DIR, item + '.soft')
    if File.exist?( cache_file )
      File.open(cache_file).read
    else
      content = Open.read(GEO_SOFT + item, :nocache => true)
      fout = File.open(cache_file,'w')
      fout.write content
      fout.close
      content
    end
  end

  #{{{ Eutils
  module Eutils
    def self.organism_platforms(org)
      name = Organism.name(org)
      Open.read("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=\"#{name}\"[Organism:exp]+AND+%22gpl%22[Filter]&retmax=10000").
        scan(/<Id>(\d+?)<\/Id>/).collect{|id| id.first}.collect{|id| "GPL#{id.sub(/^100*/,'')}"}
    end

    def self.GPL_datasets(platform)
      Open.read("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=#{platform}[Accession]&retmax=2000").
      scan(/<Id>(\d+?)<\/Id>/).collect{|id| id.first}.select{|id| !id.match(/^(1|2)000/) }.collect{|id| "GDS#{id}"}
    end

    def self.GSE_dataset?(gse)
      Open.read("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=geo&term=#{gse}[Accession]&retmax=2000").
      match(/<Id>(\d+?)<\/Id>/) != nil
    end

  end



  #{{{ Helper functions
  

  def self.consecutive?(ids)
    ids.collect{|id| id.to_i}.sort[0..19] == (1..20).to_a
  end

  def self.numerical?(ids)
    ids.compact.select{|id| ! id.match(/^\d+$/)}.uniq.length < ids.length.to_f / 10
  end

  def self.dna_sequence?(ids)
    ids.compact.select{|id| ! id.strip.match(/^[ATCG]+$/i)}.empty?
  end


  ID_FIX = {
    :mgi_unigene => proc{|gene| if gene then gene.match(/^Mm./) ? gene : "Mm." + gene end},
    :human_unigene => proc{|gene| if gene then gene.match(/^Hs./) ? gene : "Hs." + gene end},
  }

  @@formats = {}
  def self.guessIds(genes,org, name = nil)
    @@formats[org] ||= Organism.id_formats(org)
    if consecutive?(genes) || dna_sequence?(genes) || (numerical?(genes) && (name.nil? || !name.match(/entrez/i)))
      id = nil
    else
      fix = ID_FIX[(org + "_" + name.downcase).to_sym] if name
      if fix
        genes = genes.collect{|gene| fix.call(gene)}
      end
      id = Organism.guessIdFormat(@@formats[org], genes)
    end
    
    id 
  end

  @@r = nil
  def self.r
    if @@r.nil?

      # FIXME: RSruby does not install very well, this require id hidden here.
      require 'rsruby'

      RSRuby.instance.source(MARQ.rootdir + '/R/MA.R')
      RSRuby.instance.source(MARQ.rootdir + '/R/GEO.R')
      @@r = RSRuby.instance
    end
    @@r
  end


  #{{{ Process

  def self.get_GPL(name, prefix, id_field = nil)
    r.GEO_GPL_process(name, prefix, id_field, CACHE_DIR)
  end

  def self.get_GDS(name, prefix, id_field = nil, id_file = nil)
    r.GEO_GDS_process(name, prefix, id_field, id_file, CACHE_DIR)
  end

  def self.get_GSE(gsms, conditions, do_log, prefix, id_file = nil, fields= nil, title = nil, description = nil)
    r.GEO_GSE_process(gsms, conditions, prefix, do_log, id_file, fields, title, description, CACHE_DIR)
  end

  def self.GSE_info(series)
    soft = get_soft(series)
    raise "SOFT file error" if soft !~ /!/

    if match = soft.scan(/!Series_platform_id\s*=?\s*(.*)/)
      platform = match.flatten.collect{|p| p.strip}
    else
      raise "No Platform information" 
    end

    if soft.match(/!Series_title \s*=?\s*(.*)/)
      title = $1
    else
      raise "No Title information" 
    end

    if soft.match(/!Series_summary \s*=?\s*(.*)/)
      matches = soft.scan(/!Series_summary \s*=?\s*(.*)/).to_a
      description = matches.collect{|m| m.to_s.strip.sub(/!Series_summary \s*=?\s*/,'')}.join("\n")
    else
      raise "No Summary information" 
    end

    if soft.match(/!Series_sample_id \s*=?\s*(.*)/)
      matches = soft.scan(/!Series_sample_id \s*=?\s*(.*)/).to_a
      samples = matches.collect{|m| m.to_s.strip.sub(/!Series_sample_id \s*=?\s*/,'')}
    else
      raise "No Summary information" 
    end

    {
      :platform => platform.join("_"),
      :description =>description.strip,
      :title => title.strip,
      :samples => samples,
    }
  end

  def self.GSM_info(array)
    soft = get_soft(array)

    if soft.match(/!Sample_title\s*=?\s*(.*)/)
      title = $1
    else
      raise "No Title information" 
    end


    if soft.match(/!Sample_description \s*=?\s*(.*)/)
      description = $1
    else
      raise "No Description information" 
    end

    {
      
      :description =>description.strip,
      :title => title.strip,
    }
  end

  def self.GPL_id_fields(platform)
    soft = get_soft(platform)
    data = soft.split(/!platform_table_begin/s)[1].collect{|l| l.chomp.split(/\t/)}
    data.shift
    data.shift
  end
  
  def self.GPL_info(platform)
    if !File.exist?(File.join(MARQ.datadir, 'GEO', 'platforms',"#{platform}.yaml"))  &&
       !File.exist?(File.join(MARQ.datadir, 'GEO', 'platforms',"#{platform}.skip")) 
      begin
        if platform =~ /_/
          organism =  GPL_info(platform.match(/(.*?)_/)[1])[:organism]
          
          info =  {
            :organism => organism,
            :title => "Merged platforms #{ platform }",
          }
          return info
        end
        soft = get_soft(platform)


        raise "SOFT file error" if soft !~ /!/

        organisms = soft.scan(/!Platform_organism\s*=\s*(.*)/).collect{|v| v.first.strip}

        if organisms.empty?
          raise "No Organism information" 
        else
          # This might happen actually GPL2529
          organisms.delete('Schizosaccharomyces pombe') if  organisms.include?('Saccharomyces cerevisiae')
          org_name = organisms.first
        end


        title = ""
        if soft.match(/!Platform_title\s*=\s*(.*)/)
          title = $1
        end

        org = Organism.name2org(org_name)
        raise "Organism not identified" if org.nil?

        if soft.match(/!platform_table_begin/)
          data = soft.split(/!platform_table_begin/s)[1].collect{|l| l.chomp.split(/\t/)}
          data.shift
          names = data.shift
          total = data.first.length
          genes = data.sort_by{ rand }[1..1000].collect{|v| v.first}

          id = guessIds(genes,org, names.first)
          other = nil
          other_pos = 0
          other_count = 0
          other_name = 0
          if id.nil? 
            (1..total - 1).to_a.each{|num|
              genes = data.collect{|v| v[num]}
              other = guessIds(genes,org, name = names[num])

              if other && other[1] > other_count
                other_pos = num
                other_count = other[1]
                other_name = names[num]
              end
            }
          end
        else
          raise "Soft file incomplete"
        end

        info = {:organism => org, :BioMart_ID => id ? id.first : nil, :title => title }
        info[:other_ID_field] = [other_pos + 1, other_name] if other_pos > 0


        Open.write(File.join(MARQ.datadir, 'GEO', 'platforms',"#{platform}.yaml"), info.to_yaml)
      rescue Exception
        puts $!.message
        puts $!.backtrace
        Open.write(File.join(MARQ.datadir, 'GEO', 'platforms',"#{platform}.skip"), $!.message)
      end
    end

    raise "Platform info for #{ platform } is not available and could not be automatically produced." if File.exist?(File.join(MARQ.datadir, 'GEO', 'platforms',"#{platform}.skip"))

    YAML::load(File.open(File.join(MARQ.datadir, 'GEO', 'platforms',"#{platform}.yaml")))
  end

  def self.GDS_info(name)
    begin
      title, description = Open.read(dataset_path(name) + '.description').split(/\n--\n/).values_at(0,1)
      {:title => title.strip, :description => description.strip}
    rescue Exception
      puts $!.message
      {:title => "" , :description => "" }
    end

  end


  #{{{ Misc Info

  def self.clean(name)
    name.sub(/_cross_platform/,'') if name
  end

  def self.platform_path(platform)
    File.join(MARQ.datadir, "GEO/#{clean(platform)}")
  end

  def self.dataset_path(dataset, platform = nil)
    if platform
      return Dir.glob(File.join(platform_path(clean(platform)),"/*/#{ dataset }")).first.match(/(.*)\./)[1]
    else
      files = Dir.glob(File.join(MARQ.datadir, "GEO/GPL*/*/#{ dataset }.*"))
      if files.any?
        return files.first.match(/(.*)\./)[1]
      else
        return ""
      end
    end
  end

  def self.is_cross_platform?(dataset)
    dataset =~ /_cross_platform/
  end

  def self.has_cross_platform?(dataset = nil, platform = nil)
    platform = clean(platform)
    raise "Dataset #{ dataset } not found" if dataset && dataset_path(dataset, platform).nil? 
    raise "Platform #{ platform } not found" if platform && platform_path(platform).nil? 
    if dataset
      File.exists?(dataset_path(dataset, platform) + "_cross_platform.orders")
    else
      Dir.glob(File.join(platform_path(platform), '*', '*_cross_platform.orders')).any?
    end
  end


  def self.platform_datasets(platform)
    Dir.glob(File.join(platform_path(platform),"*/*.orders")).collect{|f| File.basename(f).sub(/.orders$/,'')}.select{|d| !is_cross_platform?(d)}
  end

  def self.dataset_platform(dataset)
    dataset_path(dataset).match(/(GPL\d+)/)
    $1
  end

  def self.organism_platforms(organism)
    Dir.glob(File.join(MARQ.datadir, "GEO/GPL*")).collect{|f| 
      File.basename(f)
    }.select{|platform| 
      GPL_info(platform)[:organism] == organism &&
      platform_datasets(platform).any?
    }
  end

  #{{{ Processing

  def self.process_GDS(dataset, platform, field = nil)
    puts "Processing GDS #{ dataset }. Platform #{ platform }"

    puts "-- Original"
    prefix = File.join(platform_path(platform), 'GDS', dataset.to_s)
    GEO.get_GDS(dataset, prefix, field, nil)

    # Was there an error?
    if File.exist?(prefix + '.skip')
      FileUtils.cp(prefix + '.skip', prefix + '_cross_platform.skip')
      return
    end

    if File.exist?(File.join(platform,'cross_platform'))
      puts "-- Translated to cross_platform format"
      GEO.get_GDS(dataset, prefix + '_cross_platform', field, File.join(platform_path(platform), 'translations')) 
    end
  end

  # Rearange the lines of a file with the given order. The order specifies, for
  # each position in the original file, where it should en in the final file
  def self.rearange(order, file, missing = "NA")
    orig_lines = []
    File.open(file).each_line{|l| orig_lines << l}

    return if orig_lines.empty?
    columns = orig_lines.first.split(/\t/).length
    
    lines = Array.new(order.length)

    orig_lines.each_with_index{|l,i|
      next if order[i].nil?
      lines[order[i]] = l.chomp
    }

    lines = lines.collect{|l| l || [missing]*columns*"\t"}

    fout = File.open(file, 'w')
    fout.puts(lines.join("\n"))
    fout.close
  end

  # Fix possible discrepancies in ids between series and platforms
  def self.fix_GSE_ids(platform_codes_file, prefix)
    platform_codes = File.open(platform_codes_file).collect{|l| l.chomp}
    platform_order = {}
    
    platform_codes.each_with_index{|code, i|
      platform_order[code] = i
    }

    series_codes = File.open(prefix + '.codes').collect{|l| l.chomp}

    platform_positions = platform_order.values_at(*series_codes)

    # Fill with nil for missing positions
    platform_positions[platform_codes.length - 1] ||= nil

    %w(t logratios orders pvalues).each{|ext|
      rearange(platform_positions, prefix + '.' + ext)
    }

    Open.write(prefix + '.swap', platform_positions.join("\n"))
  end



  def self.process_GSE(series, info)
    return if Dir.glob(File.join(info[:platform], 'GSE', series) + '.*').any?

    gsms = []
    conditions = {}
    info[:arrays].each{|gsm, cond|
      gsms << gsm
      cond.each{|condition, value|
        conditions[condition] ||= []
        conditions[condition] << value
      }
    }
    platform = info[:platform]
    do_log = nil
    do_log = !info[:log2] if info[:log2]
    fields = info[:fields]

    puts "Processing GSE #{ series }. Platform #{ platform }"

    prefix = File.join(platform_path(platform), 'GSE', series.to_s)
    puts "-- Original"
    GEO.get_GSE(gsms, conditions, do_log, prefix, nil, fields, info[:title], info[:description])

    # Was there an error?
    if File.exist?(prefix + '.skip')
      FileUtils.cp(prefix + '.skip', prefix + '_cross_platform.skip')
      return
    end

    if platform =~ /_/
      FileUtils.cp(prefix + '.codes', File.join(platform_path(platform),'codes'))
      codes  = Open.read(File.join(platform_path(platform), 'codes')).collect{|l| l.chomp}
      organism = GEO::GPL_info(platform.match(/(.*?)_/)[1])[:organism]
      translations = ID.translate(organism, codes) 
      Open.write(File.join(platform_path(platform), 'translations'), translations.collect{|v| v || "NO MATCH"}.join("\n"))
      Open.write(File.join(platform_path(platform), 'cross_platform'), translations.compact.sort.uniq.join("\n"))
    else
      # Are the codes of the series equivalent to the ones in the platform?
      if  File.open(File.join(platform_path(platform),'codes')).collect{|l| l.chomp} != File.open(prefix + '.codes').collect{|l| l.chomp}
        fix_GSE_ids(File.join(platform_path(platform), 'codes'),prefix);
        FileUtils.cp(File.join(platform_path(platform), 'codes'),prefix + '.codes')

      end
    end


    if File.exist?(File.join(platform,'translations'))
      FileUtils.cp(File.join(platform,'translations'), prefix + '.translations')
      if File.exist?(prefix + '.swap')
        orders = Open.read(prefix + '.swap').collect{|l| l.chomp}
        inverse_orders = Array.new(orders.length)
        orders.each_with_index{|pos,i|
          next if pos !~ /\d/
          inverse_orders[pos.to_i] = i
        }
        rearange(inverse_orders, prefix + '.translations', "NO MATCH")
      end
      puts "-- Translated to cross_platform format"
      GEO.get_GSE(gsms, conditions, do_log, prefix + '_cross_platform', prefix + '.translations',fields, info[:title], info[:description])
      fix_GSE_ids(File.join(platform_path(platform), 'cross_platform'),prefix + '_cross_platform');
      FileUtils.cp(File.join(platform_path(platform), 'cross_platform'),prefix + '_cross_platform.codes')
      FileUtils.rm(prefix + '.translations') if File.exist?(prefix + '.translations')
    end
     FileUtils.rm(prefix + '.swap') if File.exist?(prefix + '.swap')
  end

  def self.process_platform(platform)
    path = platform_path(platform)
    return if File.exist? path

    if platform =~ /_/
      FileUtils.mkdir(path)
      FileUtils.mkdir(path + '/GSE')
      FileUtils.mkdir(path + '/GDS')
      return
    end

    info = GEO::GPL_info(platform)
    organism = info[:organism]

    field = info[:other_ID_field]
    id = info[:BioMart_ID]
    org = info[:organism]
    field = nil if field == ""
    id = nil if id == ""
 

    puts "Processing Platform #{ platform }"
    [platform,
      File.join(platform_path(platform), 'GDS'),
      File.join(platform_path(platform), 'GSE'),
    ].each{|d|
      FileUtils.mkdir d unless File.exist? d
    }

    get_GPL(platform, platform_path(platform), nil)
    FileUtils.mv platform_path(platform) + '.codes', File.join(platform_path(platform), 'codes')
    

    # AILUN translations
    codes  = Open.read(File.join(platform_path(platform), 'codes')).collect{|l| l.chomp}
    ailun = ID.AILUN_translate(platform, codes)
    Open.write(File.join(platform_path(platform), 'ailun'), ailun.collect{|v| v || "NO MATCH"}.join("\n")) if ailun.compact.length > codes.length.to_f / 10

    # BioMart translations
    biomart = []
    if id || field
      if id
        codes  = Open.read(File.join(platform_path(platform), 'codes')).collect{|l| l.chomp}
      else
        if field
          get_GPL(platform, platform_path(platform), field[0])
          FileUtils.mv platform_path(platform) + '.codes', File.join(platform_path(platform), 'other')
        end

        fix = ID_FIX[(organism + "_" + field[1].downcase).to_sym]
        codes  = Open.read(File.join(platform_path(platform), 'other')).collect{|l| 
          code = l.chomp
          code = fix.call(code) if fix
          code
        }
      end

      biomart = ID.translate(organism, codes) 
      Open.write(File.join(platform_path(platform), 'biomart'), biomart.collect{|v| v || "NO MATCH"}.join("\n")) if biomart.compact.length > codes.length.to_f / 10
    end

    # Select Best and save
    translations = []
    if ailun.compact.uniq.length > biomart.compact.uniq.length
      id_type = ID::DEFAULT_FORMATS[organism] || ID::DEFAULT_FORMAT_ALL || id || field || "Entrez Gene Id"
      if id_type.to_s !~ /Entrez/i
        translations = ID.translate(org,ailun.collect{|gene| gene || "NO MATCH"}) 
      else
        translations = ailun 
      end
    else
      translations = biomart
    end

    if translations.compact.length > codes.length.to_f / 10        
      Open.write(File.join(platform_path(platform), 'translations'), translations.collect{|v| v || "NO MATCH"}.join("\n"))
      Open.write(File.join(platform_path(platform), 'cross_platform'), translations.compact.sort.uniq.join("\n"))
    end

  end


  def self.process_platform_datasets(platform, force = false)
    raise "Platform #{ platform } not ready" unless File.exist? platform_path(platform)

    info = YAML::load(File.open(File.join(MARQ.datadir, "GEO/platforms/#{platform}.yaml")))

    datasets = GEO::Eutils::GPL_datasets(platform)
    datasets.each{|dataset|
      next if Dir.glob(File.join(platform_path(platform), 'GDS', dataset) + '.*').any? && ! force
      process_GDS(dataset, platform, nil) 
    }
  end

end

if __FILE__ == $0

  p GEO.GPL_info('GPL920_GPL927')
  p GEO.GPL_id_fields('GPL920')
  puts GEO.GSE_info('GSE962')
  puts GEO.GSE_info('GSE8982')
  puts GEO::Eutils.GSE_dataset?('GSE8982')
  puts GEO::Eutils.GSE_dataset?('GSE962')

  exit

  #puts GEO::dataset_path('GDS1103').inspect 
  #puts GEO::dataset_platform('GDS1103').inspect 

  #  puts GEO.dataset_path('GDS2931')
  #  puts GEO.platform_datasets('GPL91')
  #  puts GEO.platform_datasets('GPL91').select{|d| GEO.has_cross_platform?(d)}
  #
  #  gpls = Open.read('ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/PLATFORMS.txt').collect{|l|
  #    l.chomp.split.first
  #  }
  #
  #  %w(GPL85).each{|gpl|
  #    puts gpl
  #    puts GEO::GPL_info(gpl).inspect if gpl =~ /GPL/
  #  }
  #
  #puts GEO::GSM_info('GSM70604').inspect 

  p GEO::Eutils.organism_platforms('human')

end
