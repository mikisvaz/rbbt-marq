require 'rbbt/util/open'
require 'MARQ'
require 'MARQ/ID'

module CustomDS
  @@r = nil

  def self.customdir
    File.join(MARQ.datadir,'CustomDS')
  end

  def self.r
    require 'rsruby'
    if @@r.nil?
      RSRuby.instance.source(MARQ.rootdir + '/R/MA.R')
      RSRuby.instance.source(MARQ.rootdir + '/R/CustomDS.R')
      @@r = RSRuby.instance
    end
    @@r
  end

  def self.info(dataset)
    begin
      text =  Open.read(path(dataset) + '.description')
      if text =~ /(.*)\n--\n(.*)/
        {:title => $1.strip, :description => $2.strip}
      elsif text.collect.length > 1 || text.length > 200
        {:title => "", :description => text}
      else
        {:title => text, :description => ""}
      end
    rescue Exception
      puts $!.message
      {:title => "" , :description => "" }
    end
  end

  def self.path(dataset)
    files = Dir.glob(customdir + "/*/#{ dataset }.orders")
    if files.length == 1
      files.first.sub(/.orders/,'')
    else
      Dir.glob(customdir + "/*/#{ dataset }").first
    end
  end

  def self.organism(dataset)
    path(dataset).match(/#{ customdir }\/(.*?)\//)[1]
  end

  def self.is_cross_platform?(dataset)
    dataset.match(/_cross_platform/)
  end

  def self.clean(dataset)
    dataset.sub(/_cross_platform/,'')
  end
  
  def self.has_cross_platform?(dataset)
    Dir.glob(path(clean(dataset)) + '_cross_platform.orders').any?
  end

  def self.datasets(org)
    Dir.glob(File.join(customdir, org) + '/*.orders').collect{|f| clean(File.basename(f.sub(/.orders/,'')))}.uniq
  end

  def self.process_matrix(prefix, org)
    conditions = Dir.glob(prefix + '/*').collect{|f| File.basename(f)} - %w(values codes info description cross_platform)
    description = Open.read(File.join(prefix, 'description'))

    info = YAML.load(File.open(File.join(prefix, 'info')))
    r.CustomDS_process(prefix, false, conditions, description, info["two_channel"], !info["log2"])


    codes = Open.read(File.join(prefix,'codes')).collect{|l| l.chomp}
    cross_platform = ID.translate(org, codes) 

    if cross_platform.length > codes.length / 4
      Open.write(File.join(prefix,'cross_platform'),cross_platform.collect{|c| c || "NO MATCH"}.join("\n"))
      r.CustomDS_process(prefix, true, conditions, description, info["two_channel"], !info["log2"])
    end
  end

  def self.process(name)
    puts "Processing #{ name }"
    org = organism(name)
    prefix = File.join(customdir, org, name)

    CustomDS::process_matrix(prefix, org)
  end

end


if __FILE__ == $0
  p CustomDS::datasets('sgd')
  p CustomDS::path('HaploidData')
  p CustomDS::path('HaploidData_cross_platform')

  exit

  org = 'sgd'
  process = Dir.glob(File.join(CustomDS::customdir, org) + '/*').select{|f| File.directory? f}.collect{|f| File.basename(f)} - CustomDS.datasets('sgd')
  p process
  process.each{|d| CustomDS::process(d)}

end
