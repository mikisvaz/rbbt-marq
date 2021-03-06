require 'rbbt/util/open'
require 'MARQ'
require 'MARQ/ID'

module CustomDS
  @@r = nil

  DATA_DIR = File.join(MARQ.datadir,'CustomDS')


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

  def self.organism(dataset)
    path(dataset).match(/#{ DATA_DIR }\/(.*?)\//)[1]
  end

  def self.datasets(org)
    Dir.glob(File.join(DATA_DIR, org) + '/*.orders').collect{|f| clean(File.basename(f.sub(/.orders/,'')))}.uniq
  end

  def self.process_dataset(dataset, platform)
    org    = platform_organism(platform)
    platform_path = platform_path(platform)
    prefix = File.join(DATA_DIR, org, MARQ::Name.clean(dataset))
    
    conditions  = Dir.glob(File.join(platform_path, '*')).collect{|f| File.basename(f)} - %w(values codes info description cross_platform)
    description = Open.read(File.join(platform_path, 'description'))
    info        = YAML.load(File.open(File.join(platform_path, 'info')))

    if MARQ::Dataset.is_cross_platform?(dataset)
      r.CustomDS_process(prefix, true, conditions, description, info["two_channel"], !info["log2"])
    else
      r.CustomDS_process(prefix, false, conditions, description, info["two_channel"], !info["log2"])
    end

 
  end

  def self.organisms
    Dir.glob(File.join(DATA_DIR, '*')).
      select {|path| File.directory? path}.
      collect {|path| File.basename path}
  end

  def self.dataset_path(dataset)

    organisms.each do |organism|
      if File.exists?(File.join(DATA_DIR, organism, dataset + '.orders')) || File.exists?(File.join(DATA_DIR, organism, dataset + '.skip'))
        return File.join(DATA_DIR, organism, dataset)
      end
    end

    return nil
  end

  def self.platform_path(platform)
    Dir.glob(File.join(DATA_DIR, '*', platform)).first
  end

  def self.platform_datasets(platform)
    MARQ::Name.clean(platform)
  end

  def self.platform_organism(platform)
    path = platform_path(platform)
    return nil if path.nil?
    path.match(/#{DATA_DIR}\/(.*)\/#{ platform }$/)
    return $1
  end

  def self.dataset_organism(dataset)
    platform_organism(dataset)
  end

  def self.dataset_platform(dataset)
    dataset
  end

  def self.organism_platforms(organism)
    Dir.glob(File.join(DATA_DIR, organism, '*', 'codes')).
      collect {|path| File.basename(File.dirname(path))}.uniq
  end

  def self.process_platform(platform)
    prefix = platform_path(platform)
    org    = platform_organism(platform)

    codes = Open.read(File.join(prefix,'codes')).collect{|l| l.chomp}
    cross_platform = ID.translate(org, codes) 

    if cross_platform.length > codes.length / 4
      Open.write(File.join(prefix,'cross_platform'),cross_platform.collect{|c| c || "NO MATCH"}.join("\n"))
    end
  end
end


if __FILE__ == $0
  p CustomDS::datasets('sgd')
  p CustomDS::dataset_path('HaploidData')
  p CustomDS::dataset_path('HaploidData_cross_platform')

  exit

  org = 'sgd'
  process = Dir.glob(File.join(CustomDS::DATA_DIR, org) + '/*').select{|f| File.directory? f}.collect{|f| File.basename(f)} - CustomDS.datasets('sgd')
  p process
  process.each{|d| CustomDS::process(d)}

end
