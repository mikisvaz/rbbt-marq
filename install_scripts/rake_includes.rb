require 'progress-monitor'

$expr_threshold ||= (ENV['threshold'] || 0.05).to_f
$folds          ||= (ENV['folds'] || 2.5).to_f
$nth_genes      ||= (ENV['nth_genes'] || 100).to_i

$force       = [$force, ENV['force'], false].compact.first.to_s == 'true'  
$tranlations = [$tranlations, ENV['translations'], false].compact.first.to_s == 'true'  
$update_db   = [$update_db, ENV['update_db'], false].compact.first.to_s == 'true'  
$skip_db     = [$skip_db, ENV['skip_db'], false].compact.first.to_s == 'true'  
$fdr         = [$fdr, ENV['fdr'], true].compact.first.to_s == 'true'  
$do_folds    = [$do_folds, ENV['do_folds'], true].compact.first.to_s == 'true'  


$changes = false
module GEO::Process::R
  class << self
    alias_method :GDS_old, :GDS
    def GDS(*args)
      $changes = true
      GDS_old(*args)
    end
  
    alias_method :GSE_old, :GSE
    def GSE(*args)
      $changes = true
      GSE_old(*args)
    end
  end
end

module CustomDS
  class << self
    alias_method :process_dataset_old, :process_dataset
    def process_dataset(*args)
      $changes = true
      process_dataset_old(*args)
    end
  end
end

desc "Analyze datasets"
task 'data' do 
  platforms_to_save = []

  platforms = process_list

  Progress.monitor("Processing #{platforms.keys.length} platforms") if platforms.keys.length > 1
  platforms.each{|platform, datasets|
    begin
      # Prepare the platform
      MARQ::Platform.process(platform)
    rescue
      puts "Error processing platform #{platform}"
      puts $!.message
      puts "\n" * 3
      #puts $!.backtrace.join("\n")
      next
    end
    
    next if $tranlations

    $changes = false
    # Process all datasets
    
    Progress.monitor("Processing #{datasets.length} datasets", :announcement => Proc.new{|d| "Dataset #{ d }" }) if datasets.length > 1
    datasets.each{|dataset|
      begin
        already_processed = MARQ::Dataset.exists?(dataset) || MARQ::Dataset.broken?(dataset)
        next if already_processed && ! $force

        MARQ::Dataset.process(dataset, platform)
        MARQ::Dataset.process(MARQ::Name.cross_platform(dataset), platform) if MARQ::Platform.has_cross_platform?(platform)
      rescue
        puts "Error processing dataset #{ dataset }"
        puts $!.message
        puts "\n" * 3
        #puts $!.backtrace.join("\n")
      end
    }
    
    # Mark the platform for saving in DB
    platforms_to_save << platform if $changes || $update_db
  }

  Progress.monitor("Saving #{platforms_to_save.length} platforms in DB") if platforms_to_save.length > 1
  platforms_to_save.each{|platform|
    begin
      MADB.save_platform(platform) 
    rescue
      puts "Error saving platform #{ platform }"
      puts $!.message
      puts "\n" * 3
      #puts $!.backtrace.join("\n")
    end
  }
end

def annotations(name, cross_platform = false, &block)
  platforms = process_list

  Progress.monitor("Processing #{platforms.keys.length} platforms for #{ name } annotations")
  platforms.each do |platform, datasets|
    next if ! MARQ::Platform.exists? platform
    Progress.monitor("Processing #{datasets.length} datasets", :announcement => Proc.new{|d| "Dataset #{ d }" }) if datasets.length > 1
    datasets.each do |dataset|
      begin
        next if File.exist?(File.join("annotations", name, dataset)) && ! $force
        next if MARQ::Dataset.path(dataset).nil?

        FileUtils.mkdir_p File.join("annotations", name)
        filename = File.join("annotations", name, dataset)
        dataset += MARQ::Name.cross_platform(dataset) if cross_platform && MARQ::Platform::has_cross_platform?(platform)
        next if ! MARQ::Dataset.exists?(dataset)
        terms = block.call(dataset)
        Open.write(filename, terms.to_yaml)
      rescue
        puts "Error processing dataset #{ dataset }"
        puts $!.message
        puts $!.backtrace.join("\n")
      end
    end
  end
end


task 'annotate_Words' do
  require 'MARQ/annotations'
  require 'rbbt/bow/bow'
  annotations('Words') do |dataset|
    terms = {}
    description = Open.read(MARQ::Dataset.path(dataset) + '.description')
    terms[:dataset] = [dataset] +  description.words.uniq
    Open.read(MARQ::Dataset.path(dataset) + '.experiments').collect{|name|
      name = name.strip
      terms[name] = name.sub(/.*?: /,'').sub(/\[ratio\]/,'').words.uniq
    }
    terms
  end
end


task 'annotate_UMLS' do
  require 'MARQ/annotations'
  require 'rbbt/util/misc'
  annotations('UMLS') do |dataset|
    terms = {}
    description = Open.read(MARQ::Dataset.path(dataset) + '.description')
    terms[:dataset] = Annotations::UMLS::OBA(description).uniq
    Open.read(MARQ::Dataset.path(dataset) + '.experiments').collect{|name|
      name = name.strip
      terms[name] = Annotations::UMLS::OBA(name.sub(/.*?: /,'').sub(/\[ratio\]/,'')).uniq
    }
    terms
  end
end


task 'annotate_Polysearch' do
  require 'MARQ/annotations'
  require 'rbbt/util/misc'
  require 'rbbt/sources/polysearch'
  annotations('Polysearch') do |dataset|
    terms = {}
    description = Open.read(MARQ::Dataset.path(dataset) + '.description')
    terms[:dataset] = Polysearch::match(description).values.flatten.sort.collect{|n| n.gsub(/\s+/,' ').downcase}.uniq
    Open.read(MARQ::Dataset.path(dataset) + '.experiments').collect{|name|
      name = name.strip
      terms[name] = Polysearch::match(name.sub(/.*?: /,'').sub(/\[ratio\]/,'')).values.flatten.sort.collect{|n| n.gsub(/\s+/,' ').downcase}.uniq
    }
    terms
  end

end

def goterms(org, list, slim, threshold)
  return [] if list.empty?
  results = Annotations::Genes::Genecodis::Local.analysis(org, list, slim)
  return [] if results.nil?
  results.
    select{|info| info[:s].to_i  > 2 }.
    select{|info| info[:hyp_c].to_f < threshold }.
    collect{|info| info[:items]}.collect{|id| GO::id2name(id)}
end

task 'annotate_GO' do
  require 'MARQ/annotations'
  require 'rbbt/sources/go'
  options = { :cut_off => $expr_threshold, :fdr => $fdr, :folds => $folds, :do_folds => $do_folds, :nth_genes => $nth_genes}
  annotations('GO_up', true) do |dataset|
    org = MARQ::Dataset.organism(dataset)
    genes = Annotations::Genes.get_genes(dataset, options)

    up = {}
    genes[:up] ||= []
    genes[:up].collect{|experiment,list|
      up[experiment] =  goterms(org, list, false, $expr_threshold)
    }
    up
  end

  annotations('GO_down', true) do |dataset|
    org = MARQ::Dataset.organism(dataset)
    genes = Annotations::Genes.get_genes(dataset, options)

    down = {}
    genes[:down] ||= []
    genes[:down].collect{|experiment,list|
      down[experiment] = goterms(org, list, false, $expr_threshold)
    }
    down
  end 

  annotations('GOSlim_up', true) do |dataset|
    org = MARQ::Dataset.organism(dataset)
    genes = Annotations::Genes.get_genes(dataset, options)

    up = {}
    genes[:up] ||= []
    genes[:up].collect{|experiment,list|
      up[experiment] = goterms(org, list, true, $expr_threshold)
    }
    up
  end

  annotations('GOSlim_down', true) do |dataset|
    org = MARQ::Dataset.organism(dataset)
    genes = Annotations::Genes.get_genes(dataset, options)

    down = {}
    genes[:down] ||= []
    genes[:down].collect{|experiment,list|
      down[experiment] =  goterms(org, list, true, $expr_threshold)
    }
    down
  end
end

task 'annotate_SENT' do
  require 'MARQ/annotations'
  options = { :cut_off => $expr_threshold, :fdr => $fdr, :folds => $folds, :do_folds => $do_folds, :nth_genes => $nth_genes}
  annotations('SENT') do |dataset|
    org = MARQ::Dataset.organism(dataset)
    genes = Annotations::Genes.get_genes(dataset, options)
    terms = Annotations::Genes::SENT.terms(org, genes)
    terms
  end


end

task 'default' do
  Rake::Task['data'].invoke
  Rake::Task['annotate_Words'].invoke
  #Rake::Task['annotate_UMLS'].invoke
  #Rake::Task['annotate_Polysearch'].invoke
  Rake::Task['annotate_GO'].invoke
end
