require 'MARQ'
require 'MARQ/MADB'
require 'MARQ/score'

module MARQ
  module Name
    def self.clean(name)
      name.sub(/_cross_platform/,'') unless name.nil?
    end

    def self.cross_platform(name)
      if name =~ /_cross_platform/
        name
      else
        name + "_cross_platform"
      end
    end

    def self.is_cross_platform?(name)
      ! name.match(/_cross_platform$/).nil?
    end

    def self.is_ratio?(name)
      ! name.match(/\[ratio\]$/).nil?
    end
  end

  module Platform
    def self.is_GEO?(platform)
      ! platform.match(/^GPL/).nil?
    end

    def self.is_cross_platform?(platform)
      MARQ::Name.is_cross_platform? platform
    end

    def self.path(platform)
      platform = MARQ::Name.clean(platform)
      if is_GEO? platform
        GEO.platform_path(platform)
      else
        CustomDS.platform_path(platform)
      end
    end

    def self.exists?(platform)
      path(platform) != nil
    end

    def self.has_cross_platform?(platform)
      File.exists? File.join(path(platform), 'cross_platform')
    end

    def self.datasets(platform)
      if is_GEO? platform
        GEO.platform_datasets(platform)
      else
        CustomDS.platform_datasets(platform)
      end
    end

    def self.codes(platform)
      platform = MARQ::Name.clean(platform)
      Open.read(File.join(path(platform), 'codes')).scan(/[^\s]+/)
    end

    def self.cross_platform(platform)
      platform = MARQ::Name.clean(platform)
      Open.read(File.join(path(platform), 'cross_platform')).scan(/[^\s]+/)
    end

    def self.organism(platform)
      platform = MARQ::Name.clean(platform)
      if is_GEO? platform
        GEO.platform_organism platform
      else
        CustomDS.platform_organism platform
      end
    end

    def self.process(platform)
      platform = MARQ::Name.clean(platform)
      if is_GEO? platform
        GEO.process_platform(platform)
      else
        CustomDS.process_platform(platform)
      end
    end

    def self.organism_platforms(organism)
      GEO.platforms.select {|platform|
        GEO::SOFT.GPL(platform)[:organism] == organism && MARQ::Platform.datasets(platform).any?
      } +
      CustomDS.organism_platforms(organism)
    end
  end


  module Dataset
    def self.is_GEO?(dataset)
      ! dataset.match(/^(?:GDS|GSE)/).nil?
    end

    def self.path(platform)
      if is_GEO? platform
        GEO.dataset_path(platform)
      else
        CustomDS.dataset_path(platform)
      end
    end

    def self.exists?(dataset)
      path = path(dataset)
      if path.nil?
        return false
      else
        return File.exists?(path + '.orders')
      end
    end

    def self.broken?(dataset)
      path = path(dataset)

      return false if path.nil?
      
      if File.exists?(path + '.skip')
        return true
      else
        return false
      end
    end

    def self.is_cross_platform?(dataset)
      MARQ::Name.is_cross_platform? dataset
    end

    def self.has_cross_platform?(dataset)
      File.exists?(path(dataset) + '_cross_platform.orders')
    end

    def self.info(name)
      begin
        title, description = Open.read(path(name) + '.description').split(/\n--\n/).values_at(0,1)
        {:title => title.strip, :description => description.strip}
      rescue Exception
        puts $!.message
        {:title => "" , :description => "" }
      end
    end

    def self.platform(dataset)
      if is_GEO? dataset
        GEO.dataset_platform(dataset)
      else
        CustomDS.dataset_platform(dataset)
      end
    end

    def self.organism(dataset)
      MARQ::Platform.organism(platform(dataset))
    end

    def self.process(dataset, platform = nil)
      if is_GEO? dataset
        GEO.process_dataset(dataset, platform)
      else
        CustomDS.process_dataset(dataset, platform)
      end
    end

    def self.read_file(dataset, ext)
      Open.read(path(dataset) + '.' + ext)
    end

    def self.read_values(dataset, file, integer = false)
      result = {}

      experiments = experiments(dataset)
      experiments.each{|experiment| result[experiment] = [] }
      read_file(dataset, file).split(/\n/).each do |line|
        values = line.chomp.split(/\t/)
        if integer
          values.each_with_index{|value, i| result[experiments[i]] << (value == 'NA' ? nil : value.to_i)}
        else
          values.each_with_index{|value, i| result[experiments[i]] << (value == 'NA' ? nil : value.to_f)}
        end
      end

      result
    end

    def self.read_values_t(dataset, file)
      experiments = experiments(dataset).reject{|experiment| MARQ::Name.is_ratio? experiment }
     
      return {} if experiments.empty?

      result = {}

      experiments.each{|experiment| result[experiment] = [] }

      read_file(dataset, file).split(/\n/).each do |line|
        values = line.chomp.split(/\t/)
        values.each_with_index{|value, i| result[experiments[i]] << (value == 'NA' ? nil : value.to_f) }
      end

      result
    end


    def self.experiments(dataset)
      read_file(dataset, 'experiments').split(/\n/).collect{|exp| exp.strip }
    end

    def self.codes(dataset)
      read_file(dataset, 'codes').split(/\n/)
    end

    def self.orders(dataset)
      read_values(dataset, 'orders', true)
    end

    def self.logratios(dataset)
      read_values(dataset, 'logratios')
    end

    def self.pvalues(dataset)
      read_values_t(dataset, 'pvalues')
    end

    def self.t(dataset)
      read_values_t(dataset, 't')
    end

    def self.codes_for(dataset, type, experiment)
      codes  = codes(dataset)
      values = send(type, dataset)[experiment]
      Hash[*codes.zip(values).reject{|p| p.last.nil? }.flatten]
    end

  end

  module RankQuery
    NULL_SIZE = 10000

    def self.dataset_scores(dataset, up, down)
      Score.scores_up_down(dataset, up, down)
    end

    def self.platform_scores(platform, up, down)
      scores = {}
      MARQ::Platform.datasets(platform).each do |dataset|
        dataset = MARQ::Name.cross_platform dataset if MARQ::Name.is_cross_platform?(platform)
        scores.merge!(dataset_scores(dataset, up, down))
      end

      scores
    end

    def self.organism_scores(organism, up, down)
      scores = {}
      MARQ::Platform.organism_platforms(organism).each do |platform|
        scores.merge!(platform_scores(MARQ::Name.cross_platform(platform), up, down))
      end

      scores
    end

    def self.add_pvalues(scores, up_size, down_size)
      null_scores = Score.null_scores(up_size, down_size, NULL_SIZE)
      Score.add_pvalues(scores, null_scores)
    end

  end
end
