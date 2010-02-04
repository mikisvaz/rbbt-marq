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
      result = {}

      experiments = experiments(dataset).select{|experiment| experiment !~ /\[ratio\]$/ }
      return {} if experiments.empty?
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

  end

  module RankQuery
    def self.complete_positions(positions, matched, genes)
      matched = matched.collect{|gene| gene.strip.downcase}
      genes   = genes.collect{|gene| gene.strip.downcase}

      pos = Hash[*matched.zip(positions).flatten]

      complete = genes.collect{|gene|
        if matched.include? gene
          pos[gene] || "MISSING"
        else
          "NOT IN PLATFORM"
        end
      }
      complete
    end


    def self.position_scores(up, down, positions_up, positions_down, platform_entries, matched_up, matched_down, missing_up, missing_down)
      scores = []
      positions_up.keys.each do |experiment|
        score = Score.score_up_down(positions_up[experiment], positions_down[experiment], platform_entries, missing_up, missing_down)
        score[:total_entries]    = platform_entries
        score[:positions_up]     = complete_positions(positions_up[experiment], matched_up, up) if up.any?
        score[:positions_down]   = complete_positions(positions_down[experiment], matched_down, down) if down.any?
        scores << score
      end

      pvalues = Score.pvalues(scores.collect{|s| s[:score]}, up.length, down.length, platform_entries)

      results = {}
      positions_up.keys.each_with_index{|experiment,i|
        results[experiment] = scores[i].merge(:pvalue => pvalues[i])
      }

      results
    end

    def self.dataset_scores(dataset, up, down)
      positions_up, matched_up, platform_entries     = MADB.dataset_positions(dataset, up)
      missing_up = positions_up.length - matched_up.length

      positions_down, matched_down                   = MADB.dataset_positions(dataset, down)
      missing_down = positions_down.length - matched_down.length

      position_scores(up, down, positions_up, positions_down, platform_entries, matched_up, matched_down, missing_up, missing_down)
    end

    def self.platform_scores(platform, up, down)
      positions_up, matched_up, platform_entries     = MADB.platform_positions(platform, up)
      missing_up = up.length - matched_up.length


      positions_down, matched_down                   = MADB.platform_positions(platform, down)
      missing_down = down.length - matched_down.length

      position_scores(up, down, positions_up, positions_down, platform_entries, matched_up, matched_down, missing_up, missing_down)
    end

    def self.organism_scores(organism, up, down)
      platforms = MARQ::Platform.organism_platforms(organism).
        select  {|p| MARQ::Platform.has_cross_platform? p }
        collect {|p| MARQ::name.cross_platform p }

      total_scores = {}
      platforms.each do |platform|
        scores = platform_scores(platform, up, down)
        total_scores.merge!(scores) 
      end

      total_scores
    end

  end
end

if __FILE__ == $0
  p MARQ::Dataset.platform 'GDS2791_cross_platform'
  p MARQ::Platform.organism 'GPL96'
  exit
  #puts MARQ::organism_platforms('human')
  #puts MARQ.platform_organism("HaploidData")
  #puts MARQ::platform_scores_up_down("HaploidData",%w( YMR261c YDL140c YIL122w YPL093w YHR211w YDL142c YHR106w YOR103c YDR233c YLR181c),%w()).keys

  up = %w(

  51228_at    215046_at   205009_at   204915_s_at 202707_at  
  208265_at   210618_at   201185_at   206650_at   200719_at  
  215661_at   202071_at   214408_s_at 215092_s_at 206168_at  
  212686_at   214162_at   221008_s_at 217709_at   210957_s_at

  )


  require 'MARQ/ID'
  require 'pp'
  genes = ID.translate('human',up).compact

  #pp up.zip(genes)
  #genes = Open.read("/home/miki/git/MARQ/test/GDS1375_malignant_vs_normal_down.genes").collect{|l| l.chomp.strip}
  positions =  MARQ::GEORQ.dataset_positions('GDS1231_cross_platform',genes)
  pp positions


  #MARQ::GEORQ.platform_scores_up_down('GPL96_cross_platform',genes,[]).each{|ex, r|
  #  puts ex
  #  puts r[:pvalue]
  #}

  #Score.draw_hits(positions["disease.state: malignant melanoma <=> normal"], MADB::GEORQ.experiment_entries('GPL96','GDS1375: disease.state: malignant melanoma <=> normal') , '/tmp/foo.png',:size => 1000)



end


