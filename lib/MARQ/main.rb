require 'MARQ'
require 'MARQ/MADB'
require 'MARQ/score'

module MARQ
  def self.platform_type(platform)
    if platform.match(/GPL\d+|GDS\d+|GSE\d+/)
      return :GEO
    else
      return :CustomDS
    end
  end

  def self.dataset_path(dataset)
    if platform_type(dataset) == :GEO
      GEO::dataset_path(dataset)
    else
      CustomDS::path(dataset)
    end
  end

  def self.platform_organism(platform)
    if platform_type(platform) == :GEO
      if platform.match(/^GPL/)
        GEO::SOFT::GPL(platform)[:organism]
      else
        GEO::SOFT::GPL(GEO::dataset_platform(platform))[:organism]
      end
    else
      CustomDS::organism(platform)
    end
  end

  def self.organism_platforms(org)
    GEO::organism_platforms(org) + CustomDS::datasets(org)
  end

  def self.has_cross_platform?(dataset, platform=nil)
    if platform_type(platform) == :GEO
      GEO::has_cross_platform?(dataset, platform)
    else
      CustomDS::has_cross_platform?(platform)
    end
  end

  def self.is_cross_platform?(dataset, platform=nil)
    if platform_type(platform) == :GEO
      GEO::is_cross_platform?(dataset, platform)
    else
      CustomDS::is_cross_platform?(platform)
    end
  end

  def self.complete_positions(positions, matched, genes)

    pos = Hash[*matched.zip(positions).flatten]
    complete = genes.collect{|gene|
      gene = gene.downcase.strip
      if matched.include? gene
        pos[gene] || "MISSING"
      else
          "NOT IN PLATFORM"
      end
    }

  end

  module Dataset

    def self.exists?(dataset)
      MARQ::dataset_path(dataset) != nil
    end

    def self.read_file(dataset, ext)
      Open.read(MARQ::dataset_path(dataset) + '.' + ext)
    end

    def self.read_values(dataset, file)
      result = {}

      experiments = experiments(dataset)
      experiments.each{|experiment| result[experiment] = [] }
      read_file(dataset, file).split(/\n/).each do |line|
        values = line.chomp.split(/\t/)
        values.each_with_index{|value, i| result[experiments[i]] << value.to_f }
      end

      result
    end

    def self.platform_codes(platform)
      if MARQ::is_cross_platform? platform
        file = 'cross_platform'
      else
        file = 'codes'
      end

      Open.read(File.join(MARQ::platform_path(platform), file))
    end

    def self.experiments(dataset)
      read_file(dataset, 'experiments').split(/\n/)
    end

    def self.codes(dataset)
      read_file(dataset, 'codes').split(/\n/)
    end

    def self.orders(dataset)
      read_values(dataset, 'orders')
    end

    def self.logratios(dataset)
      read_values(dataset, 'logratios')
    end
    
    def self.pvalues(dataset)
      read_values(dataset, 'pvalues')
    end
  end

  def self.platform_scores_up_down(platform, up, down)
    if platform_type(platform) == :GEO
      GEORQ.platform_scores_up_down(platform, up, down)
    else
      CustomDSRQ.scores_up_down(platform, up, down)
    end
  end

  module CustomDSRQ

    def self.scores_up_down(platform, up, down)
      matched_up = DBcache.matches(platform + '_codes', up).length
      positions_up, matched_up = *MADB::CustomDS.positions(platform, up)
      missing_up = up.length - matched_up.length



      matched_down = DBcache.matches(platform + '_codes', down).length
      positions_down, matched_down = *MADB::CustomDS.positions(platform, down)
      missing_down = down.length - matched_down.length

      experiments = positions_up.keys
      platform_entries = MADB::CustomDS.platform_entries(platform)

      scores = experiments.collect{|experiment|
        score = Score.score_up_down(positions_up[experiment], positions_down[experiment], platform_entries, missing_up, missing_down)
        score[:total_entries]    = platform_entries
        score[:positions_up]     = MARQ.complete_positions(positions_up[experiment], matched_up, up) if positions_up.any?
        score[:positions_down]   = MARQ.complete_positions(positions_down[experiment], matched_down, down) if positions_down.any?
        score
      }

      pvalues = Score.pvalues(scores.collect{|s| s[:score]}, up.length, down.length, platform_entries)

      results = {}
      experiments.each_with_index{|experiment,i|
        results[experiment] = scores[i].merge(:pvalue => pvalues[i])
      }


      results
    end
  end

  module GEORQ

    def self.platform_scores_up_down(platform, up, down)
      platform_entries = MADB::GEO.platform_entries(platform)

      positions_up, matched_up = *MADB::GEO.positions(platform, up)
      missing_up = up.length - matched_up.length

      positions_down, matched_down = *MADB::GEO.positions(platform, down)
      missing_down = down.length - matched_down.length

      return {} if positions_up.keys.empty? && positions_down.keys.empty?
      if positions_up.keys.any?
        experiments = positions_up.keys
      else 
        experiments = positions_down.keys
      end

      scores = experiments.collect{|experiment|
        score = Score.score_up_down(positions_up[experiment], positions_down[experiment], platform_entries, missing_up, missing_down)
        score[:total_entries]    = platform_entries
        score[:positions_up]     = MARQ.complete_positions(positions_up[experiment], matched_up, up) if positions_up.any?
        score[:positions_down]   = MARQ.complete_positions(positions_down[experiment], matched_down, down) if positions_down.any?
        score
      }

      pvalues = Score.pvalues(scores.collect{|s| s[:score]}, up.length, down.length, platform_entries)

      results = {}
      experiments.each_with_index{|experiment,i|
        results[experiment] = scores[i]
        results[experiment][:pvalue] = pvalues[i]
      }


      results
    end


    def self.draw_hits(platform, genes, directory)
      positions, matched = MADB::GEO.positions(platform, genes).values_at(0,1)

      platform_entries = MADB::GEO.platform_entries(platform)

      positions.each{|experiment, positions|
        Score.draw_hits(positions.compact, platform_entries, File.join(directory, experiment.hash_code + '.png'), {:size => 1000, :bg_color => :green})
      }
    end


    def self.dataset_positions(dataset, genes)
      if dataset =~ /_cross_platform/
        dataset.sub!(/_cross_platform/,'')
        platform = GEO.dataset_platform(dataset) + '_cross_platform'
      else
        platform = GEO.dataset_platform(dataset)
      end


      positions = {}
      MADB::GEO.positions(platform, genes).first.select{|k,v| k =~ /^#{ dataset }/}.each{|p|
        positions[p[0].sub(/^#{ dataset }: /,'')] = p[1]
      }
      positions
    end

    def self.dataset_scores(dataset, genes)
      total = DBcache::num_rows(platform)
      positions = dataset_positions(dataset, genes)

      scores = {}
      positions.each{|experiment, positions|
        scores[experiment] = Score.score(positions, total, genes.length - positions.length)
      }
      scores
    end

  end
end

if __FILE__ == $0
  p MARQ::Dataset.orders('GDS2419')
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


