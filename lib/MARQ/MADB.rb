require 'MARQ'
require 'MARQ/GEO'
require 'MARQ/CustomDS'

module MADB
  module CustomDS

    def self.save(dataset)
      prefix = Object::CustomDS.path(dataset)

      codes = File.open(prefix  + '.codes').collect{|l| l.chomp.downcase}
      
      DBcache.save(dataset + '_codes', codes)

      experiments = File.open(prefix + '.experiments').collect{|l| l.chomp}
      orders = File.open(prefix + '.orders').collect{|l| values = l.chomp.split(/\t/).collect{|v| v == "NA" ? nil : v.to_i };}

      data = {}
      codes.each_with_index{|code,i|
        data[code.to_sym] = orders[i]
      }
      case 
      when codes.length < 65535
        type = "SMALLINT UNSIGNED"
      when codes.length < 16777215
        type = "MEDIUMIN UNSIGNED"
      else
        type = "INT UNSIGNED"
      end

      DBcache.save(dataset + '_experiments', experiments)
      DBcache.save(dataset, data, [type] * orders.first.length)

      return unless Object::CustomDS::has_cross_platform?(dataset)
      dataset = dataset + '_cross_platform'
      prefix = Object::CustomDS.path(dataset)

      codes = File.open(prefix  + '.codes').collect{|l| l.chomp.downcase}
      
      DBcache.save(dataset + '_codes', codes)

      experiments = File.open(prefix + '.experiments').collect{|l| l.chomp}
      orders = File.open(prefix + '.orders').collect{|l| values = l.chomp.split(/\t/).collect{|v| v == "NA" ? nil : v.to_i };}

      data = {}
      codes.each_with_index{|code,i|
        data[code.to_sym] = orders[i]
      }
      case 
      when codes.length < 65535
        type = "SMALLINT UNSIGNED"
      when codes.length < 16777215
        type = "MEDIUMIN UNSIGNED"
      else
        type = "INT UNSIGNED"
      end

      DBcache.save(dataset + '_experiments', experiments)
      DBcache.save(dataset, data, [type] * orders.first.length)
      nil
    end
   
    def self.positions(dataset, genes)
      return [{},[]] if genes.empty?
      genes = genes.collect{|gene| gene.downcase.strip}

      platform_entries = platform_entries(dataset + '_codes').to_f

      data = {}
      matched = []

      gene_positions = DBcache.load(dataset, genes)
      matched ||= gene_positions.keys

      experiments = DBcache.load(dataset + '_experiments').sort{|a,b| 
        a[0].to_i <=> b[0].to_i
      }.collect{|p| 
        Object::GEO::clean(dataset) + ": " + p[1].first
      }


      matched = (matched + gene_positions.keys).uniq
      scale = (0..experiments.length - 1).collect{|i| 
        rows = DBcache.num_rows(dataset, "C#{i}"); 
        if rows > 0 
          platform_entries / rows 
        else 
          nil 
        end 
      }

      gene_x_experiment = gene_positions.values

      experiment_x_gene = gene_x_experiment.transpose

      experiments.each_with_index{|experiment, i|
        next if scale[i].nil? || experiment_x_gene[i].nil?
        values = experiment_x_gene[i].collect{|v| v.nil? ? nil : (v.to_f * scale[i]).to_i}
        data[experiment] = values
      }
      
      [data, matched]
    end

    def self.platform_entries(platform)
      DBcache.num_rows(platform)
    end
  end





  module GEO

    def self.saveGPL(platform)
      datasets = Object::GEO.platform_datasets(platform).sort
      return if datasets.empty?


      codes = File.open(File.join(Object::GEO.platform_path(platform),'codes')).collect{|l| l.chomp.downcase}
      
      DBcache.save(platform, codes)

      datasets.sort.each{|dataset|
        path = Object::GEO.dataset_path(dataset)
        experiments = File.open(path + '.experiments').collect{|l| l.chomp}
        orders = File.open(path + '.orders').collect{|l| values = l.chomp.split(/\t/).collect{|v| v == "NA" ? nil : v.to_i };}

        data = {}
        codes.each_with_index{|code,i|
          data[code.to_sym] = orders[i]
        }
        case 
        when codes.length < 65535
          type = "SMALLINT UNSIGNED"
        when codes.length < 16777215
          type = "MEDIUMINT UNSIGNED"
        else
          type = "INT UNSIGNED"
        end

        DBcache.save(dataset + '_experiments', experiments)
        DBcache.save(dataset, data, [type] * orders.first.length)
      }


      return unless File.exist?(File.join(Object::GEO.platform_path(platform),'cross_platform'))
      codes = File.open(File.join(Object::GEO.platform_path(platform),'cross_platform')).collect{|l| l.chomp.downcase}

      DBcache.save(platform + '_cross_platform', codes)

      Progress.monitor("Saving #{ platform }")
      datasets.sort.each{|dataset|
        path = Object::GEO.dataset_path(dataset)
        next unless File.exists?(path + '_cross_platform.experiments')
        experiments = File.open(path + '_cross_platform.experiments').collect{|l| l.chomp}
        orders = File.open(path + '_cross_platform.orders').collect{|l| values = l.chomp.split(/\t/).collect{|v| v == "NA" ? nil : v.to_i };}

        data = {}
        codes.each_with_index{|code,i|
          data[code.to_sym] = orders[i]
        }

        case 
        when codes.length < 65535
          type = "SMALLINT UNSIGNED"
        when codes.length < 16777215
          type = "MEDIUMIN UNSIGNED"
        else
          type = "INT UNSIGNED"
        end


        DBcache.save(dataset + '_cross_platform_experiments', experiments)
        DBcache.save(dataset + '_cross_platform', data, [type] * orders.first.length)
      }
    end

    def self.positions(platform, genes)
      return [{},[]] if genes.empty?
      genes = genes.collect{|gene| gene.downcase.strip}

      datasets = Object::GEO.platform_datasets(platform).sort
      platform_entries = platform_entries(platform).to_f

      data = {}
      matched = nil

      datasets.each{|dataset|
        dataset += '_cross_platform' if Object::GEO::is_cross_platform?(platform)
        gene_positions = DBcache.load(dataset, genes)
        matched ||= gene_positions.keys
        
        experiments = DBcache.load(dataset + '_experiments').sort{|a,b| 
          a[0].to_i <=> b[0].to_i
        }.collect{|p| 
          Object::GEO::clean(dataset) + ": " + p[1].first
        }

        scale = (0..experiments.length - 1).collect{|i| 
          rows = DBcache.num_rows(dataset, "C#{i}"); 
          if rows > 0 
            platform_entries / rows 
          else 
            nil 
          end 
        }
        
        gene_x_experiment = gene_positions.values

        experiment_x_gene = gene_x_experiment.transpose

        experiments.each_with_index{|experiment, i|
          next if scale[i].nil? || experiment_x_gene[i].nil?
          values = experiment_x_gene[i].collect{|v| v.nil? ? nil : (v.to_f * scale[i]).to_i}
          data[experiment] = values
        }
      }

      [data, matched]
    end

    def self.platform_entries(platform)
      DBcache.num_rows(platform)
    end
  end

end

if __FILE__ == $0
  #CustomDS::datasets('sgd').each{|d| MADB::CustomDS::save(d)}

  require 'pp'
  pp MADB::GEO::positions('GPL91_cross_platform', %w(2778))[0].select{|k,v| k =~ /GDS989/}.sort
  #p MADB::CustomDS::positions("HaploidData",%w( YMR261c YDL140c YIL122w YPL093w YHR211w YDL142c YHR106w YOR103c YDR233c YLR181c yomeman))
  #p MADB::CustomDS::positions("HaploidData_cross_platform",%w( S000002685 S000001149 S000003068 S000003153 S000003355 S000000127 S000004444 S000004875 S000001702 S000005843 S000000862))
end
