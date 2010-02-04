require 'MARQ'
require 'MARQ/GEO'
require 'MARQ/CustomDS'
require 'MARQ/main'

module MADB

  # {{{ Saving Positions
  
  # Save the actual data, cross_platform or not
  def self.save_dataset_instance(dataset)

    # Get info
    codes        = MARQ::Dataset.codes(dataset);
    experiments  = MARQ::Dataset.experiments(dataset);
    orders       = MARQ::Dataset.orders(dataset).values.transpose;

    # Save codes and experiments
    DBcache.save(dataset + '_codes', codes)
    DBcache.save(dataset + '_experiments', experiments)

    # Asign orders to codes
    data = {}
    codes.each_with_index{|code,i|
      data[code] = orders[i]
    }

    # Save orders
    case 
    when codes.length < 65535
      type = "SMALLINT UNSIGNED"
    when codes.length < 16777215
      type = "MEDIUMINT UNSIGNED"
    else
      type = "INT UNSIGNED"
    end
    DBcache.save(dataset, data, [type] * orders.first.length)

  end

  # Save dataset, all instances, cross_platform if available.
  def self.save_dataset(dataset)
    save_dataset_instance(dataset)
    save_dataset_instance(MARQ::Name.cross_platform(dataset)) if MARQ::Dataset.has_cross_platform?(dataset)
    nil
  end

  def self.save_platform_instance(platform)
    DBcache.save(platform + '_codes', 
                 MARQ::Platform.is_cross_platform?(platform) ? 
                   MARQ::Platform.cross_platform(platform) : 
                   MARQ::Platform.codes(platform))
  end
  
  def self.save_platform(platform)
    datasets = MARQ::Platform.datasets(platform).sort
    return if datasets.empty?

    save_platform_instance(platform)
    save_platform_instance(MARQ::Name.cross_platform(platform)) if MARQ::Platform.has_cross_platform?(platform)
    
    datasets.sort.each do |dataset|
      save_dataset(dataset)
    end
  end

  # {{{ Loading Positions

  def self.platform_entries(platform)
    DBcache.num_rows(platform + '_codes')
  end


  def self.load_positions(dataset, genes, platform_entries)
    gene_positions = DBcache.load(dataset, genes)
    data = {}
    matched = gene_positions.keys 

    # Get experiments
    experiments = DBcache.load(dataset + '_experiments').sort{|a,b| 
      a[0].to_i <=> b[0].to_i
    }.collect{|p| 
      MARQ::Name.clean(dataset) + ": " + p[1].first
    }

    # Get scale factors (to account for genes missing in the dataset)
    scale = (0..experiments.length - 1).collect{|i| 
      rows = DBcache.num_rows(dataset, "C#{i}"); 
      if rows > 0 
        platform_entries / rows 
      else 
        nil 
      end 
    }

    # Get experiment positions and scale them
    experiment_x_gene = gene_positions.values.transpose
    experiments.each_with_index{|experiment, i|
      next if scale[i].nil? || experiment_x_gene[i].nil?
      values = experiment_x_gene[i].collect{|v| v.nil? ? nil : (v.to_f * scale[i]).to_i}
      data[experiment] = values
    }

    [data, matched, platform_entries]
  end
 
  def self.dataset_positions(dataset, genes)
    return [{},[],0] if genes.empty?

    genes = genes.collect{|gene| gene.downcase.strip}
    platform_entries = platform_entries(dataset)

    load_positions(dataset, genes, platform_entries)
  end


  def self.platform_positions(platform, genes)
    return [{},[],0] if genes.empty?

    genes = genes.collect {|gene| gene.downcase.strip }
    platform_entries = platform_entries(platform)

    cross_platform = MARQ::Platform.is_cross_platform? platform
    datasets = MARQ::Platform.datasets(platform).sort

    total_data = {}
    total_matched = []

    datasets.each do |dataset|
      dataset = MARQ::Name.cross_platform dataset if cross_platform
      data, matched = load_positions(dataset, genes, platform_entries)
      total_data = total_data.merge(data)
      total_matched +=  matched
    end
    total_matched.uniq!

    [total_data, total_matched, platform_entries]
  end

end

if __FILE__ == $0
  #CustomDS::datasets('sgd').each{|d| MADB::CustomDS::save(d)}

  require 'pp'
  pp MADB::dataset_positions('GDS113', %w(2778))[0]
  pp MADB::platform_positions('GPL54', %w(2778))[0]
  #p MADB::CustomDS::positions("HaploidData",%w( YMR261c YDL140c YIL122w YPL093w YHR211w YDL142c YHR106w YOR103c YDR233c YLR181c yomeman))
  #p MADB::CustomDS::positions("HaploidData_cross_platform",%w( S000002685 S000001149 S000003068 S000003153 S000003355 S000000127 S000004444 S000004875 S000001702 S000005843 S000000862))
end
