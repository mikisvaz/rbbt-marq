require 'MARQ'
module RankProduct

  def self.ranks(dataset, experiment , from_FC = false, invert = false)
    path = MARQ.dataset_path(dataset)
    codes = Open.read(path + '.codes').collect{|line| line.strip}

    experiments = Open.read(path + '.experiments').collect{|line| line.strip}
    field = experiments.index experiment.strip

    if from_FC
      ratios = Open.read(path + '.logratios').collect{|line| 
        value = line.strip.split("\t")[field]
        case
        when value == "NA"
          nil

        # To sort decreasingly we change sign by default
        when invert
          value.to_f
        else
          - value.to_f
        end
      }
      Hash[*codes.zip(ratios).sort{|a,b| b[1] <=> a[1]}.collect{|p| p[0]}.zip((1..codes.length).to_a).flatten]
    else
      orders = Open.read(path + '.orders').collect{|line| 
        value = line.strip.split("\t")[field]
        case
        when value == "NA"
          nil
        when invert
          codes.length - line.strip.split("\t")[field].to_i + 1
        else
          line.strip.split("\t")[field].to_i
        end

      }
      Hash[*codes.zip(orders).flatten]
    end
  end

  def self.score(gene_ranks, signature_sizes)
    scores = {}
    log_sizes = signature_sizes.collect{|size| Math::log(size)}
    gene_ranks.each{|gene, positions|
      scores[gene] = positions.zip(log_sizes).
        collect{|p| Math::log(p[0]) - p[1]}.
        inject(0){|acc, v| acc += v  }
    }
    scores
  end

  def self.permutations(num_signatures, num = 1000)
    scores = []
    num.times{
       value = 0
       num_signatures.times{|size_and_log| 
         value += Math::log(rand)
       } 
       scores << value
    }
    scores
  end

  def self.permutations_full(signature_sizes)
    gene_ranks = {}
    signature_sizes.each{|size|
      (1..size).to_a.shuffle.each_with_index{|gene, pos|
        gene_ranks[gene] ||= []
        gene_ranks[gene] << pos + 1
      }
    }
    gene_ranks.delete_if{|code, positions| positions.length != signature_sizes.length}

    scores = score(gene_ranks, signature_sizes)
    scores.values
  end

  def self.rankproduct(signatures, options = {})
    invert, from_FC, cross_platform = {
      :invert => [],
      :from_FC => false,
      :cross_platform => false,
    }.merge(options).values_at(:invert, :from_FC, :cross_platform)

    ranks = {}
    signatures.each{|signature|
      dataset, experiment = signature.match(/^([^\:]*): (.*)/).values_at(1,2)
      dataset = dataset + '_cross_platform' if cross_platform
      ranks[signature] = self.ranks(dataset, experiment, from_FC, invert.include?(signature))
    }

    gene_ranks = {}
    sizes = []
    ranks.each{|signature, orders|
      sizes << orders.length
      orders.each{|code, position|
        next if position.nil?
        gene_ranks[code] ||= []
        gene_ranks[code] << position
      }
    }

    gene_ranks.delete_if{|code, positions| positions.length != signatures.uniq.length}
    
    scores = score(gene_ranks, sizes)
    num_permutations = 50000

    permutation_scores = permutations(sizes.length, num_permutations)

    permutation_scores = permutation_scores.sort


    results = {}
    scores.each{|gene, score|
      pos = permutation_scores.count_smaller(score)
      results[gene] = [score, pos.to_f / num_permutations]
    }


    num_genes = results.length
    results.sort{|a,b|
      a[1][0] <=> b[1][0]
    }.each_with_index{|p,i|
      gene = p[0]
      info = p[1]
      pvalue = info[1]
      pfp  = pvalue * num_genes / (i + 1)
      info << pfp
    }
       
    results
  end

end

if __FILE__ == $0
  scores = RankProduct.rankproduct(
    [ 'GDS2952: disease.state: rheumatoid arthritis <=> normal',
      #'GDS2952: protocol: pre-treatment <=> untreated',
      'GDS2952: protocol: post-treatment <=> untreated'], :cross_platform => true, :from_FC => true)

  p scores.sort{|a,b| a[1][0] <=> b[1][0]}[1..10]
  #scores.sort{|a,b| a[1][0] <=>  b[1][0]}.reverse
end
