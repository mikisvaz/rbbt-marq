require 'MARQ/util'
require 'MARQ/main'
module RankProduct
  def self.ranks(dataset, experiment , from_FC = false, invert = false)
    codes = MARQ::Dataset.codes(dataset)
    
    if from_FC 
      ratios    = MARQ::Dataset.logratios(dataset)[experiment.strip]
      sorted_genes = codes.zip(ratios).
        reject  {|p| p[1].nil? }.
        sort_by {|p| p[1] }.
        collect {|p| p[0] }
      sorted_genes.reverse! unless invert
      ranks  = Hash[*sorted_genes.zip((1..sorted_genes.length).to_a).flatten]
      (codes - sorted_genes).each {|gene| ranks[gene] = nil}
    else
      orders = MARQ::Dataset.orders(dataset)[experiment.strip]

      if invert
        num_genes = codes.length + 1
        orders.collect! {|pos| pos.nil? ? nil : num_genes - pos }
      end

      ranks = Hash[*codes.zip(orders).flatten]
    end
    ranks
  end

  def self.score(gene_ranks, signature_sizes)
    scores = {}
    log_sizes = signature_sizes.collect{|size| Math::log(size)}
    gene_ranks.each{|gene, positions|
      scores[gene] = positions.zip(log_sizes).
        collect{|p| Math::log(p[0]) - p[1]}.    # Take log and substract from size (normalize)
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

    # Gather gene ranks from signatures
    ranks = {}
    signatures.each{|signature|
      dataset, experiment = signature.match(/^([^\:]*): (.*)/).values_at(1,2)
      dataset = dataset + '_cross_platform' if cross_platform
      ranks[signature] = self.ranks(dataset, experiment, from_FC, invert.include?(signature))
    }

    # Invert the hash, from signature keys to gene keys
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

    # Remove incomplete genes
    gene_ranks.delete_if{|code, positions| positions.length != signatures.uniq.length}
    
    # Compute scores
    scores = score(gene_ranks, sizes)

    # Compute permutations
    num_permutations = 50000
    permutation_scores = permutations(sizes.length, num_permutations)
    permutation_scores = permutation_scores.sort


    # Compute p-values from permutations
    results = {}
    scores.each {|gene, score|
      pos = permutation_scores.count_smaller(score)
      results[gene] = [score, pos.to_f / num_permutations]
    }

    # Complete the information with pfp
    num_genes = results.length
    results.sort {|a,b|
      a[1][0] <=> b[1][0]
    }.each_with_index{|p,i|
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
