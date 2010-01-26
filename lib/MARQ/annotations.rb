require 'inline'
require 'net/http'
require 'uri'
require 'MARQ'
require 'rbbt/bow/dictionary'
require 'MARQ/fdr'
require 'digest/md5'
require 'base64'
require 'rbbt/util/open'

module Annotations
  class << self
    inline do |builder|
      builder.c_raw <<-EOC

/**
 * Compute log(k!)
 * @param k The value k.
 * @return The result.
 */
double lFactorial(double k)
{
	double r = 0;
	int i;
	for(i=2 ; i<=(int)k ; i++)
	{
		r = r + (double)(log((double)i));
	}
	return r;
}



/**
 * Compute the log(binom(n,k))
 * @param n The number of possible items.
 * @param k The number of selected items.
 * @return The result.
 */
double lBinom(double n, double k)
{
	long i;
	double r = 0;
	
	if(n > n-k){
		k = n-k;
	}

	for(i = (long)n ; i> (n-k) ; i--)
	{
		r = r + log((double)i);
	}
		
	r = r - lFactorial(k);
	
	return r;
}
  EOC
  
      builder.c <<-EOC
/**
*  * Compute the Hypergeometric accumulated value.
*  * @param total => total size
*  * @param support => total support
*  * @param list => selected list size,
*  * @param found => support
*  * @return The result
*  */
double hypergeometric(double total, double support, double list, double found)
{
	double other = total - support;

	double top = list;
	if(support < list){
		top = support;
	}

	double log_n_choose_k = lBinom(total,list);

	double lfoo = lBinom(support,top) + lBinom(other, list-top);
	
	double sum = 0;
  int i;
	for (i = (int)top; i >= found; i-- )
	{
		sum = sum + exp(lfoo - log_n_choose_k);
		if ( i > found)
		{
			lfoo = lfoo + log(i / (support - i+1)) +  log( (other - list + i) / (list-i+1)  );
		}
	}
	return sum;
}
      EOC
    end
  end

  def self.exp2gds(experiment)
    experiment =~ /(.*?):/
    $1
  end

  def self.compare(a,b)
    case 
    when a[:pvalue] < b[:pvalue]
      -1 
    when a[:pvalue] > b[:pvalue]
      1 
    when a[:pvalue] == b[:pvalue]
      b[:score].abs <=> a[:score].abs
    end
  end

  RANK_SIZE_BINS = %w(1 2 3 4 5 7 10 15 20 30 40 50 65 80 100 125 150 175 200 250 300 350 400 450 500 600 700 800 900 1000 1500 2000 2500 3000)

  def self.enrichment_rank(annotations, ranks, options = {})
    dict_options = {
      :dict_options => {:low => 0, :hi => 0.5, :limit => 100000}
    }.merge(options)[:dict_options]
    positions = {}
    found_datasets = []

    dict = Dictionary::TF_IDF.new
    ranks.each_with_index{|experiment, rank|
      info = annotations[experiment]

      dataset_terms   = info[:dataset]
      signature_terms = info[:signature]
      
      dataset = exp2gds experiment

      terms = signature_terms
      terms += dataset_terms

      term_count = {}
      terms.each{|term|
        term_count[term] ||= 0
        term_count[term] += 1
      }
      dict.add(term_count)
    }
 
    best = dict.best(dict_options).keys

    found_datasets = []
    ranks.each_with_index{|experiment, rank|
      info = annotations[experiment]

      dataset_terms   = info[:dataset]
      signature_terms = info[:signature]

      dataset = exp2gds experiment

      terms = signature_terms
      
      if ! found_datasets.include? dataset
        terms += dataset_terms
        found_datasets << dataset
      end

      terms.uniq.each{|term|
        next if not best.include? term
        positions[term] ||= []
        positions[term] <<  rank
      }
    }
 
    scores = []


    sizes = {}
    RANK_SIZE_BINS.each{|size| sizes[size.to_i] = []}


    # For each term compute the rank score. Also, place it in the closest size
    # bin for the permutations.
    best.each_with_index{|term, pos|
      if positions[term]
        list = positions[term]

        # place it on the size bin 
        found = false
        sizes.keys.sort.each_with_index{|size,i|
          next if found
          if list.length < size
            found = true
            sizes[sizes.keys.sort[i-1]] << pos
          end
        }
        sizes[sizes.keys.sort.last] << pos if !found
        
        scores << Score::score(list, ranks.length, 0)[:score]
      else # it has no score
        scores << nil
      end
    }

    info = {}

    # Go through all the size bins, run the permutations and assign the pvalues
    # to all terms in the bin.
    sizes.keys.each{|size|
      next if size == 1
      next if sizes[size].empty?

      # This are the actual scores for the terms in the bin
      sub_list_scores = sizes[size].collect{|pos| scores[pos] || 0}

      # Compute the pvalues for all the terms in the bin. The size of the
      # permutation list is that of the bin
      pvalues = Score::pvalues(sub_list_scores, size, 0, ranks.length)

      # Save the information from the terms, score, hits, and pvalues.
      sizes[size].zip(pvalues).each{|p|
        pos = p[0]
        pvalue = p[1]
        score = scores[pos]
        next if score < 0

        term = best[pos]
        hits = positions[term].nil? ? 0 : positions[term].length 

        info[term] = {:score => score, :hits => hits, :pvalue => pvalue}
      }
    }

    info
  end

  def self.enrichment_hypergeometric(annotations, relevant, options)
    dict_options = {
      :dict_options => {:low => 0, :hi => 0.5, :limit => 100000}
    }.merge(options)[:dict_options]
    positions = {}
    found_datasets = []

    dict = Dictionary::TF_IDF.new
    ranks.each_with_index{|experiment, rank|
      info = annotations[experiment]

      dataset_terms   = info[:dataset]
      signature_terms = info[:signature]
      
      dataset = exp2gds experiment

      terms = signature_terms
      terms += dataset_terms

      term_count = {}
      terms.each{|term|
        term_count[term] ||= 0
        term_count[term] += 1
      }
      dict.add(term_count)
    }
 
    best = dict.best(dict_options).keys

    terms = {}
    found_datasets = []
    annotations.each{|experiment, info|
      dataset_terms   = info[:dataset]
      signature_terms = info[:signature]

      dataset = exp2gds experiment

      signature_terms.each{|term|
        next if ! best.include? term
        terms[term] ||= {:relevant => 0, :total => 0}
        terms[term][:total]    += 1 
        terms[term][:relevant] += 1 if relevant.include? experiment
      }
      
      next if found_datasets.include? dataset
      found_datasets << dataset

      dataset_terms.each{|term|
        next if ! best.include? term
        terms[term] ||= {:relevant => 0, :total => 0}
        terms[term][:total]    += 1 
        terms[term][:relevant] += 1 if relevant.include? experiment
      }
    }
 

    total   = annotations.keys.length
    list    = relevant.length

    terms.each{|term, info|
      info[:pvalue] = Annotations.hypergeometric(total,info[:total],list, info[:relevant])
    }

    terms
  end



  def self.annotations(scores, type, pvalue = 0.05, algorithm = :rank) 
    annot = {}
    relevant = []

    dict_options = {}
    if type == "Words"
      dict_options = {:low => 0, :hi => 0.05, :limit => 100000}
    else
      dict_options = {:low => 0, :hi => 0.5, :limit => 100000}
    end

    case
    when type =~ /^(.*)_direct$/
      side = :direct
      type = $1
    when type =~ /^(.*)_inverse$/
      side = :inverse
      type = $1
    end


    terms_cache = {}
    scores.each{|experiment, info|
      dataset = experiment.match(/^(.*?): /)[1]
      name = $'.strip 
      case
      when side.nil?
        term_file = File.join(MARQ.datadir, MARQ.platform_type(dataset).to_s , 'annotations',type, dataset)
      when side == :direct && info[:score] > 0 || side == :inverse && info[:score] < 0
        term_file = File.join(MARQ.datadir, MARQ.platform_type(dataset).to_s , 'annotations',type + '_up', dataset)
      else
        term_file = File.join(MARQ.datadir, MARQ.platform_type(dataset).to_s , 'annotations',type + '_down', dataset)
      end

      if File.exist? term_file
        terms_cache[term_file] ||= YAML::load(File.open(term_file))
        terms = terms_cache[term_file] 
        annot[experiment] = {:dataset => (terms[:dataset] || []), :signature => (terms[name] || [])}
      else
        annot[experiment] = {:dataset =>  [], :signature => []}
      end

      relevant << experiment if info[:pvalue] <= pvalue
    }

    if algorithm == :rank
      ranks = scores.sort{|a,b| compare(a[1],b[1]) }.collect{|p| p[0]}
      terms = enrichment_rank(annot, ranks, dict_options)
    else
      terms = enrichment_hypergeometric(annot, relevant, dict_options)
    end

    merged_annotations = {}
    annot.each{|key, info|
      merged_annotations[key] = info[:dataset] + info[:signature]
    }
    [merged_annotations, terms]
  end

  module Genes
    module Genecodis
      ORGS = {
        'sgd' => 'Sc' ,
        'rgd' => 'Rn' ,
        'mgi' => 'Mm' ,
        'pombe' => 'Sp' ,
        'cgd' => 'Ca' ,
        'human' => 'Hs' ,
        'tair' => 'At' ,
        'worm' => 'Ce' ,
      }

      FIELDS = %w(Id Items S TS Hyp Hyp_c Genes).collect{|f| f.downcase.to_sym}

      module WS
        def self.driver
          require 'soap/wsdlDriver'
          wsdl_url = File.join('http://genecodis.dacya.ucm.es/static/wsdl/genecodisWS.wsdl')
          driver = SOAP::WSDLDriverFactory.new(wsdl_url).create_rpc_driver
          driver
        end

        def self.analysis(org, list)
          puts "GO for #{ org } #{list.length} genes"

          gc_org = ORGS[org.to_s]
          return [] if gc_org.nil?

          job_id = driver.analyze(gc_org,2,0,-1,3,list,%w(GO_Biological_Process ),[])


          while (stat = driver.status(job_id)) == 1
            sleep 1
          end

          if stat < 0
            return []
          else
            lines =  driver.results(job_id).collect{|l| l.chomp}
            lines.shift
            lines.collect{|l| Hash[*FIELDS.zip(l.chomp.split(/\t/)).flatten]}
          end
        rescue
          puts $!.message
          puts $!.backtrace
        end
      end

      module Local

        def self.analysis(org,list, slim = false)
          require 'genecodis'

          gc_org = ORGS[org.to_s]
          if slim
            groups = ['GOSlim_Process']
          else
            groups = ['GO_Biological_Process']
          end

          job_id = Object::Genecodis.analyze(gc_org,2,0,-1,3,list,groups,nil)
          return [] if job_id.nil?

          while (stat = Object::Genecodis.status(job_id)) == 1
            sleep 0.5
          end

          if stat < 0
            return []
          else
            res =  Object::Genecodis.results(job_id)
            return [] if res.nil?
            lines = res.collect{|l| l.chomp}
            lines.shift
            lines.collect{|l| Hash[*FIELDS.zip(l.chomp.split(/\t/)).flatten]}
          end
        rescue
          puts $!.message
          puts $!.backtrace
        end
      end
    end

    module SENT

      class SENTError < StandardError; end


      WSDL="http://sent.dacya.ucm.es/wsdl/SentWS.wsdl"
      def self.driver
        require 'soap/wsdlDriver'
        driver = SOAP::WSDLDriverFactory.new(WSDL).create_rpc_driver
        driver
      end

      def self.process_results(job)
        result_ids = driver.results(job)

        summary = YAML::load(driver.result(result_ids[0]))
        ccc     = driver.result(result_ids[1]).to_f
        associations = Open.to_hash(StringIO.new(driver.associations(job)), :flatten => true)

        summary.each do |group|
          group[:articles] = group[:genes].inject(0) {|acc, gene| acc += associations[gene].length}
        end

        [summary, ccc]
      end

      def self.analysis(organism, genes, factors)
        hash = Digest::MD5.hexdigest([organism, genes.sort].inspect)

        if @@jobs[hash]
          orig_job = @@jobs[hash]
          job = driver.refactor(orig_job, factors, 'MARQ')
        else
          job = driver.analyze(organism, genes, factors, 'MARQ')
        end

        while ! driver.done(job)
          sleep 5
        end

        raise SENTError "Job failed with error #{driver.messages(job).last}" if driver.error(job)
        @@jobs[hash] = job

        summary, ccc = process_results(job)
      end

      def self.terms(organism, genes, num = 20)
        factor_list = [2,3,4,5,7,10]
        terms = {}
        ccc   = {}
        factor_list.each do |factors|
          summary, ccc = analyze(organism, genes, factors)
          articles = summary.inject(0) {|acc, group| acc += group[:articles] }
          terms_per_article = num.to_f / articles
          summary.each{|group|
            num_terms = terms_per_article * group[:articles]
            terms[factors] ||= []
            terms[factors] += group[:words][num_terms]
          }
          ccc[factors] = ccc
        end
        best_k = ccc.sort_by{|p| p[1]}.first[1]

        terms[k]
      end
    end

    def self.get_genes_nth(dataset, num_genes)
      path = MARQ.dataset_path(dataset)

      experiments = File.open(path + '.experiments').collect{|l| l.chomp.strip}
      genes       = File.open(path + '.codes').collect{|l| l.chomp.strip}
      total_genes = genes.length

      genes_up = {}
      genes_down = {}
      experiments.each{|exp| genes_up[exp] = []; genes_down[exp] = []}

      File.open(path + '.orders').each_with_index{|l, i|
        values = l.chomp.split(/\t/)
        experiments.zip(values).each{|p|
          name = p.first
          value = p.last
          next if p.last == "NA"
          genes_up[name]   <<  genes[i] if value.to_i < num_genes
          genes_down[name] <<  genes[i] if value.to_i > total_genes - num_genes
        }
      }

      {:up => genes_up, :down => genes_down}
    end

    def self.get_genes(dataset, options = {})
      fdr, cut_off, folds, do_folds, nth_genes = {
        :fdr => false,
        :cut_off => 0.05,
        :folds => 2.5,
        :do_folds => true,
        :nth_genes => 0,
      }.merge(options).values_at(:fdr, :cut_off, :folds, :do_folds, :nth_genes)

      if nth_genes > 0
        return get_genes_nth(dataset, nth_genes)
      end
      

      path = MARQ.dataset_path(dataset)

      experiments = File.open(path + '.experiments').collect{|l| l.chomp.strip}
      genes       = File.open(path + '.codes').collect{|l| l.chomp.strip}


      experiments_ts   = experiments.select{|exp| exp !~ /\[ratio\]/}
      experiments_fold = experiments.select{|exp| exp =~ /\[ratio\]/}

      experiments_fold = [] if ! do_folds

      values_up = {}
      values_down = {}
      experiments.each{|exp| values_up[exp] = []; values_down[exp] = []}

      File.open(path + '.pvalues').each_with_index{|l, i|
        values = l.chomp.split(/\t/)
        experiments_ts.zip(values).each{|p|
          name = p.first
          value = p.last == "NA" ? 1.0 : p.last.to_f
          values_up[name]     << (value > 0 ? value : 1.0)
          values_down[name]   << (value < 0 ? - value : 1.0)
        }
      }

      File.open(path + '.logratios').each_with_index{|l, i|
        values = l.chomp.split(/\t/)
        experiments.zip(values).each{|p|
          name = p.first
          next unless experiments_fold.include? name
          value = p.last == "NA" ? 0 : p.last.to_f
          values_up[name]     << (p.last > 0 ? value : 0)
          values_down[name]   << (p.last < 0 ? - value : 0)
        }
      }
     
      genes_up = {}
      genes_down = {}

      threshold = cut_off
      values_up.each{|experiment, values|
        genes_up[experiment] = []
        if experiments_ts.include? experiment
          if fdr
            threshold = FDR.step_up(values.sort, cut_off)
            next if threshold == 0.0
          end
          values.each_with_index{|value, i| genes_up[experiment] << genes[i] if value < threshold}
        elsif experiment_fold.include? experiment
          values.each_with_index{|value, i| genes_up[experiment] << genes[i] if value < folds}
        end
      }
      values_down.each{|experiment, values|
        genes_down[experiment] = []
        if experiments_ts.include? experiment
          if fdr
            threshold = FDR.step_up(values.sort, cut_off)
            next if threshold == 0.0
          end
          values.each_with_index{|value, i| genes_down[experiment] << genes[i] if value < threshold}
        elsif experiment_fold.include? experiment
          values.each_with_index{|value, i| genes_down[experiment] << genes[i] if value < folds}
        end
      }


      {:up => genes_up, :down => genes_down}
    end

    def self.get_genes_old(dataset, cut_off = 0.1, fdr = false)

      path = MARQ.dataset_path(dataset)

      experiments = File.open(path + '.experiments').collect{|l| l.chomp.strip}.select{|name| !name.match(/\[ratio\]/)}
      genes = File.open(path + '.codes').collect{|l| l.chomp.strip}


      values_up = {}
      values_down = {}
      experiments.each{|exp| values_up[exp] = []; values_down[exp] = []}


      if File.exist?(path + '.pvalues')
        File.open(path + '.pvalues').each_with_index{|l, i|
          values = l.chomp.split(/\t/)
          experiments.zip(values).each{|p|
            value = p.last == "NA" ? 1.0 : p.last.to_f
            values_up[p.first]     << (value > 0 ? value : 1.0)
            values_down[p.first]   << (value < 0 ? - value : 1.0)
          }
        }
      end

      genes_up = {}
      genes_down = {}

      threshold = cut_off
      values_up.each{|experiment, values|
        genes_up[experiment] = []
        if fdr
          threshold = FDR.step_up(values.sort, cut_off)
          next if threshold == 0.0
        end
        values.each_with_index{|value, i| genes_up[experiment] << genes[i] if value < threshold}
      }
      values_down.each{|experiment, values|
        genes_down[experiment] = []
        if fdr
          threshold = FDR.step_up(values.sort, cut_off)
          next if threshold == 0.0
        end
        values.each_with_index{|value, i| genes_down[experiment] << genes[i] if value < threshold}
      }


      {:up => genes_up, :down => genes_down}
    end
  end
  module UMLS
    SEMANTIC_TYPES="T020,T100,T116,T123,T023,T118,T043,T049,T103,T200,T060,T047,T203,T126,T050,T131,T125,T129,T037,T197,T191,T114,T110,T167,T024"

    def self.OBA(text)

      res = Net::HTTP.post_form(URI.parse('http://rest.bioontology.org/obs_hibernate/annotator'),
        {
          'longestOnly'=> true, 
          'wholeWordOnly'=> true,
          'withDefaultStopWords' => true,
          'scored' => true,
          'ontologiesToExpand' => 'null', 
          'ontologiesToKeepInResult' => "",
          'levelMax' => 0,
          'levelMin' => 0,
          'textToAnnotate' => text.gsub(/\s/,' '),
          'semanticTypes' => SEMANTIC_TYPES,
          'mappingTypes' => 'null',
          'format' => 'tabDelimited'
      })

      res.body.collect{|l| l.split(/\t/)}.select{|v| v[0].to_i > 0}.collect{|v| v[2].sub(/\W*NOS/,'').downcase}.select{|w| w !~ /^\d+$/}.sort.uniq
    end

  end


end



if __FILE__ == $0
  require 'pp'

  
  
  exit
  #Annotations::GO::Genecodis::Local.init
 
  #genes = Annotations::GO::get_genes('GDS1916')
  #genes[:up].each{|exp, genes|
  #  puts exp
  #  p Annotations::GO::Genecodis::Local.analysis('mgi', genes)
  #}
  #exit
  #res = Annotations::GO::get_genes('GDS948')
  #res[:down].each{|exp, values| puts "#{ exp }\t#{ values.length }"}

  texts = []
  texts << <<-EOT

  Analysis of femurs and tibias of growth hormone (GH) deficient animals at 6
  and 24 hours following treatment with 4 mg/kg body weight GH. Results provide
  insight into the insulin-like growth factor-I dependent and independent
  pathways that mediate the action of GH in bone.
  
  EOT
  texts << <<-EOT
  
  Comparison of total transcription profiles for temperature-sensitive TOR2
  mutant strain SH121 to its isogenic wild type counterpart SH100. Results
  indicate that TOR2 inactivation leads to enhanced transcription of
  Gcn4-controlled target genes.
  

  EOT
  texts << <<-EOT

  Analysis of tissue specimens representing benign nevus, atypical nevus,
  melanoma in situ, vertical growth phase (VGP) melanoma, and metastatic growth
  phase (MGP) melanoma. Results identify expression signatures that distinguish
  benign and atypical nevi and melanomas in situ from VGPs and MGPs.
  
  
  EOT
  texts << <<-EOT

  Analysis of estrogen receptor (ER)-positive MCF7 breast cancer cells up to 48
  hours following treatment with estradiol (E2). ERs facilitate the
  transcriptional effects of hormones. These results, together with ChIP-PET
  results, suggest potential correlations between ER binding and gene
  regulation.

  EOT
  texts << <<-EOT

  
  Analysis of anaerobic chemostat cultures of Saccharomyces cerevisae exposed
  to one of several weak organic acids. Weak organic acids are used as
  preservatives in food and beverages. Yeasts are able to proliferate at the
  maximum legal dosage of such preservatives.

  EOT
  texts << <<-EOT

  Zucker diabetic fatty model of type 2 diabetes: various insulin-sensitive
  tissues Analysis of adipose, skeletal muscle, and liver tissues of Zucker
  diabetic fatty animals at pre-diabetic and diabetic stages. ZDF animals have
  a mutated leptin receptor. Results provide insight into the molecular
  mechanisms responsible for insulin resistance and progression to type 2
  diabetes.

  EOT
  texts << <<-EOT

  Strain differences in copper sulfate Each strain was grown overnight then
  diluted in fresh rich media. After three hours strains were again rediluted
  into either rich media or rich media supplemented with copper sulfate. After
  another three hours cultures were sampled, cells lysed and flash frozen using
  liquid nitrogen. RNA was extrated using hot phenol-chloroform, reverse
  transcripbed using amino-allyle dUTP and labelled with either Cy3 or Cy5
  flourescent dye. Each hybridization is of a single sample compared to a
  reference pool contructed from all the strains with the same treatment.
    Labelled probe were hybridized to DNA microarrays spotted with 6144 70 bp
  oligonucleotides obtained from Qiagen-Operon. After an overnight
  hybridization, microarrays were scanned using a GenePix 4000A scanner and
  spot intensities extracted using GenePix 4.0 software. Bad spots were flagged
  based on the image.


  EOT
  texts << <<-EOT

  MyD88-deficient macrophage response to lipopolysaccharide and E. coli (dye swap)
  Analysis of MyD88 null mutant macrophages treated with LPS or live E. coli.
  MyD88 transduces cell signaling events downstream of Toll-like receptors, a
  key component of host defense. Results suggest most of the host response to
  endotoxin or live bacteria is actually regulated independently of MyD88.

  
  EOT
  texts.reverse.each{|text|
    puts "\n\n--------------\n"
    puts "Text: "
    puts "\n" + text.strip + "\n\n\n"
    puts "Annotations: "
    puts
    puts Annotations::UMLS::OBA(text).join(", ")
  }

#

  #puts  Annotations.hypergeometric(2000,100,100,2)
  #p Annotations::GO::annotate(MARQ.platform_organism('GDS1365'),genes[:up].collect.first.last[1..100])
end

