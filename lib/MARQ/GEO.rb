require 'MARQ'
require 'MARQ/main'
require 'rbbt/sources/organism'

# Work with GEO datasets
module GEO

  CACHE_DIR = File.join(MARQ.cachedir,'GEO')
  FileUtils.mkdir_p CACHE_DIR unless File.exists? CACHE_DIR

  DATA_DIR = File.join(MARQ.datadir, 'GEO')

  PLATFORM_BLACKLIST = %w(GPL4065)

  # Get information from Entrez
  module Remote

    @@nice = 1
    def self.organism_platforms(org)
      name = Organism.name(org)
      Open.read("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=\"#{name}\"[Organism:exp]+AND+%22gpl%22[Filter]&retmax=10000", :nice => @@nice).
        scan(/<Id>(\d+?)<\/Id>/).collect{|id| id.first}.collect{|id| "GPL#{id.sub(/^100*/,'')}"}
    end

    def self.platform_datasets(platform)
      Open.read("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=#{platform}[Accession]&retmax=2000", :nice => @@nice).
      scan(/<Id>(\d+?)<\/Id>/).collect{|id| id.first}.select{|id| !id.match(/^(1|2)000/) }.collect{|id| "GDS#{id}"}
    end

    def self.dataset_platform(dataset)
      if dataset =~ /GSE/
        Open.read("http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=#{dataset}", :nice => @@nice).scan(/GPL\d+/).uniq.sort.join("_")
      else
        Open.read("http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=#{dataset}", :nice => @@nice).scan(/GPL\d+/).uniq.sort.join("_")
      end
    end

    def self.series_dataset?(gse)
      Open.read("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=geo&term=#{gse}[Accession]&retmax=2000", :nice => @@nice).
      match(/<Id>(\d+?)<\/Id>/) != nil
    end

    def self.platform_organism(platform)
      Open.read("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=#{platform}", :nice => @@nice).
        scan(%r#<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi\?mode=Info&amp;id=\d+" onmouseout="onLinkOut\('HelpMessage' , geo_empty_help\)" onmouseover="onLinkOver\('HelpMessage' , geoaxema_organismus\)">(.*?)</a>#).collect{|p| p.first}.join(', ')
    end

  end

  # Parse information in .soft files
  module SOFT

    @@nice = 1
    GEO_SOFT="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&view=full&form=text&acc="

    # Download a soft file. Uses cache
    def self.get_soft(item)
      item = item.strip
      cache_file = File.join(CACHE_DIR, item + '.soft')
      if File.exist?( cache_file )
        File.open(cache_file).read
      else
        content = Open.read(GEO_SOFT + item, :nocache => true, :nice => @@nice)
        raise "SOFT file error" if content !~ /!/
        fout = File.open(cache_file,'w')
        fout.write content
        fout.close
        content
      end
    end

    #{{{ Guess the format of the IDS

    @@formats = {}

    ID_FIX = {
      :mgi_unigene => proc{|gene| if gene then gene.match(/^Mm./) ? gene : "Mm." + gene end},
      :human_unigene => proc{|gene| if gene then gene.match(/^Hs./) ? gene : "Hs." + gene end},
    }

    # Id list is in sequence
    def self.consecutive?(ids)
      ids.collect{|id| id.to_i}.sort[0..19] == (1..20).to_a
    end

    # Id list is numerical
    def self.numerical?(ids)
      ids.compact.select{|id| ! id.match(/^\d+$/)}.uniq.length < ids.length.to_f / 10
    end

    # ID are DNA bases
    def self.dna_sequence?(ids)
      ids.compact.select{|id| ! id.strip.match(/^[ATCG]+$/i)}.empty?
    end

    # Guess the format of the id in the list. The name parameter can be used to
    # identify some exceptions
    def self.guessIds(genes,org, name = nil)
      @@formats[org] ||= Organism.id_formats(org)
      if consecutive?(genes) || dna_sequence?(genes) || (numerical?(genes) && (name.nil? || !name.match(/entrez/i)))
        id = nil
      else
        fix = ID_FIX[(org + "_" + name.downcase).to_sym] if name
        if fix
          genes = genes.collect{|gene| fix.call(gene)}
        end
        id = Organism.guessIdFormat(@@formats[org], genes)
      end

      id 
    end


    #{{{ Parse soft files for several GEO entities

    def self.GSE(series)
      soft = get_soft(series)

      # Find platform
      if match = soft.scan(/!Series_platform_id\s*=?\s*(.*)/)
        platform = match.flatten.collect{|p| p.strip}.join("_")
      else
        raise "No Platform information" 
      end

      # Find title
      if soft.match(/!Series_title \s*=?\s*(.*)/)
        title = $1
      else
        raise "No Title information" 
      end

      # Find summary
      if soft.match(/!Series_summary \s*=?\s*(.*)/)
        matches = soft.scan(/!Series_summary \s*=?\s*(.*)/).to_a
        description = matches.collect{|m| m.to_s.strip.sub(/!Series_summary \s*=?\s*/,'')}.join("\n")
      else
        raise "No Summary information" 
      end

      # Find samples
      if soft.match(/!Series_sample_id \s*=?\s*(.*)/)
        matches = soft.scan(/!Series_sample_id \s*=?\s*(.*)/).to_a
        samples = matches.collect{|m| m.to_s.strip.sub(/!Series_sample_id \s*=?\s*/,'')}
      else
        raise "No Summary information" 
      end

      {
        :platform => platform,
        :description =>description.strip,
        :title => title.strip,
        :samples => samples,
      }
    end

    def self.GSM(array)
      soft = get_soft(array)

      # Find title
      if soft.match(/!Sample_title\s*=?\s*(.*)/)
        title = $1
      else
        raise "No Title information" 
      end


      # Find description
      if soft.match(/!Sample_description \s*=?\s*(.*)/)
        description = $1
      else
        raise "No Description information" 
      end

      {

        :description =>description.strip,
        :title => title.strip,
      }
    end

    def self.GPL(platform)
      if !File.exist?(File.join(DATA_DIR, 'platforms',"#{platform}.yaml"))  &&
         !File.exist?(File.join(DATA_DIR, 'platforms',"#{platform}.skip")) 
        begin
          raise "Platform Blacklisted" if PLATFORM_BLACKLIST.include? platform

          if platform =~ /_/
            organism =  GPL(platform.match(/(.*?)_/)[1])[:organism]

            info =  {
              :organism => organism,
              :title => "Merged platforms #{ platform }",
            }
            return info
          end
          soft = get_soft(platform)


          raise "SOFT file error" if soft !~ /!/

            organisms = soft.scan(/!Platform_organism\s*=\s*(.*)/).collect{|v| v.first.strip}

          if organisms.empty?
            raise "No Organism information" 
          else
            # This might happen actually GPL2529
            organisms.delete('Schizosaccharomyces pombe') if  organisms.include?('Saccharomyces cerevisiae')
            org_name = organisms.first
          end


          title = ""
          if soft.match(/!Platform_title\s*=\s*(.*)/)
            title = $1
          end

          org = Organism.name2org(org_name)
          raise "Organism not identified: #{org_name}" if org.nil?

          if soft.match(/!platform_table_begin/)
            data = soft.split(/!platform_table_begin/s)[1].collect{|l| l.chomp.split(/\t/)}
            data.shift
            names = data.shift
            total = data.first.length
            genes = data.sort_by{ rand }[1..1000].collect{|v| v.first}

            id = guessIds(genes,org, names.first)
            other = nil
            other_pos = 0
            other_count = 0
            other_name = 0
            if id.nil? 
              (1..total - 1).to_a.each{|num|
                genes = data.collect{|v| v[num]}
                other = guessIds(genes,org, name = names[num])

                if other && other[1] > other_count
                  other_pos = num
                  other_count = other[1]
                  other_name = names[num]
                end
              }
            end
          else
            raise "Soft file incomplete"
          end

          info = {:organism => org, :BioMart_ID => id ? id.first : nil, :title => title }
          info[:other_ID_field] = [other_pos + 1, other_name] if other_pos > 0


          Open.write(File.join(DATA_DIR, 'platforms',"#{platform}.yaml"), info.to_yaml)
        rescue Exception
          puts $!.message
          puts $!.backtrace
          Open.write(File.join(DATA_DIR, 'platforms',"#{platform}.skip"), $!.message)
        end
      end

      raise "Platform info for #{ platform } is not available and could not be automatically produced." if File.exist?(File.join(MARQ.datadir, 'GEO', 'platforms',"#{platform}.skip"))

      YAML::load(File.open(File.join(DATA_DIR, 'platforms',"#{platform}.yaml")))
    end


  end


  #{{{ Process

  # Use R to load and process the datasets 
  module Process

    class PlatformNotProcessedError < StandardError; end
    class AdhocPlatformCollisionError < StandardError; end

    # R library wrapper
    module R
      @@r = nil

      # Get the R instance
      def self.r
        if @@r.nil?

          # FIXME: RSruby does not install very well, this require id hidden here.
          require 'rsruby'

          RSRuby.instance.source(MARQ.rootdir + '/R/MA.R')
          RSRuby.instance.source(MARQ.rootdir + '/R/GEO.R')
          RSRuby.instance.source(MARQ.rootdir + '/R/GEOquery_patch.R')
          @@r = RSRuby.instance
        end
        @@r
      end

      # Use R to load GPL info
      def self.GPL(name, prefix, id_field = nil)
        r.GEO_GPL_process(name, prefix, id_field, CACHE_DIR)
      end

      # Use R to load process the dataset
      def self.GDS(name, prefix, id_field = nil, id_file = nil)
        r.GEO_GDS_process(name, prefix, id_field, id_file, CACHE_DIR)
      end

      # Use R to load process the series
      def self.GSE(gsms, conditions, do_log, prefix, id_file = nil, fields= nil, title = nil, description = nil)
        r.GEO_GSE_process(gsms, conditions, prefix, do_log, id_file, fields, title, description, CACHE_DIR)
      end
    end

    def self.translate(org, list)
      begin
        ID.translate_DB(org, list)
      rescue
        puts "DB translation failed, resorting to index"
        ID.translate_index(org, list)
      end
    end

    # Rearange the lines of a file with the given order. The order specifies, for
    # each position in the original file, where it should en in the final file
    def self.rearange(order, file, missing = "NA")
      orig_lines = []
      File.open(file).each_line{|l| orig_lines << l}

      return if orig_lines.empty?
      columns = orig_lines.first.split(/\t/).length

      lines = Array.new(order.length)

      orig_lines.each_with_index{|l,i|
        next if order[i].nil?
        lines[order[i]] = l.chomp
      }

      lines = lines.collect{|l| l || [missing]*columns*"\t"}

      fout = File.open(file, 'w')
      fout.puts(lines.join("\n"))
      fout.close
    end

    # Fix possible discrepancies in ids between series and platforms
    def self.fix_GSE_ids(platform_codes_file, prefix)
      platform_codes = File.open(platform_codes_file).collect{|l| l.chomp}
      platform_order = {}

      platform_codes.each_with_index{|code, i|
        platform_order[code] = i
      }

      series_codes = File.open(prefix + '.codes').collect{|l| l.chomp}

      platform_positions = platform_order.values_at(*series_codes)

      # Fill with nil for missing positions
      platform_positions[platform_codes.length - 1] ||= nil

      %w(t logratios orders pvalues).each{|ext|
        rearange(platform_positions, prefix + '.' + ext)
      }

      FileUtils.cp(platform_codes_file, prefix + '.codes')
      Open.write(prefix + '.swap', platform_positions.join("\n"))
    end

    def self.GSE(series, info)
      platform = info[:platform]
      do_log   = info[:log2].nil? ? nil : !info[:log2]
      fields   = info[:fields]
 
      # Determine samples and sample conditions
      gsms = []
      conditions = {}
      info[:arrays].each{|gsm, cond|
        gsms << gsm
        cond.each{|type, value|
          conditions[type] ||= []
          conditions[type] << value
        }
      }
     
      # Adhoc platforms are for series with samples from different platforms.
      # They are created when the series is processed
      adhoc_platform = platform.match(/_/) != nil

      raise PlatformNotProcessedError if ! adhoc_platform && ! MARQ::Platform.exists?(platform)
      raise AdhocPlatformCollisionError if adhoc_platform && MARQ::Platform.exists?(platform) && ! MARQ::Name.is_cross_platform?(series)

      platform_path = GEO.platform_path(platform) || File.join(MARQ.datadir, 'GEO', platform)

      prefix = File.join(platform_path, 'GSE', series)
      
      FileUtils.mkdir_p File.join(platform_path, 'GSE') unless File.exists? File.join(platform_path, 'GSE')

      FileUtils.rm(prefix + '.skip') if File.exist?(prefix + '.skip')

      if ! MARQ::Name.is_cross_platform?(series)
        R.GSE(gsms, conditions, do_log, prefix, nil, fields, info[:title], info[:description])

        # Set up codes and cross_platform for adhoc platforms
        if adhoc_platform 
          codes        = Open.read(prefix + '.codes').collect{|l| l.chomp}
          organism     = GEO.platform_organism(platform.split(/_/)[0])
          translations = translate(organism, codes) 
          FileUtils.cp(prefix + '.codes', File.join(platform_path,'codes'))
          Open.write(File.join(platform_path, 'translations'), translations.collect{|v| v || "NO MATCH"}.join("\n"))
          Open.write(File.join(platform_path, 'cross_platform'), translations.compact.sort.uniq.join("\n"))
        else
          fix_GSE_ids(File.join(platform_path, 'codes'), prefix);
        end

      else
        FileUtils.cp(File.join(platform,'translations'), prefix + '.translations')

        swap_file = File.join(MARQ::Dataset.path(MARQ::Name.clean series) + '.swap')
        if File.exist?(swap_file)
          orders = Open.read(swap_file).collect{|l| l.chomp}
          inverse_orders = Array.new(orders.length)
          orders.each_with_index{|pos,i|
            next if pos !~ /\d/
              inverse_orders[pos.to_i] = i
          }
          rearange(inverse_orders, prefix + '.translations', "NO MATCH")
        end

        R.GSE(gsms, conditions, do_log, prefix, prefix +  '.translations', fields, info[:title], info[:description])

        fix_GSE_ids(File.join(platform_path, 'cross_platform'),prefix);

        FileUtils.rm(prefix + '.translations') if File.exist?(prefix + '.translations')
      end

    end

    # Process a dataset. Need to specify the platform. The field parameter can
    # be used to use a different column for the field. 
    #
    # Deprecated in favor of using the original firt column and using a
    # different one only for translation
    def self.GDS(dataset, platform, field = nil)
      raise PlatformNotProcessedError if ! MARQ::Platform.exists? platform

      cross_platform = MARQ::Name.is_cross_platform? dataset

      platform_path = GEO.platform_path(platform)
      prefix = File.join(platform_path, 'GDS', dataset)

      FileUtils.rm(prefix + '.skip') if File.exist?(prefix + '.skip')

      if cross_platform
        R.GDS(MARQ::Name.clean(dataset), prefix, field, File.join(platform_path, 'translations'))
      else
        R.GDS(dataset, prefix, field, nil) 
      end
    end


    # Load GPL data. Translates IDS of the platform probes using AILUN and our
    # system (called biomart for clarity)
    def self.GPL(platform)
      path = File.join(DATA_DIR, platform)
      return if File.exist?(path)

      if platform =~ /_/
        FileUtils.mkdir(path)
        FileUtils.mkdir(path + '/GSE')
        FileUtils.mkdir(path + '/GDS')
        return
      end

      info = SOFT.GPL(platform)
      organism = info[:organism]

      field = info[:other_ID_field]
      id = info[:BioMart_ID]
      org = info[:organism]
      field = nil if field == ""
      id = nil if id == ""


      puts "Processing Platform #{ platform }"
      [platform,
        File.join(path, 'GDS'),
        File.join(path, 'GSE'),
      ].each{|d|
        FileUtils.mkdir d unless File.exist? d
      }

      R.GPL(platform, path, nil)
      FileUtils.mv path + '.codes', File.join(path, 'codes')


      # AILUN translations
      codes  = Open.read(File.join(path, 'codes')).collect{|l| l.chomp}
      ailun = ID.AILUN_translate(platform, codes)
      Open.write(File.join(path, 'ailun'), ailun.collect{|v| v || "NO MATCH"}.join("\n")) if ailun.compact.length > codes.length.to_f / 10

      # BioMart translations
      biomart = []
      if id || field
        if id
          codes  = Open.read(File.join(path, 'codes')).collect{|l| l.chomp}
        else
          if field
            R.GPL(platform, path, field[0])
            FileUtils.mv path + '.codes', File.join(path, 'other')
          end

          fix = GEO::SOFT::ID_FIX[(organism + "_" + field[1].downcase).to_sym]
          codes  = Open.read(File.join(path, 'other')).collect{|l| 
            code = l.chomp
            code = fix.call(code) if fix
            code
          }
        end

        biomart = translate(organism, codes) 
        Open.write(File.join(path, 'biomart'), biomart.collect{|v| v || "NO MATCH"}.join("\n")) if biomart.compact.length > codes.length.to_f / 10
      end

      # Select Best and save
      translations = []
      if ailun.compact.uniq.length > biomart.compact.uniq.length
        id_type = ID::DEFAULT_FORMATS[organism] || ID::DEFAULT_FORMAT_ALL || id || field || "Entrez Gene Id"
        if id_type.to_s !~ /Entrez/i
          translations = translate(org,ailun.collect{|gene| gene || "NO MATCH"}) 
        else
          translations = ailun 
        end
      else
        translations = biomart
      end

      if translations.compact.length > codes.length.to_f / 10        
        Open.write(File.join(path, 'translations'), translations.collect{|v| v || "NO MATCH"}.join("\n"))
        Open.write(File.join(path, 'cross_platform'), translations.compact.sort.uniq.join("\n"))
      end

    end
  end

  def self.platforms
    Dir.glob(File.join(DATA_DIR, "GPL*")).collect {|path| File.basename(path) }
  end


  def self.dataset_type(dataset)
    case
    when dataset =~ /^GDS/
      :GDS
    when dataset =~ /^GSE/
      :GSE
    end
  end

  def self.platform_path(platform)
    path = File.join(DATA_DIR, platform)
    path = nil unless File.exists? File.join(path, 'codes')
    path
  end

  def self.dataset_path(dataset, platform = nil)
    if platform
      platforms = [platform]
    else
      platforms = self.platforms
    end

    platforms.each do |platform|
      platform_path = platform_path(platform)

      next if platform_path.nil?

      prefix = File.join(platform_path, dataset_type(dataset).to_s, dataset)
      
      if File.exists?(prefix + '.orders') || File.exists?(prefix + '.skip')
        return File.join(platform_path, dataset_type(dataset).to_s, dataset)
      end

    end

    return nil
  end

  def self.platform_datasets(platform)
    cross_platform = MARQ::Platform.is_cross_platform? platform
    
    path = platform_path(MARQ::Name.clean(platform))
    return [] if path.nil?

    datasets = Dir.glob(File.join(path, '*', '*.orders')).
      collect {|path| File.basename(path).sub(/\.orders$/,'')}

    if cross_platform
      datasets.select  {|dataset| MARQ::Dataset.is_cross_platform? dataset }.
               collect {|dataset| MARQ::Name.clean(dataset) }
    else
      datasets.select  {|dataset| ! MARQ::Dataset.is_cross_platform? dataset }
    end
  end

  def self.dataset_platform(dataset)
    path = dataset_path(dataset)
    return nil if path.nil?
    path.match(/(GPL\d+)/)
    return $1
  end

  def self.platform_organism(platform)
    GEO::SOFT.GPL(platform)[:organism]
  end

  def self.dataset_organism(dataset)
    platform_organism(dataset_platform(dataset))
  end

  def self.process_platform(platform)
    GEO::Process.GPL(platform) unless platform =~ /_/
  end

  def self.process_dataset(dataset, platform)
    case dataset_type(dataset)
    when :GDS
      GEO::Process.GDS(dataset, platform)
    when :GSE
      info = YAML::load(File.open("series/#{ MARQ::Name.clean dataset }.yaml"))
      FileUtils.rm("platforms/#{ info[:platform] }.skip") if File.exist?  "platforms/#{ info[:platform] }.skip"
      GEO::Process.GSE(dataset, info)
    end
  end

end


if __FILE__ == $0
  p GEO.dataset_path 'GDS2791_cross_platform', 'GPL96'

end

