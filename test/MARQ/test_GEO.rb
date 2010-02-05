require File.dirname(__FILE__) + '/../helper'
require 'rbbt/util/tmpfile'
require 'rbbt/util/open'
require 'lib/MARQ/GEO.rb'

class TestRbbtMarq < Test::Unit::TestCase
  def test_rearrage
    text_orig = <<-EOF
line0
line1
line2
line3
line4
    EOF
    filename = TmpFile.tmp_file
    Open.write(filename, text_orig)

    order = [0,2,4,3,1]
    GEO::Process.rearange(order,filename)

    text_new = <<-EOF
line0
line4
line1
line3
line2
    EOF

    assert_equal(text_new, Open.read(filename))

  end

  def test_platform_organism
    assert_equal "Homo sapiens", GEO::Remote.platform_organism('GPL570')
    assert_equal "Schizosaccharomyces pombe, Saccharomyces cerevisiae", GEO::Remote.platform_organism('GPL2529')
  end

  def test_misc
    assert GEO::platform_datasets('GPL54').include?('GDS113')
    assert GEO::platform_datasets('GPL54_cross_platform').include?('GDS113')
  end

  def valid_orders(dataset)
    logratios = MARQ::Dataset.logratios(dataset)
    t         = MARQ::Dataset.t(dataset)
    orders    = MARQ::Dataset.orders(dataset)
    
    MARQ::Dataset.experiments(dataset).each do |experiment|
      codes     = MARQ::Dataset.codes(dataset)
      values    = MARQ::Name.is_ratio?(experiment) ? logratios[experiment] : t[experiment]
      ord       = orders[experiment]

      codes_val = codes.zip(values).reject{|p| p.last.nil? }.sort_by{|p| p.last }.collect {|p| p.first }
      codes_ord = codes.zip(ord).reject{|p| p.last.nil? }.sort_by{|p| p.last }.collect {|p| p.first }
     
      assert similar_lists(codes_val.reverse[1..100], codes_ord[1..100])
      assert similar_lists(codes_val[1..100], codes_ord.reverse[1..100])
      assert_equal ord.compact.length, ord.compact.max
    end

  end

  def test_process_gds
    dataset    = 'GDS962'
    platform   = GEO::Remote.dataset_platform dataset

    GEO::Process.GDS(dataset, platform)
    valid_orders(dataset)
    
    if MARQ::Platform.has_cross_platform? platform
      GEO::Process.GDS(MARQ::Name.cross_platform(dataset), platform) 
      valid_orders(MARQ::Name.cross_platform(dataset))
    end
  end

  def test_process_gse
    dataset    = 'GSE966'
    info       = YAML::load <<-EOF
--- 
:arrays: 
  GSM15048: 
    condition: wt2hZYM
  GSM15049: 
    condition: wt2hZYM
  GSM15050: 
    condition: wt2hZYM
  GSM15051: 
    condition: wt2hZYM
:description: "S. cerevisiae was grown on YEPD. For Zymolyase experiments, yeast cells were grown overnight at 24 \xC2\xB0C to an optical density 0.8 - 1 (A600). The culture was refreshed to 0.2 O.D and grown at 24 \xC2\xB0C for 2h 30min. Next, culture was divided into two parts. One continues growing under same conditions (non-treated culture) while the other was supplemented with 5units/ml Zymolyase 100T . Cells were collected at 2 hours of growth, frozen at -80 \xC2\xB0C and processed for RNA extraction.\n\
  Keywords: repeat sample"
:title: wild type two hours exposure to zymolyase
:platform: GPL764
    EOF

    GEO::Process.GSE(dataset, info)
    valid_orders(dataset)
    
    if MARQ::Platform.has_cross_platform? info[:platform]
      GEO::Process.GSE(dataset, info) 
      valid_orders(MARQ::Name.cross_platform dataset)
    end
  end

end

