require File.dirname(__FILE__) + '/../helper'
require 'lib/MARQ/MADB'

class TestMARQ < Test::Unit::TestCase
  def test_misc
    assert MADB::dataset_positions('GDS113', %w(2778))[0].keys.select{|exp| exp =~ /GDS113/}.any?
    assert MADB::platform_positions('GPL54', %w(2778))[0].keys.select{|exp| exp =~ /GDS113/}.any?

    assert MADB::dataset_positions('GDS113_cross_platform', %w(S000006132))[0].keys.select{|exp| exp =~ /GDS113/}.any?
    assert MADB::platform_positions('GPL54_cross_platform', %w(S000006132))[0].keys.select{|exp| exp =~ /GDS113/}.any?
  end

  def test_platform_entries
    platform = 'GPL999'

    MARQ::Platform.process(platform) unless MARQ::Platform.exists? platform
    MADB::save_platform_instance(platform)
    MADB::save_platform_instance(MARQ::Name.cross_platform platform)

    assert_equal MARQ::Platform.codes(platform).length, MADB.platform_entries(platform)
    assert_equal MARQ::Platform.cross_platform(platform).length, MADB.platform_entries(MARQ::Name.cross_platform platform)
  end

  def test_positions
    org        = 'sgd'
    dataset    = 'GDS30'
    experiment = 'time: 90 minute [ratio]'
    signature  = "#{ dataset }: #{ experiment }"
    platform   = MARQ::Dataset.platform dataset
    dataset  = MARQ::Name.cross_platform dataset

    GEO::Process.GDS(dataset, platform) unless MARQ::Dataset.exists? dataset

    
    MADB.save_dataset_instance(dataset)

    genes  = load_data(org)[:down]
    native = ID.translate(org, genes)
    translations = Hash[*native.zip(genes).flatten].reject{|k,v| v.nil?}
    native = ["S000001922","S000003506","S000001258","S000001181","S000005815"]
    
    positions, matched, entries = MADB.platform_positions(MARQ::Name.cross_platform(platform), native.compact)
    
    assert matched.any?

    pos       = Hash[*matched.zip(positions[signature]).flatten]
    logratios = MARQ::Dataset.codes_for(dataset, 'logratios', experiment).reject{|k,v| ! native.compact.include? k }
    orders    = MARQ::Dataset.codes_for(dataset, 'orders', experiment).reject{|k,v| ! native.compact.include? k }

    orders_pos = pos.reject{|k,v| v.nil?}.sort_by       {|p| p.last}.collect {|p| p.first}
    orders_log = logratios.reject{|k,v| v.nil?}.sort_by {|p| p.last}.collect {|p| p.first}
    orders_ord = orders.reject{|k,v| v.nil?}.sort_by    {|p| p.last}.collect {|p| p.first}

    assert_equal orders_log, orders_pos.reverse
    assert_equal orders_ord, orders_pos
  end
end

