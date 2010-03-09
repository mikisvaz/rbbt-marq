require File.dirname(__FILE__) + '/../helper'
require 'lib/MARQ/MADB'

class TestMARQ < Test::Unit::TestCase
  def test_positions
    MADB::load_positions('GDS113', %w(2778)).each do |exp, list|
      assert list.any?
    end
  end

  def test_num_codes
    platform = 'GPL999'

    MARQ::Platform.process(platform) unless MARQ::Platform.exists? platform
    MADB::save_platform_instance(platform)
    MADB::save_platform_instance(MARQ::Name.cross_platform platform)

    assert_equal MARQ::Platform.codes(platform).length, MADB.num_codes(platform)
    assert_equal MARQ::Platform.cross_platform(platform).length, MADB.num_codes(MARQ::Name.cross_platform platform)
  end

  def test_positions
    org        = 'Sce'
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

    positions = MADB.load_positions(MARQ::Name.cross_platform(dataset), native.compact)

    assert positions[signature].compact.any?

    pos       = Hash[*native.compact.zip(positions[signature]).flatten]
    logratios = MARQ::Dataset.codes_for(dataset, 'logratios', experiment).reject{|k,v| ! native.compact.include? k }
    orders    = MARQ::Dataset.codes_for(dataset, 'orders', experiment).reject{|k,v| ! native.compact.include? k }

    orders_pos = pos.reject{|k,v| v.nil?}.sort_by       {|p| p.last}.collect {|p| p.first}
    orders_log = logratios.reject{|k,v| v.nil?}.sort_by {|p| p.last}.collect {|p| p.first}
    orders_ord = orders.reject{|k,v| v.nil?}.sort_by    {|p| p.last}.collect {|p| p.first}

    assert_equal orders_ord, orders_pos
    assert_equal orders_log, orders_pos.reverse
  end

  def test_num_values
    dataset = 'GDS113'

    assert_equal MADB::num_codes(dataset), MADB::num_codes(MARQ::Dataset.platform dataset)

    assert_equal MADB::num_values(dataset).keys.sort, 
                 MARQ::Dataset.experiments(dataset).collect{|exp| "#{ dataset }: #{ exp }"}.sort

    assert MADB::num_values(dataset).values.first <= MADB::num_codes(dataset)
    
  end
end

