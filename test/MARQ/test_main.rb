require File.dirname(__FILE__) + '/../helper'
require 'lib/MARQ/main'

class TestMARQ < Test::Unit::TestCase
  
  def test_basic
    assert ! MARQ::Dataset.exists?('FAKE')
    assert ! MARQ::Dataset.broken?('FAKE')
    
    assert MARQ::Dataset.exists?('GDS113')
    assert ! MARQ::Dataset.broken?('GDS113')

    assert MARQ::Dataset.exists?('GDS113_cross_platform')
    assert ! MARQ::Dataset.broken?('GDS113_cross_platform')
    
    assert MARQ::Platform.has_cross_platform?(MARQ::Dataset.platform('GDS113_cross_platform'))

    assert MARQ::Name.is_ratio? 'GDS750: genotype/variation: ADH1proHAC1 [ratio]'
  end
  
  def test_misc
    assert_equal('sgd', MARQ::Dataset.organism('GDS113'))
    assert_equal('GPL54', MARQ::Dataset.platform('GDS113'))
  end

  def test_score
    assert MARQ::RankQuery.dataset_scores('GDS113_cross_platform',%w(),%w()).empty?
    assert MARQ::RankQuery.platform_scores('GPL54',%w(),%w()).empty?
  end

  def test_codes_for
    dataset = 'GDS750_cross_platform'
    experiment = 'agent: tunicamycin [ratio]'
    MARQ::Dataser.process(dataset) unless MARQ::Dataset.exists? dataset
    assert_equal 64.165, MARQ::Dataset.codes_for(dataset, 'logratios', experiment)['S000001922']
  end
end

