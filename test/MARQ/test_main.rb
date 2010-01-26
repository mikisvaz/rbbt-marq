require File.dirname(__FILE__) + '/../helper'
require 'lib/MARQ/main'

class TestMARQ < Test::Unit::TestCase
  def test_misc
    assert_equal('sgd', MARQ::Dataset.organism('GDS113'))
    assert_equal('GPL54', MARQ::Dataset.platform('GDS113'))
  end

  def test_score
    assert MARQ::RankQuery.dataset_scores('GDS113_cross_platform',%w(),%w()).empty?
    assert MARQ::RankQuery.platform_scores('GPL54',%w(),%w()).empty?
  end
end

