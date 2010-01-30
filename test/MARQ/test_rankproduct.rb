require File.dirname(__FILE__) + '/../helper'
require 'MARQ/main.rb'
require 'MARQ/rankproduct.rb'

class TestRbbtMarq < Test::Unit::TestCase
  def test_rearrage
    dataset = 'GDS113_cross_platform'
    experiment = 'time: 240 minute [ratio]'
    gene = "S000001880"

    assert_equal RankProduct.ranks(dataset, experiment, false)[gene], 
                 RankProduct.ranks(dataset, experiment, true)[gene]
  end
end

