require File.dirname(__FILE__) + '/../helper'
require 'MARQ'
require 'MARQ/main'
require 'MARQ/ID'

class TestMARQ < Test::Unit::TestCase
  def check_results(scores, positive, negative)
    negative.each do |signature|
      assert scores[signature][:pvalue] < 0.01
      assert scores[signature][:score]  < 0
    end
    positive.each do |signature|
      assert scores[signature][:pvalue] < 0.01
      assert scores[signature][:score]  > 0
    end

  end

  def check_organism(organism)
    info = load_data(organism)
    up = ID.translate(organism, info[:up]).compact
    down = ID.translate(organism, info[:down]).compact

    require 'progress-monitor'

    Progress.monitor("Querying #{ organism }", :stack_depth => 1, :skip => 1, :announcement => Proc.new{|p| "Platform #{ p }"} )
    scores = MARQ::RankQuery.organism_scores(organism, up, down)

    check_results(scores, info[:positive], info[:negative])
  end

  def test_organisms
    %w(mgi sgd tair rgd human).each do |organism|
      puts "Testing #{ organism }"
      check_organism(organism)
    end
  end
end

