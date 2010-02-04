require File.dirname(__FILE__) + '/../helper'
require 'MARQ'
require 'MARQ/main'
require 'MARQ/ID'

class TestMARQ < Test::Unit::TestCase
  DATADIR = File.join(File.dirname(__FILE__), 'data')
  def load_data(org)
    info = YAML::load(File.open(File.join(DATADIR, org + '.yaml')))
  end

  def check_results(scores, positive, negative)
    positive.each do |signature|
      p signature
      assert scores[signature][:pvalue] < 0.01
      assert scores[signature][:score]  > 0
    end
    negative.each do |signature|
      p signature
      assert scores[signature][:pvalue] < 0.01
      assert scores[signature][:score]  < 0
    end
  end

  def check_organism(organism)
    info = load_data(organism)
    up = ID.translate(organism, info[:up]).compact
    down = ID.translate(organism, info[:down]).compact

    require 'progress-monitor'

    Progress.monitor("Querying #{ organism }", :stack_depth => 1, :announcement => Proc.new{|p| "Platform #{ p }"} )
    scores = MARQ::RankQuery.organism_scores(organism, up, down)

    check_results(scores, info[:positive], info[:negative])
  end

  def test_organisms
    %w(tair sgd  rgd human mgi).each do |organism|
      puts "Testing #{ organism }"
      check_organism(organism)
    end
  end
end

