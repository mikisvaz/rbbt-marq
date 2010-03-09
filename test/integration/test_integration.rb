require File.dirname(__FILE__) + '/../helper'
require 'MARQ'
require 'MARQ/main'
require 'MARQ/ID'
require 'progress-monitor'

class TestMARQ < Test::Unit::TestCase
  def check_results(scores, positive, negative)
    negative.each do |signature|
      puts "--- (-)"
      p signature
      p scores[signature]
      puts
      puts
      assert scores[signature][:pvalue] < 0.05
      assert scores[signature][:score]  < 0
    end
    positive.each do |signature|
      puts "--- (+)"
      p signature
      p scores[signature]
      puts
      puts
      assert scores[signature][:pvalue] < 0.05
      assert scores[signature][:score]  > 0
    end

  end

  def check_organism(organism)
    info = load_data(organism)
    up   = ID.translate(organism, info[:up]).compact
    down = ID.translate(organism, info[:down]).compact


    Progress.monitor("Querying #{ organism }", :stack_depth => 1, :announcement => Proc.new{|p| "Platform #{ p }"} )
    scores = MARQ::RankQuery.organism_scores(organism, up, down)
    scores = MARQ::RankQuery.add_pvalues(scores, up.length, down.length)

    check_results(scores, info[:positive], info[:negative])
  end

  def test_organisms
    %w(Sce Ath Hsa Mmu Rno).each do |organism|
      puts "Testing #{ organism }"
      check_organism(organism)
    end
  end
end

