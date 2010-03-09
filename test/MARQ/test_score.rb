require File.dirname(__FILE__) + '/../helper'
require 'MARQ/score'
require 'MARQ/MADB'
require 'rand'

class TestID < Test::Unit::TestCase

  def test_fast_score
    org        = 'Hsa'
    dataset    = 'GDS2057'
    experiment = 'agent: RTI-018 <=> control'
    signature  = "#{ dataset }: #{ experiment }"

    genes  = load_data(org)[:down]
    native = ID.translate(org, genes)
    
    scores = Score.scores(MARQ::Name.cross_platform(dataset), native)
  end

  def _test_score
    org        = 'Sce'
    dataset    = 'GDS30'
    experiment = 'time: 90 minute [ratio]'
    signature  = "#{ dataset }: #{ experiment }"

    dataset    = MARQ::Name.cross_platform dataset

    orders = MARQ::Dataset.orders(dataset)[experiment]
    genes  = MARQ::Dataset.codes(dataset)

    sorted_codes = genes.zip(orders).reject{|p| p.last.nil?}.sort_by{|p| p.last}.collect{|p| p.first}
    up   = sorted_codes[0..100]
    down = sorted_codes[(sorted_codes.length - 100)..-1]

    assert Score.scores(dataset, up)[signature][:score] > 0
    assert Score.scores(dataset, down)[signature][:score] < 0

    assert Score.scores_up_down(dataset, up, down)[signature][:score] > 0
    assert Score.scores_up_down(dataset, down, up)[signature][:score] < 0

    null_scores = Score.null_scores(up.length, down.length)

    assert Score.add_pvalues(Score.scores_up_down(dataset, up, down), null_scores)[signature][:pvalue] < 0.05
    assert Score.add_pvalues(
      Score.scores_up_down(dataset, 
                           sorted_codes.shuffle[0..100], 
                           sorted_codes.shuffle[0..100]), 
                           null_scores
    )[signature][:pvalue] > 0.05
  end

  def _test_draw
    org       = 'Sce'
    dataset   = 'GDS750'
    signature = 'GDS750: agent: tunicamycin [ratio]'

    dataset  = MARQ::Name.cross_platform dataset

    genes  = load_data(org)[:down]
    native = ID.translate(org, genes)
    
    positions = MADB.load_positions(dataset, native.compact)
    num_values = MADB.num_values(dataset)
    
    Score.draw_hits(positions[signature].compact, num_values[signature], '/tmp/hits.jpg')
  end
end
