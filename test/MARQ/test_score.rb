require File.dirname(__FILE__) + '/../helper'
require 'MARQ/score'
require 'MARQ/MADB'

class TestID < Test::Unit::TestCase

  def test_draw
    org       = 'sgd'
    dataset   = 'GDS750'
    signature = 'GDS750: agent: tunicamycin [ratio]'

    dataset  = MARQ::Name.cross_platform dataset
    platform = MARQ::Name.cross_platform MARQ::Dataset.platform dataset

    genes  = load_data(org)[:down]
    native = ID.translate(org, genes)
    translations = Hash[*native.zip(genes).flatten]
    
    positions, matched, entries = MADB.dataset_positions(dataset, native.compact)
    
    Score.draw_hits(positions[signature].compact, entries, '/tmp/hits.jpg')
  end
end
