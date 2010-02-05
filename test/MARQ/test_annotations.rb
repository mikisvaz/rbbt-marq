require File.dirname(__FILE__) + '/../helper'
require 'lib/MARQ/annotations'

class TestMARQ < Test::Unit::TestCase
  def test_dataset_annotations
    dataset = 'GDS113_cross_platform'
    experiment = "time: 15 minute [ratio]"

    assert Annotations.dataset_annotations(dataset, 'Words', experiment)[:dataset].any?
    assert Annotations.dataset_annotations(dataset, 'GO_up', experiment)[:signature].any?
    assert Annotations.dataset_annotations(dataset, 'GO_up', experiment)[:dataset].empty?
  end


  def test_annotations
  end

  def test_GO_local
    assert Annotations::Genes::Genecodis::Local.analysis('sgd', %w(YHL004W YBR251W YGL123W YNL178W YLR441C YML063W YHR203C)).any?
  end

  def _test_SENT
    #summary, ccc = Annotations::Genes::SENT::process_results('Example_Human')
    #assert ccc > 0.5
    #assert summary.length >= 2


    dataset = 'GDS113_cross_platform'
    codes  = MARQ::Dataset.codes(dataset)
    orders = MARQ::Dataset.orders(dataset).values.first
    genes  = codes.zip(orders).
      reject  {|p| p[1].nil? }.
      sort_by {|p| p[1]}[0..99].
      collect  {|p| p[0]}

    terms  = Annotations::Genes::SENT.terms(MARQ::Dataset.organism(dataset), genes)
  end
end

