require File.dirname(__FILE__) + '/../helper'
require 'lib/MARQ/annotations'

class TestMARQ < Test::Unit::TestCase
  def test_misc
    assert ! Annotations.dataset_annotations('GSE966', 'Words', nil).empty?
  end

  def _test_GO_local
    assert Annotations::Genes::Genecodis::Local.analysis('sgd', %w(YHL004W YBR251W YGL123W YNL178W YLR441C YML063W YHR203C)).any?
  end

  def test_SENT
    #summary, ccc = Annotations::Genes::SENT::process_results('Example_Human')
    #assert ccc > 0.5
    #assert summary.length >= 2


    dataset = 'GDS113_cross_platform'
    codes  = MARQ::Dataset.codes(dataset)
    orders = MARQ::Dataset.orders(dataset)
    genes  = codes.zip(orders).sort_by{|p| p[1]}[1..100].select {|p| p[0]}

    terms  = Annotations::Genes::SENT.terms(MARQ::Dataset.organism(dataset), genes)
    p terms
  end
end

