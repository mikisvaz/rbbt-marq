require 'rubygems'
require 'test/unit'

$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))

class Test::Unit::TestCase
  DATADIR = File.join(File.dirname(__FILE__), 'data')
  def load_data(org)
    info = YAML::load(File.open(File.join(DATADIR, org + '.yaml')))
  end

  def similar_lists(list1, list2)
    (list1 & list2).length > list1.length * 0.9
  end

end
