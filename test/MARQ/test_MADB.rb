require File.dirname(__FILE__) + '/../helper'
require 'lib/MARQ/MADB'

class TestMARQ < Test::Unit::TestCase
  def test_misc
    assert MADB::dataset_positions('GDS113', %w(2778))[0].keys.select{|exp| exp =~ /GDS113/}.any?
    assert MADB::platform_positions('GPL54', %w(2778))[0].keys.select{|exp| exp =~ /GDS113/}.any?

    assert MADB::dataset_positions('GDS113_cross_platform', %w(S000006132))[0].keys.select{|exp| exp =~ /GDS113/}.any?
    assert MADB::platform_positions('GPL54_cross_platform', %w(S000006132))[0].keys.select{|exp| exp =~ /GDS113/}.any?
  end
end

