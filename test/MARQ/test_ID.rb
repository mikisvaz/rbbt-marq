require File.dirname(__FILE__) + '/../helper'
require 'MARQ/ID'
require 'test/unit'

class TestID < Test::Unit::TestCase

  def test_index
    assert_equal(%w(1020), ID.translate('human', %w(CDK5)))
  end

end
