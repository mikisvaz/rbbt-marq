require File.dirname(__FILE__) + '/../helper'
require 'rbbt/util/tmpfile'
require 'rbbt/util/open'
require 'lib/MARQ/GEO.rb'

class TestRbbtMarq < Test::Unit::TestCase
  def test_rearrage
    text_orig = <<-EOF
line0
line1
line2
line3
line4
    EOF
    filename = TmpFile.tmp_file
    Open.write(filename, text_orig)

    order = [0,2,4,3,1]
    GEO::Process.rearange(order,filename)

    text_new = <<-EOF
line0
line4
line1
line3
line2
    EOF

    assert_equal(text_new, Open.read(filename))

  end

  def test_platform_organism
    assert_equal "Homo sapiens", GEO::Remote.platform_organism('GPL570')
    assert_equal "Schizosaccharomyces pombe, Saccharomyces cerevisiae", GEO::Remote.platform_organism('GPL2529')
  end

  def test_misc
    assert GEO::platform_datasets('GPL54').include?('GDS113')
    assert GEO::platform_datasets('GPL54_cross_platform').include?('GDS113')
  end

end

