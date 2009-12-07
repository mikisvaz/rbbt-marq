require 'test/unit'
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
    GEO.rearange(order,filename)

    text_new = <<-EOF
line0
line4
line1
line3
line2
    EOF

    assert_equal(text_new, Open.read(filename))

  end
end

