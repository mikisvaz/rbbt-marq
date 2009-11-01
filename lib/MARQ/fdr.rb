require 'inline'

module FDR

  # values should be sorted
  def self.step_up_native(values, rate)
    total = values.length

    last = 0
    values.each_with_index  do |value, i|
      if value > rate * (i + 1).to_f / total
        return last
      end
      last = value
    end
    return last
  end

  # values should be sorted
  def self.adjust_native(values)
    total = values.length.to_f

    adjusted = []
    last = 1
    values.reverse.each_with_index do |value, i|
      adj = [last, value * total / (total - i )].min
      last = adj
      adjusted << adj
    end

    adjusted.reverse
  end

  class << self
    inline do |builder|
      builder.c <<-EOC
        double step_up_fast(VALUE ps, double rate){
           long idx;
           int total = RARRAY(ps)->len;
  
           double last_value = 0;
           for (idx = 0; idx < total; idx++){
             double p  = (double) RFLOAT(rb_ary_entry(ps, idx))->value;
             
             if (p > rate * (double) (idx + 1) / (double) total){
                return last_value;
             }
             last_value = p;
           }

           return last_value;
        }

      EOC


      builder.c <<-EOC
          
         VALUE adjust_fast_self(VALUE ps){
           long idx;
      
           int total = RARRAY(ps)->len;

           VALUE new = rb_ary_new();

           double last = 1;
           for (idx = total - 1; idx >= 0 ; idx-- ){
             double p  = (double) RFLOAT(rb_ary_entry(ps, idx))->value;

      
             p = p * (double) total / (double) (idx + 1);
             if (p > last) p = last;
             last = p;

             RFLOAT(rb_ary_entry(ps, idx))->value = p;
           }

          return ps;
         }
      EOC
         
      builder.c <<-EOC
         VALUE adjust_fast(VALUE ps){
           long idx;
      
           int total = RARRAY(ps)->len;

           VALUE new = rb_ary_new();

           double last = 1;
           for (idx = total - 1; idx >= 0 ; idx-- ){
             double p  = (double) RFLOAT(rb_ary_entry(ps, idx))->value;

      
             p = p * (double) total / (double) (idx + 1);
             if (p > last) p = last;
             last = p;

             rb_ary_unshift(new,rb_float_new(p));
           }

           return new;
         }
      EOC
    end
  end

  class << self
    alias_method :adjust, :adjust_fast
    alias_method :adjust!, :adjust_fast_self
    alias_method :step_up, :step_up_fast
  end

  # This will change the values of the floats in situ
  def self.adjust_hash!(data, field = nil)
    keys = []
    values = []

    data.collect{|key, value| [key, field.nil? ? value : value[field]] }.sort{|a,b| 
      a[1] <=> b[1]
    }.each{|p|
      keys << p[0]
      values << p[1]
    }

    FDR.adjust!(values)

    data
  end

end

if __FILE__ == $0

  require 'test/unit'
  require 'rsruby'

  class FDRTest < Test::Unit::TestCase
    def clean(values)
      if Array === values
        values.collect{|v| (v * 10000).to_i.to_f / 10000}
      else
        (values * 10000).to_i.to_f / 10000
      end
    end

    def copy(values)
      values.collect{|v| v + 0.0}
    end

    def setup
      @r = RSRuby.instance
      @values = [0.001, 0.002, 0.003, 0.003, 0.003, 0.004, 0.006, 0.07, 0.09]
      @threshold = 0.01
      @r_adj = @r.p_adjust(@values,'BH')

    end

    def test_step_up
      assert_equal(0.006, clean(FDR.step_up(@values, @threshold)))
      assert_equal(clean(FDR.step_up_native(@values, @threshold)), clean(FDR.step_up_fast(@values,@threshold)))
      assert_equal(@r_adj.select{|v| v <= @threshold}.length, @values.select{|v| v <= FDR.step_up(@values, @threshold)}.length)
    end

    def test_adjust
      assert_equal(clean(@r_adj), clean(FDR.adjust_native(@values)))
      assert_equal(clean(FDR.adjust_fast(@values)), clean(FDR.adjust_native(@values)))

      assert_equal(clean(@r_adj), clean(FDR.adjust_fast_self(copy(@values))))
    end

  end

  require 'test/unit/ui/console/testrunner'
  Test::Unit::UI::Console::TestRunner.run(FDRTest)

end
