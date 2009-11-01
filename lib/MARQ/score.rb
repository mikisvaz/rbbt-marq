require 'png'
require 'inline'

module Score
  def self.combine(up, down)
    return down if up == 0
    return up if down == 0

    return up - down
    if (up > 0) == (down > 0)
      return 0
    else
      up - down
    end
  end

  def self.average(list)
    clean = list.compact
    clean.inject(0){|acc, e| acc += e}.to_f / clean.length
  end

  def self.score_area(positions, platform_entries, missing = 0)
    return {:score => 0, :top => 0, :bottom => 0} if positions.nil? || positions.empty? || positions.compact.empty?

    clean_positions = positions.compact.sort

    total_tags      = positions.length + missing
    extra           = total_tags - clean_positions.length

    top = 0
    bottom = 0

    clean_positions.each_with_index{|p,i|
      rel_qt = (i + 1).to_f / total_tags
      rel_qb = ( i + extra ).to_f / total_tags
      rel_p  = p.to_f / platform_entries


      top    += rel_qt - rel_p if rel_qt > rel_p
      bottom += rel_p - rel_qb if rel_p > rel_qb
    }


    {
      :top => top,
      :bottom => bottom,
      :score => top > bottom ? top.to_f / total_tags : - bottom.to_f / total_tags, 
    }
  end

  def self.score_max_norm(positions, platform_entries, missing = 0)
    return {:score => 0, :top => 0, :bottom => 0} if positions.nil? || positions.empty? || positions.compact.empty?

    clean_positions = positions.compact.sort

    extra           = missing + (positions.length - clean_positions.length)
    total_tags      = extra + clean_positions.length

    mean = platform_entries / 2

    values_top = [0]
    values_bottom = [0]

    clean_positions.each_with_index{|p,i| 
      rel_qt = (i + 1).to_f / total_tags
      rel_qb = ( i + extra ).to_f / total_tags
      rel_p  = p.to_f / platform_entries


      values_top << (rel_qt - rel_p) * ((p - mean).abs.to_f / mean)**2
      values_bottom << (rel_p - rel_qb) * ((p - mean).abs.to_f / mean)**2
    }

    top    = values_top.max
    bottom = values_bottom.max 


    {
      :score => top > bottom ? top : -bottom, 
    }


  end

  def self.scale_score1(positions, platform_entries)

    mean = platform_entries/2
    max_top = 0
    max_bottom = 0

    top_list = []
    bottom_list = []

    weights = positions.sort.collect{|position| 
      rel_pos = ((position - mean).abs.to_f / mean);
      0.3 * rel_pos + 0.7 * Math::exp(30*rel_pos)/Math::exp(30)
    }
    weights.unshift(0)
    total_weights = weights.inject(0){|v,acc| acc += v}
    weights.collect!{|v| v / total_weights}

    rel_qt = 0
    rel_qb = 0
    positions.sort.each_with_index{|position, idx|

      rel_qt += weights[idx + 1] 
      rel_qb += weights[idx]     
      rel_p   = position.to_f    / platform_entries
 
      top    = (rel_qt - rel_p);
      bottom = (rel_p - rel_qb);

      top_list << top
      bottom_list << bottom

      if (top > max_top) 
        max_top = top;
      end
      if (bottom > max_bottom) 
        max_bottom = bottom;
      end
    }

    p [top_list, bottom_list]
    if (max_top > max_bottom)
      return max_top;
    else
      return -max_bottom;
    end           
  end
  
  class << self
      inline do |builder|
        builder.c_raw <<-'EOC'
    double weight(int position, int mean){
        double rel_pos = (double) abs(position - mean) / mean; 
        double weight =  0.3 *  0.5 * rel_pos +  0.7 * (exp(30*rel_pos)/exp(30));
        return(weight);
    }
        EOC
        builder.c <<-'EOC'
    double fast_score_scale( VALUE positions, int platform_entries, double missing){
      int idx;
    
      int mean = platform_entries / 2;

      VALUE rel_q = rb_ary_new();
      VALUE rel_l = rb_ary_new();
      
      rb_ary_push(rel_q,rb_float_new(0));

      double total_weights = 0;
      for (idx = 0; idx < RARRAY(positions)->len; idx++){
        int position = FIX2INT(rb_ary_entry(positions, idx));

        rb_ary_push(rel_l, rb_float_new((double) position / platform_entries));

        total_weights = total_weights + weight(position, mean);
        rb_ary_push(rel_q,rb_float_new(total_weights));
      }

      // Add penalty for missing genes

      double penalty = missing * weight( mean * 0.8,mean);
      total_weights = total_weights + penalty;
      
      double max_top, max_bottom;
      max_top = max_bottom = 0;
      for (idx = 0; idx < RARRAY(positions)->len; idx++){
        double top    = RFLOAT(rb_ary_entry(rel_q,idx + 1))->value / total_weights -
                        RFLOAT(rb_ary_entry(rel_l,idx))->value;
        double bottom = - (penalty + RFLOAT(rb_ary_entry(rel_q,idx))->value) / total_weights +
                        RFLOAT(rb_ary_entry(rel_l,idx))->value;

        if (top > max_top) max_top = top;
        if (bottom > max_bottom) max_bottom = bottom;
      }
        
     if (max_top > max_bottom) return max_top;
     else                      return -max_bottom;
    }

        EOC



        builder.c <<-'EOC'
    double fast_norm_score( VALUE positions, int total, int extra, int platform_entries){
      int idx;
    
      double mean = (double) platform_entries / 2;
      double max_top, max_bottom;
      max_top = max_bottom = 0;
    
      for (idx = 0; idx < RARRAY(positions)->len; idx++){
        double position  = (double) FIX2INT(rb_ary_entry(positions, (long) idx));


        double rel_qt = (double) (idx + 1)     / total;
        double rel_qb = (double) (idx + extra) / total;
        double rel_p  = position    / platform_entries;
        
        double scale = (abs(position - mean) / mean);
        scale = scale * scale;

        double top    = (rel_qt - rel_p) * scale;
        double bottom = (rel_p - rel_qb) * scale;


        if (top > max_top) max_top = top;
        if (bottom > max_bottom) max_bottom = bottom;
      }

      if (max_top > max_bottom) return max_top;
      else                      return -max_bottom;
    }

        EOC
      end
  end
  def self.score_scale_fast(positions, platform_entries, missing=0)
    return {:score => 0, :top => 0, :bottom => 0} if positions.nil? || positions.empty? || positions.compact.empty?

    clean_positions = positions.compact.sort
    missing = missing + positions.length - clean_positions.length

    {
      :score => fast_score_scale(clean_positions,  platform_entries, missing)
    }
  end


  def self.score_norm_fast(positions, platform_entries, missing = 0)
    return {:score => 0, :top => 0, :bottom => 0} if positions.nil? || positions.empty? || positions.compact.empty?

    clean_positions = positions.compact.sort

    extra           = missing + (positions.length - clean_positions.length)
    total_tags      = extra + clean_positions.length

    {
      :score => fast_norm_score(clean_positions, total_tags, extra, platform_entries)
    }
  end

  def self.score_max(positions, platform_entries, missing = 0)
    return {:score => 0, :top => 0, :bottom => 0} if positions.nil? || positions.empty? || positions.compact.empty?

    clean_positions = positions.compact.sort

    extra           = missing + (positions.length - clean_positions.length)
    total_tags      = extra + clean_positions.length


    values = [0]

    clean_positions.each_with_index{|p,i| 
      
      up   = (i + 1).to_f / total_tags
      down = p.to_f       / platform_entries

      values << up - down
    }

    top    = values.max
    bottom = values.min + ((extra - 1).to_f/ total_tags)


    {
      :score => top.abs > bottom.abs ? top : bottom, 
    }
  end

  class << self
    alias_method :score, :score_scale_fast
  end


  def self.score_up_down(up, down, total, missing_up = 0, missing_down = 0)
    up = score(up, total, missing_up)
    down = score(down, total, missing_down)

    {:up =>  up[:score], :down => down[:score], :score => combine(up[:score], down[:score])}
  end

  def self.permutations(genes, total, times = 10000)
    scores = []
    times.times{
      positions = Array.new(genes){ (rand * total).to_i }
      scores << score(positions, total)[:score]
    }
    scores
  end

  def self.pvalues(scores, up, down, total, options = {})
    times = options[:times]|| 1000

    permutations_up = permutations(up, total, times)
    permutations_down = permutations(down, total, times )
    permutations = permutations_up.zip(permutations_down).collect{|p| combine(*p).abs }.sort

    scores.collect{|score| 
      num = permutations.count_smaller(score.abs)
      (times - num).to_f / times
    }
  end

  COLORS = {
    :red =>   PNG::Color::Red,
    :green => PNG::Color::Green,
    :white => PNG::Color::White,
    :black => PNG::Color::Black,

  }

  def self.draw_hits(hits, total, filename = nil, options = {})

    size = options[:size] || total
    bg_color = options[:bg_color] || :white
    width = options[:width] || 20
    sections = options[:sections] || []

    size = [size, total].min

    hits = hits.collect{|h| h - 1}
    if size < total
      hits = hits.collect{|h| (h.to_f * size / total).to_i}
    end

    canvas = PNG::Canvas.new size, width, COLORS[bg_color]

    sections.each{|color, info|
      start = info[0]
      finish = info[1]
      (start..finish).each{|x|
        (0..width - 1).each{|y|
          canvas[x,y] = COLORS[color]
        }
      }
    }

    hits.each{|hit|
      canvas.line hit, 0, hit , width - 1, PNG::Color::Black
    }

    png = PNG.new canvas

    if filename
      png.save filename
    else
      png.to_blob
    end
  end
end

if __FILE__ == $0
  size = 1000
  positions=%w(10 30 200).collect{|v| v.to_i}
  np = positions.collect{|p| size - p}
  p Score.score(positions, size )
  p Score.score(np, size )


  p Score.scale_score1(positions,  size )
  p Score.scale_score1(np,  size )

  require 'benchmark'


  p = (0..100).collect{ (rand * 1000).to_i}
  puts Benchmark.measure{
    1000.times{|i|
      Score.score_max_norm(p, 1000);
    }
  }
  puts Benchmark.measure{
    1000.times{|i|
      Score.score_max_norm_fast(p, 1000);
    }
  }


  per_list = []
  1000.times{
    per_list << Array.new(200){(rand * 1000).to_i}
  }

  require 'benchmark'
  puts Benchmark.measure{
    per_list.each{|p|
      Score.score_max_norm(p, 1000);
    }
  }

end
