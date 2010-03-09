require 'png'
require 'inline'

module Score

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
    double fast_score_scale(VALUE positions, int total, int missing){
      int idx;
    
      int mean = total / 2;

      VALUE rel_q = rb_ary_new();
      VALUE rel_l = rb_ary_new();
      
      rb_ary_push(rel_q,rb_float_new(0));

      // Rescale positions and accumulate weights
      double total_weights = 0;
      for (idx = 0; idx < RARRAY(positions)->len; idx++){
        int position = FIX2INT(rb_ary_entry(positions, idx));

        rb_ary_push(rel_l, rb_float_new((double) position / total));

        total_weights += weight(position, mean);
        rb_ary_push(rel_q, rb_float_new(total_weights));
      }

      // Add penalty for missing genes
      double penalty = missing * weight(mean * 0.8, mean);
      total_weights  = total_weights + penalty;
      
      // Traverse list and get extreme values
      double max_top, max_bottom;
      max_top = max_bottom = 0;
      for (idx = 0; idx < RARRAY(positions)->len; idx++){
        double top    = RFLOAT(rb_ary_entry(rel_q, idx + 1))->value / total_weights -
                        RFLOAT(rb_ary_entry(rel_l, idx))->value;
        double bottom = - (penalty + RFLOAT(rb_ary_entry(rel_q, idx))->value) / total_weights +
                        RFLOAT(rb_ary_entry(rel_l, idx))->value;

        if (top > max_top)       max_top    = top;
        if (bottom > max_bottom) max_bottom = bottom;
      }
        
     if (max_top > max_bottom) return max_top;
     else                      return -max_bottom;
    }
        EOC

      end

  end


  def self.score(*args)
    self.fast_score_scale(*args)
  end

  def self.scores(dataset, genes)
    positions   = MADB.load_positions(dataset, genes)
    values      = MADB.num_values(dataset)

    experiments = positions.keys

    scores      = {}
    experiments.each do |experiment|
      hits  = positions[experiment].compact
      total = values[experiment]
      if hits.nil? || hits.empty?
        score = 0
      else
        missing = genes.length - hits.length
        score = self.fast_score_scale(hits.sort, total, missing) 
      end
      scores[experiment] = {
        :positions => positions[experiment],
        :score     => score,
        :total     => total,
      }
    end

    scores
  end

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

  def self.scores_up_down(dataset, up, down)
    scores_up   = scores(dataset, up)
    scores_down = scores(dataset, down)

    scores      = {}
    scores_up.keys.each do |experiment|
      scores[experiment] = {}
      scores[experiment][:up]    = scores_up[experiment]
      scores[experiment][:down]  = scores_down[experiment]
      scores[experiment][:score] = combine(scores_up[experiment][:score], scores_down[experiment][:score])
    end

    scores
  end

  def self.permutations(size, times)
    total = 10000
    if size == 0
      [0] * times
    else
      (1..times).collect do
        fast_score_scale(Array.new(size){ (rand * total).to_i }.sort, total, 0)
      end
    end
  end

  def self.null_scores(up_size, down_size, times = 10000)
    up_perm   = permutations(up_size, times)
    down_perm = permutations(down_size, times)

    up_perm.zip(down_perm).collect{|p| up, down = p; combine(up, down).abs}
  end

  def self.add_pvalues(scores, null_scores)
    null_scores = null_scores.sort
    times       = null_scores.length

    scores.each do |experiment, info|
      info[:pvalue] = (times - null_scores.count_smaller(info[:score].abs)).to_f / times
    end

    scores
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
