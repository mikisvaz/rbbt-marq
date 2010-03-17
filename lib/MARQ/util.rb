class Array
  def count_smaller(value)
    size = self.length

    low = 0
    last = size - 1
    hi = size - 1
    pos = size / 2
    while true
      case
      when value == self[pos]
        return pos 
      when pos == last && self[pos] < value
        return size
      when pos == 0 && value <= self[pos] 
        return 0 
      when  self[pos] < value && value < self[pos + 1]
        return pos + 1 
      when  value < self[pos]
        hi = pos 
        pos = (pos - low) / 2 + low
      when self[pos] < value
        low = pos
        pos = (hi - pos) / 2 + pos + 1
      end
    end
  end
end


