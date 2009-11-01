require 'fileutils'
require 'yaml'
require 'digest/md5'

module MARQ

  class MARQ::NoConfigError < StandardError; end

  def self.rootdir
    File.dirname(File.dirname(__FILE__))
  end

  config_file = [
    File.join(ENV['HOME'], ".#{ self }"),
    File.join(rootdir, "#{self}.config")
  ].select{|f| File.exist? f}.first

  raise MARQ::NoConfig if config_file.nil?

  @@config = YAML::load(File.open(config_file))

  class << self
    %w(datadir workdir cachedir dbhost dbname dbuser dbpass).each{|dir|
      define_method(dir, proc{@@config[dir]})
    }
    %w(datadir workdir cachedir).each{|dir|
      FileUtils.mkdir_p(@@config[dir]) unless (@@config[dir].nil? || File.exists?(@@config[dir].to_s))
    }

  end

end

class String
  def hash_code
    Digest::MD5.hexdigest(self)
  end
end

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


require 'DBcache'

DBcache.config({
  :dbuser => MARQ.dbuser,
  :dbpass => MARQ.dbpass,
  :dbhost => MARQ.dbhost,
  :dbname => MARQ.dbname,
})


