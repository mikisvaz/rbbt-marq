require 'MARQ'

$datadir = MARQ.datadir
$scriptdir = File.join(MARQ.rootdir, '/install_scripts')

task 'rake-includes' do |t|
  FileUtils.mkdir_p($datadir)
  FileUtils.cp File.join($scriptdir, "rake_includes.rb"), $datadir
end

task 'GEO' => 'rake-includes' do |t|
    directory = "#{$datadir}/GEO"
    FileUtils.mkdir_p(directory)
    %w(Rakefile series).each{|f|
      FileUtils.cp_r File.join($scriptdir, "GEO/#{ f }"), directory 
    }
    FileUtils.mkdir(File.join(directory, "platforms")) unless File.exist? File.join(directory, "platforms")
end

task 'CustomDS' => 'rake-includes' do |t|
    directory = "#{$datadir}/CustomDS"
    FileUtils.mkdir_p(directory)
    Dir.glob($scriptdir + '/CustomDS/*').each{|f|
      FileUtils.cp_r  f, directory 
    }
end
