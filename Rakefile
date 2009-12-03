require 'rubygems'
require 'rake'

begin
  require 'jeweler'
  Jeweler::Tasks.new do |gem|
    gem.name = "rbbt-marq"
    gem.summary = %Q{MicroArray Rank Query}
    gem.description = %Q{Find microarray experiments with similar or opposite signature to a given query. A SOAP interface and a merb portal can be found in rbbt-marq-www.}
    gem.email = "miguel.vazquez@fdi.ucm.es"
    gem.homepage = "http://github.com/mikisvaz/rbbt-marq"
    gem.authors = ["Miguel Vazquez"]

      

    gem.files = Dir['lib/**/*.rb','bin/marq_config','tasks/install.rake', 'install_scripts/**/*','R/*']
    gem.test_files = Dir['lib/**/test_*.rb']

    gem.add_dependency('rbbt')
    gem.add_dependency('rake')
    gem.add_dependency('simpleconsole')
    gem.add_dependency('DBcache')
    gem.add_dependency('DRbServe')
    gem.add_dependency('png')


    # FIXME: This wont work with ruby 1.9
    #
    # gem.add_dependency('RubyInline')
     

    # FIXME: Rubygems install does not work right away
    #
    # gem.add_dependency('rsruby')


    # gem is a Gem::Specification... see http://www.rubygems.org/read/chapter/20 for additional settings
  end
rescue LoadError
  puts "Jeweler (or a dependency) not available. Install it with: sudo gem install jeweler"
end

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  test.libs << 'lib' << 'test'
  test.pattern = 'test/**/test_*.rb'
  test.verbose = true
end

begin
  require 'rcov/rcovtask'
  Rcov::RcovTask.new do |test|
    test.libs << 'test'
    test.pattern = 'test/**/test_*.rb'
    test.verbose = true
  end
rescue LoadError
  task :rcov do
    abort "RCov is not available. In order to run rcov, you must: sudo gem install spicycode-rcov"
  end
end

task :test => :check_dependencies

task :default => :test

require 'rake/rdoctask'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "rbbt-marq #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
