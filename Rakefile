require 'bundler/gem_tasks'
require 'rspec/core'
require 'rspec/core/rake_task'

task default: [:build]

desc 'Builds and installs'
task install: [:installRdependencies, :build] do
  require_relative 'lib/geodiver/version'
  sh "gem install ./geodiver-#{GeoDiver::VERSION}.gem"
end

desc 'Runs tests and builds gem (default)'
task build: [:test] do
  sh 'gem build geodiver.gemspec'
end

desc 'Install R dependencies'
task :installRdependencies do
  sh 'Rscript RCore/Installations.R'
end

task test: :spec
RSpec::Core::RakeTask.new(:spec) do |spec|
  spec.pattern = FileList['spec/**/*_spec.rb']
end
