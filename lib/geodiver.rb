require 'yaml'
require 'fileutils'

require 'geodiver/config'
require 'geodiver/exceptions'
require 'geodiver/logger'
require 'geodiver/routes'
require 'geodiver/server'
require 'geodiver/version'

# GeoDiver NameSpace
module GeoDiver
  class << self
    def environment
      ENV['RACK_ENV']
    end

    def verbose?
      @verbose ||= (environment == 'development')
    end

    def root
      File.dirname(File.dirname(__FILE__))
    end

    def logger
      @logger ||= Logger.new(STDOUT)
    end

    # Setting up the environment before running the app...
    def init(config = {})
      @config = Config.new(config)

      init_dirs
      check_num_threads

      self

      # We don't validate port and host settings. If GeoDiver is run
      # self-hosted, bind will fail on incorrect values. If GeoDiver
      # is run via Apache/Nginx + Passenger, we don't need to worry.
    end

    attr_reader :config, :temp_dir, :public_dir, :users_dir, :db_dir

    # Starting the app manually
    def run
      check_host
      Server.run(self)
    rescue Errno::EADDRINUSE
      puts "** Could not bind to port #{config[:port]}."
      puts "   Is GeoDiver already accessible at #{server_url}?"
      puts '   No? Try running GeoDiver on another port, like so:'
      puts
      puts '       geodiver -p 4570.'
    rescue Errno::EACCES
      puts "** Need root privilege to bind to port #{config[:port]}."
      puts '   It is not advisable to run GeoDiver as root.'
      puts '   Please use Apache/Nginx to bind to a privileged port.'
    end

    def on_start
      puts '** GeoDiver is ready.'
      puts "   Go to #{server_url} in your browser and start analysing GEO datasets!"
      puts '   Press CTRL+C to quit.'
      open_in_browser(server_url)
    end

    def on_stop
      puts
      puts '** Thank you for using GeoDiver :).'
    end

    # Rack-interface.
    #
    # Inject our logger in the env and dispatch request to our controller.
    def call(env)
      env['rack.logger'] = logger
      Routes.call(env)
    end

    private

    # Set up the directory structure in @config[:gd_public_dir]
    def init_dirs
      config[:gd_public_dir] = File.expand_path config[:gd_public_dir]
      logger.debug "GeoDiver Directory: #{config[:gd_public_dir]}"
      @public_dir = create_public_dir
      @users_dir  = File.expand_path('../Users', @public_dir)
      FileUtils.mkdir_p @users_dir unless Dir.exist? @users_dir
      @db_dir     = File.expand_path('../DBs', @public_dir)
      FileUtils.mkdir_p @db_dir unless Dir.exist? @db_dir
    end

    # Create the public dir, if already created and in development environment,
    # remove the existing assets and copy over the new assets
    def create_public_dir
      public_dir = File.join(config[:gd_public_dir], 'public')
      if Dir.exist?(public_dir) && verbose?
        FileUtils.rm_r File.join(public_dir, 'assets')
        FileUtils.cp_r(File.join(GeoDiver.root, 'public/assets'), public_dir)
      else
        FileUtils.mkdir_p public_dir
        FileUtils.cp_r(File.join(GeoDiver.root, 'public/GeoDiver'), public_dir)
        FileUtils.cp_r(File.join(GeoDiver.root, 'public/assets'), public_dir)
      end
      logger.debug "Public Directory: #{public_dir}"
      public_dir
    end

    def check_num_threads
      config[:num_threads] = Integer(config[:num_threads])
      fail NUM_THREADS_INCORRECT unless config[:num_threads] > 0

      logger.debug "Will use #{config[:num_threads]} threads to run GeoDiver."
      if config[:num_threads] > 256
        logger.warn "Number of threads set at #{config[:num_threads]} is" \
                    ' unusually high.'
      end
    rescue
      raise NUM_THREADS_INCORRECT
    end

    # Check and warn user if host is 0.0.0.0 (default).
    def check_host
      return unless config[:host] == '0.0.0.0'
      logger.warn 'Will listen on all interfaces (0.0.0.0).' \
                  ' Consider using 127.0.0.1 (--host option).'
    end

    def server_url
      host = config[:host]
      host = 'localhost' if host == '127.0.0.1' || host == '0.0.0.0'
      "http://#{host}:#{config[:port]}"
    end

    def open_in_browser(server_url)
      return if using_ssh? || verbose?
      if RUBY_PLATFORM =~ /linux/ && xdg?
        system "xdg-open #{server_url}"
      elsif RUBY_PLATFORM =~ /darwin/
        system "open #{server_url}"
      end
    end

    def using_ssh?
      true if ENV['SSH_CLIENT'] || ENV['SSH_TTY'] || ENV['SSH_CONNECTION']
    end

    def xdg?
      true if ENV['DISPLAY'] && system('which xdg-open > /dev/null 2>&1')
    end
  end
end
