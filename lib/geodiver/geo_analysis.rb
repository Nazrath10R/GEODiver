require 'json'
# GeoDiver NameSpace
module GeoDiver
  # Module to run the GEO analysis
  module GeoAnalysis
    # To signal error in query sequence or options.
    #
    # ArgumentError is raised when ... exit status is 1; see [1].
    class ArgumentError < ArgumentError
    end

    # To signal internal errors.
    #
    # RuntimeError is raised when there is a problem in writing the input file,
    # running R Script, writing the output etc. These are rare, infrastructure
    # errors, used internally, and of concern only to the admins/developers.
    # One example of a RuntimeError would be R libraries not installed.
    class RuntimeError < RuntimeError
    end

    class << self
      extend Forwardable

      def_delegators GeoDiver, :logger, :public_dir, :users_dir, :db_dir

      def_delegators GeoDiver::LoadGeoData, :load_geo_db
      #
      def init(params, user)
        @user = user
        @params = params
        assert_params
        setup_run_and_public_dir
      end

      #
      def run
        run_dgea
        soft_link_output_dir_to_public_dir
        generate_relative_results_link
      end

      private

      #
      def assert_params
        assert_geo_db_present
      end

      #
      def assert_geo_db_present
        logger.debug('Asserting GEO db is present.')
        return unless @params['geo_db'].nil? || @params['geo_db'].empty?
        fail ArgumentError, 'No GEO database provided.'
      end

      #
      def setup_run_and_public_dir
        @uniq_time = Time.new.strftime('%Y-%m-%d_%H-%M-%S_%L-%N').to_s
        @run_dir = File.join(users_dir, @user.info['email'], @params['geo_db'],
                             @uniq_time)
        logger.debug("Creating Run Directory: #{@run_dir}")
        FileUtils.mkdir_p(@run_dir)
      end

      #
      def run_dgea
        return unless @params['dgea'] == 'on'
        # wait until geo db has been loaded in background thread
        load_geo_db.join unless load_geo_db.nil?
        logger.debug("Running CMD: #{dgea_cmd}")
        system(dgea_cmd)
        assert_dgea_output
      end

      #
      def dgea_cmd
        dbrdata = File.join(db_dir, @params['geo_db'],
                            "#{@params['geo_db']}.Rdata")
        "Rscript #{File.join(GeoDiver.root, 'RCore/DGEA.R')}" \
        " --dbrdata #{dbrdata} --accession #{@params['geo_db']} " \
        " --factor '#{@params['factor']}'" \
        " --popA '#{@params['groupa'].join(',')}'" \
        " --popB '#{@params['groupb'].join(',')}'" \
        " --popname1 'Dengue' --popname2 'Normal'" \
        ' --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005' \
        " --outputdir '#{@run_dir}/'"
      end

      #
      def assert_dgea_output
        true
      end

      def soft_link_output_dir_to_public_dir
        public_user_dir = File.join(public_dir, 'GeoDiver/Users',
                                    @user.info['email'], @params['geo_db'])
        logger.debug("Creating a Soft Link: #{@run_dir} ==> #{public_user_dir}")
        FileUtils.mkdir_p(public_user_dir) unless Dir.exist? public_user_dir
        FileUtils.ln_s(@run_dir, public_user_dir)
      end

      def generate_relative_results_link
        File.join('GeoDiver/Users/', @user.info['email'], @params['geo_db'],
                  @uniq_time)
      end
    end
  end
end
