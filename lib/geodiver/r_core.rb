require 'json'
# GeoDiver NameSpace
module GeoDiver
  # module to run the R core.
  module RCore
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

      def_delegators GeoDiver, :config, :logger, :public_dir, :users_dir,
                     :db_dir

      #
      def init_load_db(params)
        logger.debug('Loading Database')
        assert_load_db_params(params)
      end

      #
      def run_load_db
        copy_exemplar_geo_meta_json
        soft_link_exemplar_geo_to_public_dir
        JSON.parse(IO.read(@public_meta_json))
      end

      ## => Analysis

      #
      def init_analysis(params, user)
        @user = user
        assert_analysis_params(params)
        setup_run_and_public_dir
      end

      #
      def run_analysis
        run_dgea
        soft_link_output_dir_to_public_dir
        generate_relative_results_link
      end

      private

      #
      def assert_load_db_params(params)
        @params = params
        logger.debug("Params: #{@params}")
        assert_geo_db_present
      end

      #
      def assert_geo_db_present
        return unless @params['geo_db'].nil? || @params['geo_db'].empty?
        fail ArgumentError, 'No GEO database provided.'
      end

      # Copy the exemplar meta json file to the Users directory
      def copy_exemplar_geo_meta_json
        exemplar_meta_json_file = File.join(GeoDiver.root, 'exemplar/geo.json')
        @meta_json_file = File.join(db_dir, @params['geo_db'],
                                    "#{@params['geo_db']}.json")
        FileUtils.mkdir(File.join(db_dir, @params['geo_db']))
        FileUtils.cp(exemplar_meta_json_file, @meta_json_file)
      end

      #
      def soft_link_exemplar_geo_to_public_dir
        @public_meta_json = File.join(public_dir, 'GeoDiver/DBs/',
                                      "#{@params['geo_db']}.json")
        FileUtils.ln_s(@meta_json_file, @public_meta_json)
      end

      # => Analysis

      #
      def assert_analysis_params(params)
        @params = params
        assert_geo_db_present
        # {"geo_db"=>"GDS9832", "factor"=>"disease.state",
        # "groupa"=>["Convalescent", "Dengue Fever"],
        # "groupb"=>["Dengue Hemorrhagic Fever", "healthy control"],
        # "dgea"=>"on", "gsea"=>"on"}
      end

      #
      def setup_run_and_public_dir
        @uniq_time = Time.new.strftime('%Y-%m-%d_%H-%M-%S_%L-%N').to_s
        @run_dir = File.join(users_dir, @user.info['email'], @params['geo_db'],
                             @uniq_time)
        FileUtils.mkdir_p(@run_dir)
      end

      #
      def run_dgea
        return unless @params['dgea'] == 'on'
        system(dgea_cmd)
        assert_dgea_output
      end

      #
      def dgea_cmd
        "Rscript #{File.join(GeoDiver.root, 'RCore/DGEA.R')}" \
        " --accession GDS5093 --factor '#{@params['factor']}'" \
        " --popA '#{@params['groupa'].join(',')}'" \
        " --popB '#{@params['groupb'].join(',')}'" \
        " --popname1 'Dengue' --popname2 'Normal'" \
        ' --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005' \
        " --working_dir '#{@run_dir}/'"
      end

      #
      def assert_dgea_output
        true
      end

      def soft_link_output_dir_to_public_dir
        public_user_dir = File.join(public_dir, 'GeoDiver/Users',
                                    @user.info['email'], @params['geo_db'])
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
