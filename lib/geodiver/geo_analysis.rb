require 'forwardable'
require 'json'

require 'geodiver/pool'

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

      def_delegators GeoDiver, :config, :logger, :public_dir, :users_dir,
                     :db_dir

      def_delegators GeoDiver::LoadGeoData, :load_geo_db
      #
      def init(params, user)
        @user = user
        @params = params
        assert_params
        setup_run_and_public_dir
        save_params
        @uniq_time
      end

      #
      def run
        # wait until geo db has been loaded in background thread
        load_geo_db.join unless load_geo_db.nil?
        run_analysis
        soft_link_output_dir_to_public_dir
        generate_relative_results_link(@uniq_time)
      end

      def get_expression_json(params)
        @params = params
        assert_gene_id_present
        run_dir = run_expression_analysis
        File.join(run_dir, "dgea_#{@params[:gene_id]}.json")
      end

      def create_interactions(params)
        @params = params 
        run_interaction_analysis 
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

      def assert_gene_id_present
        logger.debug('Asserting Gene ID is present.')
        return unless @params['gene_id'].nil? || @params['gene_id'].empty?
        fail ArgumentError, 'No Gene Id provided.'
      end

      #
      def setup_run_and_public_dir
        @uniq_time = Time.new.strftime('%Y-%m-%d_%H-%M-%S_%L-%N').to_s
        @run_dir = File.join(users_dir, @user.info['email'], @params['geo_db'],
                             @uniq_time)
        logger.debug("Creating Run Directory: #{@run_dir}")
        FileUtils.mkdir_p(@run_dir)
      end

      def save_params
        File.open(File.join(@run_dir, 'params.json'), 'w') do |f|
          f.puts @params.to_json
        end
      end

      def run_analysis
        p = Pool.new(config[:num_threads]) if config[:num_threads] > 1
        if config[:num_threads] > 1
          p.schedule { run_overview }
          p.schedule { run_dgea } if @params['dgea'] == 'on'
          p.schedule { run_gage } if @params['gsea'] == 'on'
        else
          run_overview
          run_dgea if @params['dgea'] == 'on'
          run_gage if @params['gsea'] == 'on'
        end
        assert_overview_output
        assert_dgea_output if @params['dgea'] == 'on'
        assert_gsea_output if @params['gsea'] == 'on'
      ensure
        p.shutdown if config[:num_threads] > 1
      end

      def run_overview
        logger.debug("Running CMD: #{overview_cmd}")
        system(overview_cmd)
      end

      def overview_cmd
        "Rscript #{File.join(GeoDiver.root, 'RCore/overview.R')}" \
        " --dbrdata #{dbrdata} --rundir '#{@run_dir}/'" \
        " --analyse 'Boxplot,PCA'" \
        " --accession #{@params['geo_db']} --factor '#{@params['factor']}'" \
        " --popA '#{@params['groupa'].join(',')}'" \
        " --popB '#{@params['groupb'].join(',')}'" \
        " --popname1 'Dengue' --popname2 'Normal'" \
        " --distance '#{@params['heatmap_distance_method']}'" \
        " --clustering '#{@params['heatmap_clustering_method']}'" \
        ' --dev TRUE'
      end

      #
      def assert_overview_output
        true
      end

      #
      def run_dgea
        logger.debug("Running CMD: #{dgea_cmd}")
        system(dgea_cmd)
      end

      #
      def dgea_cmd
        "Rscript #{File.join(GeoDiver.root, 'RCore/dgea.R')}" \
        " --dbrdata #{dbrdata} --rundir '#{@run_dir}/'" \
        " --analyse 'Toptable,Heatmap,Volcano'" \
        " --accession #{@params['geo_db']} --factor '#{@params['factor']}'" \
        " --popA '#{@params['groupa'].join(',')}'" \
        " --popB '#{@params['groupb'].join(',')}'" \
        " --popname1 'Dengue' --popname2 'Normal'" \
        " --topgenecount #{@params['number_top_genes']} " \
        ' --foldchange 0.3 --thresholdvalue 0.005' \
        " --distance '#{@params['heatmap_distance_method']}'" \
        " --clustering '#{@params['heatmap_clustering_method']}'" \
        " --heatmaprows #{@params['heatmap_rows']} " \
        " --adjmethod '#{@params['volcano_pValue_cutoff']}'" \
        " --dendrow #{(@params['cluster_by_genes'] == 'on')} "\
        " --dendcol #{(@params['cluster_by_samples'] == 'on')} "\
        ' --dev TRUE'
      end

      #
      def assert_dgea_output
        true
      end

      def run_gage
        logger.debug("Running CMD: #{gsea_cmd}")
        system(gsea_cmd)
      end

      def gsea_cmd
        "Rscript #{File.join(GeoDiver.root, 'RCore/gage.R')}" \
        " --dbrdata #{dbrdata} --rundir '#{@run_dir}/'" \
        " --accession #{@params['geo_db']} --factor '#{@params['factor']}'" \
        " --popA '#{@params['groupa'].join(',')}'" \
        " --popB '#{@params['groupb'].join(',')}'" \
        " --comparisontype 'ExpVsCtrl' --genesettype 'KEGG' --geotype 'BP'" \
        ' --dev TRUE'
      end

      def assert_gsea_output
        true
      end

      def dbrdata
        File.join(db_dir, @params['geo_db'], "#{@params['geo_db']}.RData")
      end

      def soft_link_output_dir_to_public_dir
        public_user_dir = File.join(public_dir, 'GeoDiver/Users',
                                    @user.info['email'], @params['geo_db'])
        logger.debug("Creating a Soft Link: #{@run_dir} ==> #{public_user_dir}")
        FileUtils.mkdir_p(public_user_dir) unless Dir.exist? public_user_dir
        FileUtils.ln_s(@run_dir, public_user_dir)
      end

      def generate_relative_results_link(uniq_time)
        File.join('GeoDiver/Users/', @user.info['email'], @params['geo_db'],
                  uniq_time)
      end

      def run_expression_analysis
        run_dir = File.join(users_dir, @user.info['email'], @params['geo_db'],
                            @params['result_id'])
        cmd = expression_cmd(run_dir)
        logger.debug("Running CMD: #{cmd}")
        system(cmd)
        assert_expression_output
        run_dir
      end

      def expression_cmd(run_dir)
        "Rscript #{File.join(GeoDiver.root, 'RCore/dgea_expression.R')}" \
        " --rundir '#{run_dir}/' --geneid '#{@params[:gene_id]}'"
      end

      def assert_expression_output
        true
      end

      def run_interaction_analysis
        run_dir = File.join(users_dir, @user.info['email'], @params['geo_db'],
                            @params['result_id'])
        cmd = interaction(run_dir)
        logger.debug("Running CMD: #{cmd}")
        system(cmd)
        assert_expression_output
        generate_relative_results_link(@params['result_id'])
      end

      def interaction(run_dir)
        "Rscript #{File.join(GeoDiver.root, 'RCore/gage_interaction_network.R')}" \
        " --rundir '#{run_dir}/' --pathid '#{@params[:path_id]}'"
      end
    end
  end
end
