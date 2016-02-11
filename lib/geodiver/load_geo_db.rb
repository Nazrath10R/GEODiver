require 'json'
# GeoDiver NameSpace
module GeoDiver
  # module to run the Load the GEO dataset.
  module LoadGeoData
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

      def_delegators GeoDiver, :logger, :public_dir, :db_dir

      attr_accessor :load_geo_db

      #
      def init(params)
        logger.debug('Loading Database')
        assert_params(params)
      end

      #
      def run
        @meta_json = File.join(db_dir, @params['geo_db'],
                               "#{@params['geo_db']}.json")
        meta_data = (File.exist?(@meta_json)) ? parse_meta : download_meta_data
        soft_link_meta_json_to_public_dir
        logger.debug("Meta Data: #{meta_data}")
        convert_geo_db_into_r_objects
        meta_data
      end

      private

      #
      def assert_params(params)
        @params = params
        assert_geo_db_present
      end

      #
      def assert_geo_db_present
        logger.debug('Checking if the GEO DB parameter is present.')
        return unless @params['geo_db'].nil? || @params['geo_db'].empty?
        fail ArgumentError, 'No GEO database provided.'
      end

      def parse_meta
        logger.debug("Parse the Meta JSON file at: #{@meta_json}")
        JSON.parse(IO.read(@meta_json))
      end
      #
      def download_meta_data
        file = download_geo_file
        data = read_geo_file(file)
        data = parse_geo_db(data)
        write_to_json(data, @meta_json)
        data
      end

      #
      def download_geo_file
        remote_dir = generate_remote_url
        output_dir = File.join(db_dir, @params['geo_db'])
        FileUtils.mkdir(output_dir) unless Dir.exist? output_dir
        compressed = File.join(output_dir, "#{@params['geo_db']}.soft.gz")
        logger.debug("Downloading from: #{remote_dir} ==> #{compressed}")
        system "wget #{remote_dir} --output-document #{compressed} " \
               '>/dev/null 2>&1'
        logger.debug("Uncompressing file: #{compressed.gsub('.gz', '')}")
        system "gunzip --force -c #{compressed} > #{compressed.gsub('.gz', '')}"
        compressed.gsub('.gz', '')
      end

      #
      def generate_remote_url
        if @params['geo_db'].length == 6
          remote_dir = 'ftp://ftp.ncbi.nlm.nih.gov//geo/datasets/GDSnnn/' \
                       "#{@params['geo_db']}/soft/#{@params['geo_db']}.soft.gz"
        else
          dir_number = @params['geo_db'].match(/GDS(\d)\d+/)[1]
          remote_dir = 'ftp://ftp.ncbi.nlm.nih.gov//geo/datasets/' \
                       "GDS#{dir_number}nnn/#{@params['geo_db']}/soft/" \
                       "#{@params['geo_db']}.soft.gz"
        end
        remote_dir
      end

      #
      def read_geo_file(file)
        data = []
        IO.foreach(file) do |line|
          break if line =~ /^#ID_REF/
          data << line
        end
        data.join
      end

      #
      def parse_geo_db(data)
        {
          'Accession' => data.match(/\^DATASET = (.*)/)[1],
          'Title' => data.match(/!dataset_title = (.*)/)[1],
          'Description' => data.match(/!dataset_description = (.*)/)[1],
          'Sample_Organism' => data.match(/!dataset_platform_organism = (.*)/)[1],
          'Factors' => parse_factors(data),
          'Reference' => data.match(/!Database_ref = (.*)/)[1],
          'Update_Date' => data.match(/!dataset_update_date = (.*)/)[1]
        }
      end

      #
      def parse_factors(data)
        subsets = data.gsub(/\^DATA.*\n/, '').gsub(/\![dD]ata.*\n/, '')
        results = {}
        subsets.lines.each_slice(5) do |subset|
          desc = subset[2].match(/\!subset_description = (.*)/)[1]
          type = subset[4].match(/\!subset_type = (.*)/)[1].gsub(' ', '.')
          # samples = subset[3].match(/\!subset_sample_id = (.*)/)[1]
          results[type] ||= []
          results[type] << desc
        end
        results
      end

      #
      def write_to_json(hash, output_json)
        logger.debug("Writing meta data to file: #{output_json}")
        File.open(output_json, 'w') { |f| f.puts hash.to_json }
      end

      #
      def soft_link_meta_json_to_public_dir
        public_meta_json = File.join(public_dir, 'GeoDiver/DBs/',
                                     "#{@params['geo_db']}.json")
        logger.debug("Creating a Soft Link from: #{@meta_json} ==>" \
                     " #{public_meta_json}")
        return if File.exist? public_meta_json
        FileUtils.ln_s(@meta_json, public_meta_json)
      end

      def convert_geo_db_into_r_objects
        return if File.exist?(File.join(db_dir, @params['geo_db'],
                                        "#{@params['geo_db']}.Rdata"))
        logger.debug("Running: #{load_geo_db_cmd}")
        @load_geo_db = Thread.new { system(load_geo_db_cmd) }
      end

      def load_geo_db_cmd
        geo_db_dir = File.join(db_dir, @params['geo_db'])
        "Rscript #{File.join(GeoDiver.root, 'RCore/download_GEO.R')}" \
        " --accession #{@params['geo_db']}" \
        " --geodbpath #{File.join(geo_db_dir, "#{@params['geo_db']}.soft.gz")}"\
        " --outrdata  #{File.join(geo_db_dir, "#{@params['geo_db']}.RData")}" \
        " && echo 'Finished creating Rdata file:" \
        " #{File.join(geo_db_dir, "#{@params['geo_db']}.Rdata")}'"
      end
    end
  end
end
