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

      def_delegators GeoDiver, :config, :logger, :public_dir, :users_dir,
                     :db_dir

      #
      def init(params)
        logger.debug('Loading Database')
        assert_params(params)
      end

      #
      def run
        meta_json = File.join(db_dir, @params['geo_db'],
                              "#{@params['geo_db']}.json")
        if File.exist?(meta_json)
          meta_data = JSON.parse(IO.read(meta_json))
        else
          meta_data = download_geo_meta_data(meta_json)
        end
        soft_link_meta_json_to_public_dir(meta_json)
        logger.debug("Meta Data: #{meta_data}")
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
        return unless @params['geo_db'].nil? || @params['geo_db'].empty?
        fail ArgumentError, 'No GEO database provided.'
      end

      def download_geo_meta_data(output_json)
        file    = download_geo_file
        content = read_geo_file(file)
        data    = parse_geo_db(content)
        write_to_json(data, output_json)
        data
      end

      def download_geo_file
        remote_dir = generate_remote_url
        output_dir = File.join(db_dir, @params['geo_db'])
        FileUtils.mkdir(output_dir) unless Dir.exist? output_dir
        output_file = File.join(output_dir, "#{@params['geo_db']}.soft.gz")
        system "wget #{remote_dir} --output-document #{output_file}"
        system "gunzip --force #{output_file}" # force to overwrite file
        output_file.gsub(/.gz$/, '')
      end

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

      def read_geo_file(file)
        content = []
        IO.foreach(file) do |line|
          break if line =~ /^#ID_REF/
          content << line
        end
        content.join
      end

      def parse_geo_db(content)
        {
          "Accession": content.match(/\^DATASET = (.*)/)[1],
          "Title": content.match(/!dataset_title = (.*)/)[1],
          "Description": content.match(/!dataset_description = (.*)/)[1],
          "Sample_Organism": content.match(/!dataset_platform_organism = (.*)/)[1],
          "Factors": parse_factors(content),
          "Reference": content.match(/!Database_ref = (.*)/)[1],
          "Update_Date": content.match(/!dataset_update_date = (.*)/)[1]
        }
      end

      def parse_factors(content)
        subsets = content.gsub(/\^DATA.*\n/, '').gsub(/\![dD]ata.*\n/, '')
        results = {}
        subsets.lines.each_slice(5) do |subset|
          description = subset[2].match(/\!subset_description = (.*)/)[1]
          type = subset[4].match(/\!subset_type = (.*)/)[1].gsub(' ', '.')
          # samples = subset[3].match(/\!subset_sample_id = (.*)/)[1]
          results[type] ||= []
          results[type] << description
        end
        results
      end

      def write_to_json(hash, output_json)
        File.open(output_json, 'w') { |f| f.puts hash.to_json }
      end

      #
      def soft_link_meta_json_to_public_dir(meta_json)
        public_meta_json = File.join(public_dir, 'GeoDiver/DBs/',
                                     "#{@params['geo_db']}.json")
        return if File.exist? public_meta_json
        FileUtils.ln_s(meta_json, public_meta_json)
      end
    end
  end
end
