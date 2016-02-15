require 'forwardable'
require "base64"

# GeoDiver Namespace
module GeoDiver
  # Class to create the history
  class History
    class << self
      extend Forwardable

      def_delegators GeoDiver, :config, :logger, :public_dir, :users_dir,
                     :db_dir

      def run(user)
        generate_history(user)
      end

      def generate_history(user)
        user_dir = Pathname.new(users_dir) + Pathname.new(user.info['email'])
        return [] unless user_dir.exist?
        data = []
        user_dir.children(with_directory=false).each do |accession|
          next unless accession.to_s =~ /^GDS/
          (user_dir + accession).children(with_directory=false).each do |time|
            next unless time.to_s.length == 33
            data << generate_data_hash(user_dir, user, accession, time)
          end
        end
        data.sort_by { |d| d[:time] }.reverse
      end

      def generate_data_hash(user_dir, user, accession, time)
        json_file = user_dir + accession + time + 'params.json'
        data = { geo_db: accession.to_s,
          user: encode_email(user),
          uniq_result: time.to_s,
          uniq_url: generate_uniq_url(user, accession, time),
          share_url: generate_share_url(user, accession, time),
          full_path: (user_dir + accession + time).to_s,
          relative_path: (Pathname.new('../GeoDiver/Users') +
                          Pathname.new(user.info['email']) + accession + time).to_s,
          time: Time.strptime(time.to_s, '%Y-%m-%d_%H-%M-%S_%L-%N'),
          params_json: json_file.to_s,
          params: JSON.parse(IO.read(json_file.to_s))
        }
        data[:share] = true if File.exist? ( user_dir + accession + time + '.share' )
        data
      end

      def generate_uniq_url(user, accession, time)
        'result/' + encode_email(user) + '/' + accession.to_s + '/' + time.to_s
      end

      def generate_share_url(user, accession, time)
        'sh/' + encode_email(user) + '/' + accession.to_s + '/' + time.to_s
      end

      def encode_email(user)
        Base64.encode64(user.info['email']).chomp
      end
    end
  end
end
