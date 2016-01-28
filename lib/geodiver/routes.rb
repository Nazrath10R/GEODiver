require 'omniauth'
require 'omniauth-google-oauth2'
require 'sinatra/base'
require 'sinatra/cross_origin'
require 'slim'

require 'geodiver/r_core'
require 'geodiver/version'

module GeoDiver
  # The Sinatra Routes
  class Routes < Sinatra::Base
    register Sinatra::CrossOrigin

    configure do
      # We don't need Rack::MethodOverride. Let's avoid the overhead.
      disable :method_override

      # Ensure exceptions never leak out of the app. Exceptions raised within
      # the app must be handled by the app. We do this by attaching error
      # blocks to exceptions we know how to handle and attaching to Exception
      # as fallback.
      disable :show_exceptions, :raise_errors

      # Make it a policy to dump to 'rack.errors' any exception raised by the
      # app so that error handlers don't have to do it themselves. But for it
      # to always work, Exceptions defined by us should not respond to `code`
      # or http_status` methods. Error blocks errors must explicitly set http
      # status, if needed, by calling `status` method.
      enable :dump_errors

      # We don't want Sinatra do setup any loggers for us. We will use our own.
      set :logging, nil

      # Uer Rack::Session::Pool over Sinatra default sessions as the Pool saves
      # the session info as a instance variable (as compared to within a
      # cookie). Therefore, Rack::Session::Pool should be faster.
      use Rack::Session::Pool

      # Pass OmniAuth the Google Key and Secret Key for Authentication
      use OmniAuth::Builder do
        provider :google_oauth2, ENV['GOOGLE_KEY'], ENV['GOOGLE_SECRET'], {}
      end

      # This is the app root...
      set :root, -> { GeoDiver.root }

      # This is the full path to the public folder...
      set :public_folder, -> { GeoDiver.public_dir }
    end

    get '/' do
      slim :home, layout: false
    end

    get '/analyse' do
      redirect '/auth/google_oauth2' if session[:uid].nil?
      slim :analyse
    end

    get '/history' do
      redirect '/auth/google_oauth2' if session[:uid].nil?
      slim :history, layout: false
    end

    post '/load_geo_db' do
      redirect '/auth/google_oauth2' if session[:uid].nil?
      RCore.init_load_db(params)
      @geo_db_results = RCore.run_load_db
      slim :load_db, layout: false
    end

    post '/analyse' do
      redirect '/auth/google_oauth2' if session[:uid].nil?
      RCore.init_analysis(params, session[:user])
      @results_link = RCore.run_analysis
      slim :results, layout: false
    end

    get '/auth/:provider/callback' do
      content_type 'text/plain'
      session[:uid] = env['omniauth.auth']['uid']
      session[:user] = env['omniauth.auth']
      redirect '/analyse'
    end

    get '/logout' do
      session[:uid] = nil
      session[:user] = nil
      redirect '/'
    end

    get '/auth/failure' do
      content_type 'text/plain'
      request.env['omniauth.auth'].to_hash.inspect
    end

    # This error block will only ever be hit if the user gives us a funny
    # sequence or incorrect advanced parameter. Well, we could hit this block
    # if someone is playing around with our HTTP API too.
    error RCore::ArgumentError do
      status 400
      slim :"500", layout: false
    end

    # This will catch any unhandled error and some very special errors. Ideally
    # we will never hit this block. If we do, there's a bug in GeneValidatorApp
    # or something really weird going on.
    # TODO: If we hit this error block we show the stacktrace to the user
    # requesting them to post the same to our Google Group.
    error Exception, RCore::RuntimeError do
      status 500
      slim :"500", layout: false
    end

    not_found do
      status 404
      slim :"500" # TODO: Create another Template
    end
  end
end
