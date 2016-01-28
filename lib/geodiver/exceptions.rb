# This file defines all possible exceptions that can be thrown by
# GeoDiver on startup.
#
# Exceptions only ever inform another entity (downstream code or users) of an
# issue. Exceptions may or may not be recoverable.
#
# Error classes should be seen as: the error code (class name), human readable
# message (to_s method), and necessary attributes to act on the error.
#
# We define as many error classes as needed to be precise about the issue, thus
# making it easy for downstream code (bin/genevalidatorapp or config.ru) to act
# on them.

module GeoDiver
  # Error in config file.
  class CONFIG_FILE_ERROR < StandardError
    def initialize(ent, err)
      @ent = ent
      @err = err
    end

    attr_reader :ent, :err

    def to_s
      <<MSG
Error reading config file: #{ent}.
#{err}
MSG
    end
  end

  ## NUM THREADS ##

  # Raised if num_threads set by the user is incorrect.
  class NUM_THREADS_INCORRECT < StandardError
    def to_s
      'Number of threads should be a number greater than or equal to 1.'
    end
  end
end
