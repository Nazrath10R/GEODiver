# GeoDiver

[![Build Status](https://travis-ci.org/SBCS-Bioinformatics/GEODiver.svg?branch=master)](https://travis-ci.org/SBCS-Bioinformatics/GEODiver)





## Introduction

GeoDiver is a web app that allows users to easily analyse GEO datasets.







## Installation
### Installation Requirements
* Ruby (>= 2.0.0)


### Installation
Simply run the following command in the terminal.

```bash
gem install geodiver 
```

If that doesn't work, try `sudo gem install geodiver` instead.

##### Running From Source (Not Recommended)
It is also possible to run from source. However, this is not recommended.

```bash
# Clone the repository.
git clone https://github.com/ ...

# Move into GeoDiver source directory.
cd GeoDiver

# Install bundler
gem install bundler

# Use bundler to install dependencies
bundle install

# Optional: run tests and build the gem from source
bundle exec rake

# Run GeneValidator.
bundle exec geodiver -h
# note that `bundle exec` executes GeoDiver in the context of the bundle

# Alternativaly, install GeoDiver as a gem
bundle exec rake install
geodiver -h
```




## Launch GeneValidator

To configure and launch Geodiver, run the following from a command line.

```bash
geodiver
```

Geodiver will automatically guide you through an interactive setup process to help set up your installation.

That's it! Open http://localhost:9292/ and start using GeoDiver!






## Advanced Usage

See `$ geodiver -h` for more information on all the options available when running GeoDiver.

```bash

SUMMARY:
  GeoDiver - A easy to use web tool for analysing GEO datasets.

USAGE:
  $ geodiver [options]

Examples:
  # Launch GeoDiver with the given config file
  $ geodiver --config ~/.geodiver.conf

  # Launch GeoDiver with 8 threads at port 8888
  $ geodiver --num_threads 8 --port 8888


    -c, --config_file        Use the given configuration file
    -g, --gd_public_dir      The public directory that is served to the web application.
    -n, --num_threads        Number of threads to use to run a BLAST search
    -H, --host               Host to run GeoDiver on
    -p, --port               Port to run GeoDiver on (remember to update the port in Google Developer API as well)
    -s, --set                Set configuration value in default or given config file
    -D, --devel              Start GeoDiver in development mode
    -v, --version            Print version number of GeoDiver that will be loaded
    -h, --help               Display this help message.

```


<hr>

This program was developed at [QMUL](http://sbcs.qmul.ac.uk) as part of the [Bioinformatics Masters Course](http://www.qmul.ac.uk/postgraduate/taught/coursefinder/courses/121410.html).
