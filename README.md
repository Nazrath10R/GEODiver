# GeoDiver

[![Build Status](https://travis-ci.org/SBCS-Bioinformatics/GEODiver.svg?branch=master)](https://travis-ci.org/SBCS-Bioinformatics/GEODiver)





## Introduction

GeoDiver is a web app that allows users to easily analyse GEO datasets.







## Installation
### Installation Requirements
* Ruby (>= 2.0.0)
  * Recommended to use rvm to install ruby
* R (=3.2.2)
  * Recommended to use R to install R



### Installation
Simply run the following command in the terminal.

```bash
# Clone the repository.
git clone https://github.com/SBCS-Bioinformatics/GEODiver

# Move into GeoDiver source directory.
cd GEODiver

# Install R dependencies & Build and install the latest version of the webapp.
rake install 

# Start the web app
passenger start --envvar GOOGLE_KEY=113114282317-af1ph7hqm7uvhbc289gpu5fteuo8i4a3.apps.googleusercontent.com --envvar GOOGLE_SECRET=_QrM-_WoTNjcreNRAw6MXfZE -p 9292 -e production
```

##### Running From Source (Not Recommended)
It is also possible to run from source. However, this is not recommended.

```bash
# After cloning the web app and moving into the source directory 
# Install bundler
gem install bundler

# Use bundler to install dependencies
bundle install

# Optional: run tests and build the gem from source
bundle exec rake

# Run GeoDiver
bundle exec passenger start -h
# note that `bundle exec` executes GeoDiver in the context of the bundle

# Alternatively run Geodiver using the command line interface
bundle exec geodiver -h
```




## Launch GeoDiver

To configure and launch Geodiver, run the following from a command line from the GeoDiver root folder.

```bash
bundle exec passenger start -h

```
That's it! Open http://localhost:9292/ and start using GeoDiver!






## Advanced Usage

See `$ passenger start -h` for more information on all the options available when running GeoDiver.

# Config file
A Config file can be used to specify arguments - the default location of this file is in the home directory at `~/.geodiver.conf`. An examplar of the config file can be seen below.


```yaml
---
:num_threads: 8
:port: '9292'
:host: 0.0.0.0
:gd_public_dir: "/Users/ismailm/.geodiver"
:devel: true
```


<hr>

This program was developed at [QMUL](http://sbcs.qmul.ac.uk) as part of the [Bioinformatics Masters Course](http://www.qmul.ac.uk/postgraduate/taught/coursefinder/courses/121410.html).
