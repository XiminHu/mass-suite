# mass-suite
[![Build Status](https://travis-ci.com/XiminHu/mass-suite.svg?branch=master)](https://travis-ci.com/XiminHu/mass-suite)
[![Coverage Status](https://coveralls.io/repos/github/XiminHu/mass-suite/badge.svg?branch=master)](https://coveralls.io/github/XiminHu/mass-suite?branch=master)

This package is initiated from the University of Washington eScience capstone project.

## Mass-suite: Comprehensive Mass Spectrometry data analysis tool

Mass-suite is a package used for general MS data process and advanced data analysis.

It provided a full capacity including data import, alignment, noise removal, end member extraction, visualization and fragment searching using online database.

#### Contributers: Ximin Hu, Derek Mar, Nozomi Suzuki, Bowei Zhang
#### Release date: 

## Installation & major dependencies

## Organization of the project
The project has the following structure:
   
   
    mass-suite/
      |- README.md
      |- mss/
         |- __init__.py
         |- visreader.py
         |- tests/
            |- __init__.py
            |- test.py
         |- dev/  
            |- *.ipynb
      |- doc/
      |- example_data/
      |- LICENSE

## Featured modules

<font color = 'red'>(waiting for updates..)</font>

## Project Data

All the data used in the project is giving credit to [Center of Urban Water](https://www.urbanwaters.org/).

## Model training and testing

All the data used for model training is under 'example_data' folder and can be repeatly trained by users according to the settings.

## Testing and continuous integration

The majority of code is tested under unittest and pushed through [travis CI](https://travis-ci.com/github/XiminHu/mass-suite).

## Licensing

The package is open source and can be utilized under MIT license. Please find the detail in `licence` file.