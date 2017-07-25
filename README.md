# Introduction
Performs correlation analysis of BurningMan Census and MOOP map data.

# Installation
Dependencies (other than C++ compiler):
* CMake
* ROOT (www.root.cern.ch)
* libSILLY
* LibXml2
* LaTeX (optional)

# Data
Data is not included. 

The MOOP maps can be downloaded here:
* [2016](http://o7fe62guj6g73vlj30xpogpm.wpengine.netdna-cdn.com/wp-content/uploads/2016/10/Moop-Map_Day-9_Placement-Updated_hi-res.jpg)

Census data only available to Census team.

Quick start guide after all dependencies are met and the repository has been cloned.
```
$ mkdir build
$ cd build
$ cmake ../
$ make
```

All executables and libraries will appear in the `build` directory.

## Code structure

A quick guid to the code structure:
* README.md -> this help
* CMakeLists.txt -> input file for cmake
* inc/ -> header files
* src/ -> source files for libraries
* util/ -> source files for executables
* macros/ -> ROOT macros for final results analysis (makes nice summary plots)
* cmake/ -> contains extra CMake modules, if needed
* data/ -> contains data (e.g. MOOP maps, ...)
* doc/ -> contains more detailed documentation, especially on the statistical treatment 

# Analysis description

This code include tools to perform various steps of the analysis, in particular to:
* Digitize and validate MOOP map data
* Test correlation extraction algorithm on fake (toy) data (also used to estimate statistical uncertainties)
* Analyze data and extract correlations (only binary categories supported so far)

## Digitize MOOP map data
The information on each MOOP map is summarized in an XML file (see `data/MoopMap-2016.xml` for an example).
The actual data is then loaded from an image that is referenced in the XML file itself.
Two reference points are needed to translate the image positions (pixels) in an absolute scale (radius and angle, see inc/MoopMap.h for a description of the reference frame).

Various checks can be performed on the MOOP map to validate it. The executable `checkMap` contains:
* interactive navigation of the map
* calculation of MOOP areas
* random-sampling to re-create the map and check its global consistency

Just execute `checkMap -h` for more help.

## Toy data
Toy data can be generated using the `generateToys` executable. Invoking `generateToys -h` will provide more help.

## Analyze data
To analyze toy data it is convenient to use the `analyzeToyData` executable.
The macro `macros/drawToyResults.C` can be helpful to plot the response for various input asymmetries and determine the sensitivity as well as the expected statistical uncertainty of a given estimation.

To analyze real data, it first has to be transformed into a ROOT readable file, using `parseCensusCSV.R`.
Then data can be analyzed using the `analyzeCensusData` executable.

