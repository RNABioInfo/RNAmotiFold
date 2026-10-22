# Introduction

``RNAmotiFold`` is an RNA secondary structure prediction algorithm that parallely recognizes RNA 3D Motifs through sequence matching. It comes pre-set with a number of motifs build into the algorithm, which are taken from the RNA 3D Motif Atlas. Also provided is a manually extracted set of RNA 3D Motifs from the Rfam database available as well. You can add your own custom motifs as csv files, the required formatting for hairpin, internal and bulge loops is shown in the respective commandline arguments when calling -h.
This repository contains the code for same named Python Package, which makes the algorithm more accessible, implements convenience features, and allows incorporation directly into python scripts. If you want to work with the compiled algorithm binaries, I recommend taking a look at the [RNALoops](https://github.com/RNABioInfo/RNALoops) repo.
The underlying classified dynamic programming algorithms for RNA secondary structure prediction algorithms for RNAmotiFold, RNAmotiCes and RNAmotiAlign are written in Bellman's GAP, translated into C++, and then compiled. The prediction itself is not implemented in python to maintainperformance. You can either clone this repo, or directly install it to your python environment from the python package index with:
`pip install RNAmotiFold`

If you have any questions please reach out to me through my public email address: marius.sebeke@ibmg.uni-stuttgart.de

--- 
# Dependencies
Tested on the following dependencies:
 + Modified Bellman's GAP compiler gapcM == 2026.10.22
 + perl interpreter == perl v5.38.2
 + Python interpreter == Python 3.12.3
 + ubuntu == 24.04

---

# Setup
The RNAmotiFold python wrapper requires the modified Bellman's GAP compiler [gapcM](https://github.com/RNABioInfo/gapcM), along with it's dependencies and a perl interpreter, to be installed. If the modified Bellman's GAP compiler or perl interpreter is not installed system wide it's location can be supplied using the --gapc_path/--perl_path option. Also if your GNU make is not installed system wide it can also be supplied with --make_path. If any of these are not supplied, the script will try to find them. RNAmotiFold is compatible with most UNIX systems, but it is not windows compatible. Unfortunately it is currently not available for Mac Systems, though a port might come in the future. No admin right are required for either the modified Bellman's GAP compiler or RNAmotiFold. These dependencies are only required on first use, they are no longer required after the installation of the algorithm binaries.
Installation Guide:

0. Either directly install RNAmotiFold from PyPI or clone this repository with `git clone https://github.com/RNABioInfo/RNAmotiFold.git --recurse-submodules` and install it into your virtual environment with `pip install .` from the base `RNAmotiFold` directory. 
1. The python wrapper will now be installed, but the compiled RNA secondary structure prediction algorithms are still missing. To set them up, run `RNAmotiFold_installer` on the command line  and supply the dependencies through their associated commandline arguments (e.g. path to your modified GAPC to `--gapcM` ). All three install dependencies are also searched for on your system if they are not supplied, using shutil.which and command -v {dependency}. Installation might take a minute, as all necessaryalgorithms are compiled at once. This process has to be done once whenever installing RNAmotiFold into a new environment, as everything is installed and kept locally.
2. RNAmotiFold is now ready to use in the virtual environment it was installed in, either on the command line by directly calling `RNAmotiFold`, `RNAmotices`, `RNAmoSh`, or `RNAmotiAlign` or in your python scripts. Simply `import RNAmotiFold` and you will be able to call `RNAmotiFold.rnamotifold`, `RNAmotiFold.rnamotices`, `RNAmotiFold.rnamotices_pfc` and all the other algorithms directly from your script. 

---

# Notes

+ Be sure to check out all the available commandline parameters with ``RNAmotiFold -h``!
+ Calling ``RNAmotiFold.py`` without an input will start an interactive session with your set commandline arguments, letting you do individual calculations. Calling it with an input calls your chosen algorithm, does the predictions and returns the output. Interactive sessions are not available through the python API.
+ On the topic of inputs, RNAmotiFold accepts all kinds of inputs, from an RNA sequences as a simple string, to Biopython SeqRecords, to Paths to files or even entire folders. It is also compatible with most of the common file types like .txt, .fa, .sto, and many more.
+ When using RNAmotiFold in a Python script, the return type of the algorithm will be of type AlgorithmOutput which features a method (to_df) to convert your results directly into pandas dataframes! This will always return a list of DataFrames, though it will usually only contain one. The only case where there are multiple DataFrames here is when using any partition function algorithm with the Single Motif Mode. To read more on this refer to our publication: [INSERT RELEASE PUBLICATION HERE]
+ By default, results are piped to stdout and log messages to stderr. You can either pipe them manually with ``>`` and ``2>`` or with the ``RNAmotiFold.py`` commandline parameters ``-o [path to result destination file]`` and ``--logfile [path to log destination file]``
+ Beware `RNAmotiAlign` has not yet been peer reviewed or published and it is still WIP!
+ RNAmotiFold v1.0.5 and later are only compatible with gapcM version 2026.10.22 and later. If you are using RNAmotiFold v1.0.4 or earlier use gapcM version from 2026.09.30 or earlier!
