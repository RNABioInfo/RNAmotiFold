# Introduction

``RNAmotiFold`` is an RNA secondary structure prediction algorithm that parallely recognizes RNA 3D Motifs through sequence matching. It comes pre-set with a number of motifs build into the algorithm, which are taken from the RNA 3D Motif Atlas. Also provided is a manually extracted set of RNA 3D Motifs from the Rfam database available as well. You can add your own custom motifs as csv files, the required formatting for hairpin, internal and bulge loops is shown in the respective commandline arguments when calling -h.
This repository contains the code for same named Python Package, which I wrote to make the algorithm more accessible and to implement it into python scripts.
The underlying algorithms for RNAmotiFold, RNAmotiCes and RNAmotiAlign are written in Bellman's GAP, translated into C++, and then compiled. The prediction is thus separate from python, keeping it as efficient as possible. You can either clone this repo, or directly install it from PyPI.
This package comes with a python wrapper (``RNAmotiFold.py``) for the underlying algorithms implemented in Bellman's GAP. If you want to directly work with the compiled algorithm binaries, I recommend taking a look at the [RNALoops](https://github.com/RNABioInfo/RNALoops) repo. These lack several features though, including the single motif mode, which is implemented as part of this python package.
If you have any questions please reach out to me through my public email address: marius.sebeke@ibmg.uni-stuttgart.de

--- 
# Dependencies
Tested on the following dependencies:
 + Modified Bellman's GAP compiler gapcM == 2026.09.26
 + perl interpreter == perl v5.38.2
 + GNU make ==  3.81
 + Python interpreter == Python 3.12.3
 + ubuntu == 22.04

---

# Setup
RNAmotiFold comes with one (technically three) runtime dependencies: The modified Bellman's GAP compiler [gapcM](https://github.com/RNABioInfo/gapcM), a perl interpreter and GNU make, both of which are dependencies for the compiler as well. Your python interpreter has to be at least version 3.12 to use this wrapper. RNAmotiFold should be compatible with most UNIX systems, but it is not windows compatible. Unfortunately it is currently not available for Mac Systems, though a port might come in the future. No admin right are required for either the modified Bellman's GAP compiler or RNAmotiFold.
Installation Guide:

0. Either directly install RNAmotiFold from PyPI or clone this repository with `git clone https://github.com/RNABioInfo/RNAmotiFold.git --recurse-submodules` and install it into your virtual environment with `pip install .` from the base `RNAmotiFold` directory. 
1. The python wrapper will now be installed, but the compiled algorithms are still missing. To set them up, run any of the installed commandline interfaces like `RNAmotiFold` or `RNAmotiCes` and supply the dependencies through their associated commandline arguments (e.g. path to your modified GAPC to `--gapc_path` ). All three install dependencies are also searched for on your system if they are not supplied, using shutil.which and command -v {dependency}. Installation might take a second as all the binaries are installed during this process, but it only has be done once (unless updates happen).
2. RNAmotiFold is now ready to use in the virtual environment it was installed in, either through the command line by directly calling `RNAmotiFold`, `RNAmotices`, `RNAmoSh`, or `RNAmotiAlign` or in your scripts. Simply `import RNAmotiFold` and you will be able to call `RNAmotiFold.rnamotifold`, `RNAmotiFold.rnamotices`, `RNAmotiFold.rnamotices_pfc` and all the other algorithms directly from your script. 

---

# Notes

+ Be sure to check out all the available commandline parameters with ``RNAmotiFold -h``!
+ Calling ``RNAmotiFold.py`` without an input will start an interactive session with your set commandline arguments, letting you do individual calculations. Calling it with an input calls your chosen algorithm, does the predictions and returns the output. Interactive sessions are not available from within python scripts.
+ On the topic of inputs, RNAmotiFold accepts all kinds of inputs, a RNA sequences as a simple string, to Biopython SeqRecords, to Paths to files or even entire folders. It is also compatible with most of the common file types like .txt, .fa, .sto, and many more.
+ When using RNAmotiFold in a Python script, the return type of the algorithm will be of type AlgorithmOutput which features a method (to_df) to convert your results directly into pandas dataframes! This will always return a list of DataFrames, though it will usually only contain one. The only case where there are multiple DataFrames here is when using any partition function algorithm with the Single Motif Mode. To read more on this refer to our publication: [INSERT RELEASE PUBLICATION HERE]
+ By default, results are piped to stdout and log messages to stderr. You can either pipe them manually with ``>`` and ``2>`` or with the ``RNAmotiFold.py`` commandline parameters ``-o [path to result destination file]`` and ``--logfile [path to log destination file]``
+ Beware `RNAmotiAlign` has not yet been peer reviewed or published and it and its parameters are still WIP!
