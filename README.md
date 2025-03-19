# Introduction

``RNAmotiFold`` is an RNA secondary structure prediction algorithm that parallely recognizes RNA 3D Motifs through sequence matching. It comes pre-set with a number of motifs build into the algorithm. These are taken from the RNA 3D Motif Atlas, which is updated semi regularly, and a manually extracted set of RNA 3D Motifs from the Rfam database. You can add your own custom motifs as csv files, the required formatting for hairpin, internal and bulge loops is shown when calling -h in the respective commandline arguments. 
This package comes with a python wrapper (``RNAmotiFold.py``) for the underlying algorithms implemented in Bellman's GAP. If you want to directly work with the compiled algorithm binaries they will be located under `/RNAmotiFold/Build/bin` after installation.

For further information you read: [INSERT FINISHED PAPER HERE]
If you have any questions please reach out to me through my public email address: marius.sebeke@ibmg.uni-stuttgart.de

--- 

# Dependencies
Tested on the following dependencies:
 + CMake ==  3.16
 + m4 == 1.4.19
 + perl interpreter == perl v5.38.2
 + C++ compiler ==  GCC g++ 13.2.0
 + C compiler == GCC 13.2.0
 + GNU make ==  3.81

 + Python packages:
    + requests == 2.32.3
    + bio == 1.7.1
 ### Optional dependencies (setup installs them locally if you don't already have them installed):
 + GNU bison == 3.81
 + Flex == 2.5.34
 + GNU scientific library (GSL) == 2.7
 + boost == 1.86
    + unit-test framework
    + pool
    + program options
    + system
    + filesystem
    + serialization

---

# Setup
RNAmotiFold should be compatible with most UNIX system, but it is not windows compatible.
You do not need admin rights for your system to install it, everything will be installed locally.



0. Clone this repository with `git clone https://github.com/RNABioInfo/RNAmotiFold.git --recurse-submodules`
1. Run `pip install -r requirements.txt` to install non-default python packages into your python environment.
2. Run `python3 setup.py` (this may take a couple minutes).
3. You can now use ``RNAmotiFold.py`` as a wrapper for ``RNAmotiFold``, ``RNAmoSh`` and ``RNAmotiCes`` or run any of them directly from ``/RNAmotifold/Build/bin/``.

### Usage
1. Be sure to check out all the available commandline parameters with ``python RNAmotiFold.py -h``!
2. Calling ``RNAmotiFold.py`` without an input ( ``-i [input]``) will start an interactive session with your set commandline arguments, letting you do individual calculations. Calling it with an input calls your chosen algorithm, does the predictions and returns the output.

---

# Notes

+ By default, results are piped to stdout and log messages to stderr. You can either pipe them manually with ``>`` and ``2>`` or with the ``RNAmotiFold.py`` commandline parameters ``-o [path to result destination file]`` and ``--logfile [path to log destination file]``
+ Beware that custom motif identifiers can NOT be used with ``RNAmoSh``, unless you add the specific letters to the shape alphabet of the gapcM, which you can find under /RNAmotiFold/Build/gapc-prefix/src/gapcM/rtlib/shape_alph.hh. After changing the file you will need to recompile the gapcM by re-running the cMake build script (`cd` into the Build folder and run `cmake .. && cmake --build .` ).
    + On the other hand you can also just reuse the already implemented motif abbreviations: G, U, T, L, D, A, M, K, S, C (and their lower case versions) for you own sequences. Only do this when fully replacing the Rfam/RNA 3D Motif Atlas sequences with your own motif sequence set to avoid ambiguities (use the ``-X [path to hairpin loop csv]``/``-Y [path to internal loop csv]``/``-Z [path to bulge loop csv]`` and ``-L``/``-E``/``-G`` commandline parameters)!
+ The ``RNAmotiFold.py`` wrapper takes arguments from the commandline or from a config.ini file, a premade config file is provided under ``RNAmotiFold/src/config.ini``. In case you want to change the defaults I set for all the parameters, you can change them by editing the ``/RNAmotiFold/src/data/defaults.ini``.
