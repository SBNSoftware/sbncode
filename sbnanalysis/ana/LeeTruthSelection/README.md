Truth-Based 1l1p Selection
==========================
A. Mastbaum <mastbaum@uchicago.edu>, 2017/08/13

Updated: G. Putnam <grayputnam@uchicago.edu> on 2018/01/22

Loop through MCTracks and MCShowers, building a list of tracks on
showers using truth variables which can be configured to be distorted,
then selected based on that list.  

Selection code is built as a standalone executable which uses Gallery as
a library. Scans through Monte Carlo events and produces an output root
file which stores both truth level information and the distorted
information. Usage information below.

Code to generate histograms and covariance plots from this output is
build as a gallery-framework module and is run using a python script.
Further code to generate plots live in a number of python scripts. Usage
information below.

Tools
-----
### TSSelection: Truth-Based Selection ###
This module is the main truth-based selection code. It scans through Monte
Carlo events to extract MCTracks and MCShowers, distorts the truth
information based on configuration, applies selection cuts 
(e.g. kinetic energy), and classifies events as 1eip and 1muip. 
The exact signal regions (0p/1p/Np/Ntrack) can also be configured.
The selected events are written into an output analysis tree.

Note on Selection: I have recently been updating the way that selection
writes data to the output TTree's, and I have added in some more
configuration options. I haven't been able to run the code recently
however, and so these most recent features are untested, and may be a
little buggy. I'll remove this section when I've had a chance to test
stuff.

### TSProcessor: Connects to SBNAnalysis ###

TSProcessor inherits from `SelectionBase` and connects stuff up with
`SBNAnalysis`.

**Selection Usage:**

See README for `SBNAnalysis`


*Configuration Arguments*

Arguments are set w/ a config json file. See source for arguments.

