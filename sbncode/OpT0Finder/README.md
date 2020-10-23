# OpT0Finder

This repository contains a flash-matching code in the directory `flashmatch`. The `flashmatch`
directory should be kept up-to-date with the one in https://github.com/drinkingkazu/OpT0Finder.
This software can be run either inside LArSoft, or in a standalone mode. This is dictaded by the
compiler directive `USE_LARSOFT=1` or `=0`.

The LArSoft configuration file is in the `job` directory, while the LArSoft plugins are supposed
to be added to the experiment's respository. See, for example, directory `OpT0Finder` in `sbndcode`.