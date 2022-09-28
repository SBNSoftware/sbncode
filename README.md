SBN Code
========
This package contains code for SBN-wide simulation and analysis. This includes:

* `fcl`: Production FHiCL configurations
* `sbncode`: Common LArSoft modules and other utilities
* `env`: Scripts to set up joint SBND/MicroBooNE/ICARUS environments (to be
  automated in future versions)

Cloning submodules
------------------
`mrb gitCheckout` should take care of git submodules if the origin repository is the
default one, i.e. SBNSoftware.
If `sbncode` is cloned from any repository other than the official one (SBNSoftware),
git submodules need to be initalized after `mrb gitCheckout` with the following command:

```
$ git submodule update --init --recursive
```
