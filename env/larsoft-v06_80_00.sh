source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup icaruscode v06_80_00 -q e15:prof
unsetup cetpkgsupport
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode v06_80_00 -q e15:prof
unsetup mrb
unsetup git
unsetup gitflow
unsetup cetpkgsupport
unsetup pycurl
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v06_80_00 -q e15:prof

