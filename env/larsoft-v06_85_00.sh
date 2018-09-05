LAR_VERSION="v06_85_00"
LAR_QUAL="e15:prof"
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode ${LAR_VERSION} -q ${LAR_QUAL}
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode ${LAR_VERSION} -q ${LAR_QUAL}
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup icaruscode ${LAR_VERSION} -q ${LAR_QUAL}
