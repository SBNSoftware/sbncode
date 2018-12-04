#!/bin/bash

###############################################################################
# Set up the environment and install sbncode on a GPVM
#
# A. Mastbaum <mastbaum@uchicago.edu>, 2018/06/12
###############################################################################

LAR_VERSION="v07_11_00"
LAR_QUAL="e17:prof"

# Set up the environment
TOPDIR=${PWD}
cat << EOF > ${TOPDIR}/setup_sbncode-${LAR_VERSION}.sh
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode ${LAR_VERSION} -q ${LAR_QUAL}
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode ${LAR_VERSION} -q ${LAR_QUAL}
unsetup larbatch
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup icaruscode ${LAR_VERSION} -q ${LAR_QUAL}
EOF

source setup_sbncode-${LAR_VERSION}.sh

# Create working area
WDIR="${PWD}/sbncode-${LAR_VERSION}"
mkdir $WDIR

if [ $? -ne 0 ]; then
  echo "Unable to create directory ${WDIR}"
  exit
fi

cd ${WDIR}
mrb newDev -v ${LAR_VERSION} -q ${LAR_QUAL}
LAR_QUAL_STR=${LAR_QUAL//:/_}
MRB_SETUP=${WDIR}/localProducts_larsoft_${LAR_VERSION}_${LAR_QUAL_STR}/setup
source ${MRB_SETUP}

# Get sbncode
cd srcs
mrb gitCheckout -t ${LAR_VERSION} sbncode

# Build sbncode
mrbsetenv
mrb i -j4

cat << EOF >> ${TOPDIR}/setup_sbncode-${LAR_VERSION}.sh
source ${MRB_SETUP}
mrbslp
mrbsetenv
EOF

# Build sbnanalysis
cd ${WDIR}/srcs/sbncode/sbnanalysis
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make install -j4
source bin/setup_sbnanalysis.sh

echo "source ${PWD}/bin/setup_sbnanalysis.sh" >> ${TOPDIR}/setup_sbncode-${LAR_VERSION}.sh

echo -e "Setup finished. To use this environment, run\n\n  source setup_sbncode-${LAR_VERSION}.sh"

