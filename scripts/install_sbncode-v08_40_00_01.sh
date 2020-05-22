#!/bin/bash

###############################################################################
# Set up the environment and install sbncode on a GPVM
#
# A. Mastbaum <mastbaum@uchicago.edu>, 2018/06/12
###############################################################################

# Check we're in a good location
#ALLOWED="/sbnd/app/users/${USER}"
#if [[ ${PWD} != *"${ALLOWED}"* ]]; then
#  echo "Please move to a location inside ${ALLOWED}."
#  exit
#fi

if [ $# -eq 0 ]; then
  LAR_VERSION="v08_40_00"
  LAR_QUAL="e19:prof"
  SBN_VERSION="v08_40_00_01"
else
  LAR_VERSION=$1
  LAR_QUAL=$2
  SBN_VERSION=$LAR_VERSION
fi

if [ -z "$LAR_VERSION" ] || [ -z "$LAR_QUAL" ]; then
  echo "Invalid LArSoft version: $LAR_VERSION -q $LAR_QUAL"
  exit 1
else
  echo "Using LArSoft version: $LAR_VERSION -q $LAR_QUAL"
fi

# Set up the environment
TOPDIR=${PWD}
cat << EOF > ${TOPDIR}/setup_sbncode-${SBN_VERSION}.sh
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode ${LAR_VERSION} -q ${LAR_QUAL}
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode ${LAR_VERSION} -q ${LAR_QUAL}
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup icaruscode ${LAR_VERSION} -q ${LAR_QUAL}
EOF

source ./setup_sbncode-${SBN_VERSION}.sh

# Create working area
WDIR="${PWD}/sbncode-${SBN_VERSION}"
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
mrb gitCheckout -t ${SBN_VERSION} sbncode

#mrb gitCheckout -t master sbncode
#echo "Bailing early..."
#exit

# Build sbncode
mrbsetenv
mrb i -j4

cat << EOF >> ${TOPDIR}/setup_sbncode-${SBN_VERSION}.sh
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

echo "source ${PWD}/bin/setup_sbnanalysis.sh" >> ${TOPDIR}/setup_sbncode-${SBN_VERSION}.sh

echo -e "Setup finished. To use this environment, run\n\n  source setup_sbncode-${SBN_VERSION}.sh"

