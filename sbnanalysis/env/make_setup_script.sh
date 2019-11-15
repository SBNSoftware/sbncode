#!/bin/bash

###############################################################################
# Set up the environment and install sbncode on a GPVM
#
# A. Mastbaum <mastbaum@uchicago.edu>, 2018/06/12
###############################################################################

EOFD="EOF"

LAR_VERSION="v08_30_00"
LAR_QUAL="e17:prof"

# Check we're in a good location
#ALLOWED="/sbnd/app/users/*"
#if [[ ${PWD} != *"${ALLOWED}"* ]]; then
#  echo "Please move to a location inside appp."
#  exit
#fi

# Set up the environment
TOPDIR=${PWD}
cat << EOF > ${TOPDIR}/setup_sbncode-${LAR_VERSION}.sh
echo "Your Experiment is ${HOSTNAME%%gpvm*}"

unsetup cetpkgsupport

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode ${LAR_VERSION} -q ${LAR_QUAL}
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode ${LAR_VERSION} -q ${LAR_QUAL}
unsetup larbatch
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup icaruscode ${LAR_VERSION} -q ${LAR_QUAL}

if [ "${HOSTNAME%%gpvm*}"  == "sbnd" ]; then
    source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
elif [ "${HOSTNAME%%gpvm*}"  == "uboone" ]; then
    source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
elif [ "${HOSTNAME%%gpvm*}"  == "icarus" ]; then
    source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
else
    echo "Sorry. Are you not on a gpvm?"
    exit
fi

unsetup cetpkgsupport
setup cetpkgsupport

EOF

source ${TOPDIR}/setup_sbncode-${LAR_VERSION}.sh



cat << EOF >> ${TOPDIR}/setup_sbncode-${LAR_VERSION}.sh

setup larbatch

source $MRB_INSTALL/setup
mrbslp
mrbsetenv
mrbslp
source $MRB_SOURCE/sbncode/sbnanalysis/build/bin/setup_sbnanalysis.sh
source $MRB_SOURCE/sbncode/sbnanalysis/bash_completion 
unsetup larbatch
setup larbatch 

#For project.py 
echo "Setting up grid settings"

#Add the the experiment_utilities to the python path so all three detectors are set up.
PYTHONPATH=/sbnd/app/users/dbarker/larsoft_sbncode/sbncode/srcs/sbncode/sbnutil/python:$PYTHONPATH


#Create the setup script for the grid nodes.
sbndcode_version=$(grep "sbndcode" $MRB_SOURCE/sbncode/ups/product_deps | head -1 | cut -c11- | xargs) 
icaruscode_version=$(grep "icaruscode" $MRB_SOURCE/sbncode/ups/product_deps | head -1 | cut -c11- | xargs)
uboonecode_version=$(grep "uboonecode" $MRB_SOURCE/sbncode/ups/product_deps | head -1 | cut -c11- | xargs)

quals=$MRB_QUALS

cat << EOF >  $MRB_SOURCE/sbncode/sbnutil/setup_sbncode_grid.sh

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
echo "setting up uboonecode version: $uboonecode_version"
setup uboonecode \$uboonecode_version -q \$quals
unsetup larbatch

source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh    
echo "setting up sbndcode version: $sbndcode_version"
setup sbndcode \$sbndcode_version -q \$quals
unsetup larbatch

source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
echo "setting up icaruscode version: $icaruscode_version"
setup icaruscode  \$icaruscode_version -q \$quals 

if [ "sbnd"  == "sbnd" ]; then
    source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
elif [ "sbnd"  == "uboone" ]; then
    source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
elif [ "sbnd"  == "icarus" ]; then
    source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
else
    echo "Sorry. Are you not on a gpvm?"
    exit
fi

$EOFD

#Setup the init source
ifdhc_version=$(ups active | grep ifdhc | head -1 |cut -c 6-34 |xargs)               
ifdhc_qual=$(ups active | grep ifdhc | head -1 | sed -n 's/.* -q \([^ ]*\).*/\1/p' |xargs)

cat << EOF > $MRB_SOURCE/sbncode/sbnutil/initsource.sh
unsetup ifdhc
setup ifdhc \$ifdhc_version -q \$ifdhc_qual

echo "Setting up SBNanalysis" 
export SBN_ANALYSIS_DIR=\\\$MRB_INSTALL/sbncode/sbnanalysis/build
export SBN_LIB_DIR= \\\$SBN_ANALYSIS_DIR/lib
export LD_LIBRARY_PATH=\\\$SBN_LIB_DIR:\\\$LD_LIBRARY_PATH
export PATH=\\\$SBN_ANALYSIS_DIR/bin:\\\$PATH
export FHICL_FILE_PATH=.:\\\$MRB_INSTALL/sbncode/sbnanalysis/build/fcl:\\\$FHICL_FILE_PATH
function run_sbn_grid_SBNOsc_NueSelection() {
    sbn-grid -m SBNOsc_NueSelection \\\$@
}
$EOFD


make_sbnanalysis_tar(){

    rm -rf \$MRB_INSTALL/sbncode/sbnanalysis/
    mkdir \$MRB_INSTALL/sbncode/sbnanalysis/
    mkdir \$MRB_INSTALL/sbncode/sbnanalysis/build 
    mkdir \$MRB_INSTALL/sbncode/sbnanalysis/core
    cp    \$SBN_ANALYSIS_DIR/../core/Event.hh $MRB_INSTALL/sbncode/sbnanalysis/core
    cp    \$SBN_ANALYSIS_DIR/../core/SubRun.hh $MRB_INSTALL/sbncode/sbnanalysis/core
    cp -r \$SBN_ANALYSIS_DIR/*  $MRB_INSTALL/sbncode/sbnanalysis/build
    tar -C \$MRB_INSTALL -czf sbncode_sbnanalysis_$LARSOFT_VERSION.tar .
}



EOF


echo "source \$MRB_SOURCE/sbncode/sbnanalysis/build/bin/setup_sbnanalysis.sh" >> ${TOPDIR}/setup_sbncode-${LAR_VERSION}.sh


