#!/bin/bash

###########################################################
# Submit sbncode LArSoft jobs to the grid.
#
# Note: Please check the OUTDIR path!
#
# Script arguments are passed to "lar"
#
# From Gray Putnam, taken from script by Andy Mastbaum
###########################################################

OUTDIR=$1
NPROCESSES=$2
FILELIST="sbncode/sbnanalysis/files/$3"
CONFIG=$4

OUTDIR="/pnfs/sbnd/scratch/users/gputnam/numu-selection/${OUTDIR}"
TARDIR="/pnfs/sbnd/resilient/users/gputnam/tars"
SBNCODE_TAR="${TARDIR}/sbncode-selection.tar.gz"

# Log file
LOG="${CLUSTER}_${PROCESS}.log"
echo "Running ${0} on ${HOSTNAME}" >>${LOG} 2>&1
echo "Cluster: ${CLUSTER}" >>${LOG} 2>&1
echo "Process: ${PROCESS}" >>${LOG} 2>&1
echo "sbncode: ${SBNCODE_TAR}" >>${LOG} 2>&1
echo "lar options: $@" >>${LOG} 2>&1
date >>${LOG} 2>&1

# Environment
source /cvmfs/uboone.opensciencegrid.org/products/setup
unsetup mrb
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup mrb
setup larsoft v08_30_02 -q e17:prof

# Get sbncode tarball, unpack, and configure
ifdh cp -D ${SBNCODE_TAR} . >>${LOG} 2>&1
tar xvf $(basename ${SBNCODE_TAR}) >>${LOG} 2>&1

export MRB_PROJECT="larsoft"
export MRB_PROJECT_VERSION="v08_30_02"
export MRB_QUALS="e17:prof"
export MRB_TOP="${PWD}"
export MRB_SOURCE="${MRB_TOP}"
export MRB_INSTALL="${MRB_TOP}"
export PRODUCTS="${MRB_INSTALL}:${PRODUCTS}"

echo "local products: $(ups list -aK+ -z .)" >>${LOG} 2>&1

mrbsetenv >>${LOG} 2>&1
mrbslp >>${LOG} 2>&1

. sbncode/sbnanalysis/bin/setup_sbnanalysis.sh

NFILES=`wc -l < $FILELIST`
echo "Nfiles in list $FILELIST: $NFILES" >>${LOG} 2>&1
NRUN=$((($NFILES + $NPROCESSES - 1) / $NPROCESSES))
FIRST_FILE=$(($NRUN * $PROCESS))
CUT_FILE=$(($NFILES - $FIRST_FILE))
THISFILES=`tail -n $CUT_FILE $FILELIST | head -n $NRUN` 
THISFILES=$(for f in $THISFILES; do pnfsToXRootD $f; done)
echo "Processing files: ${THISFILES}"  >>${LOG} 2>&1
sbn -m SBNOscReco_NumuReco -c $CONFIG ${THISFILES} >>${LOG} 2>&1

LASTROOT=`ls -t | grep root | tail -1`

# Transfer output files to dCache
OUT="${CLUSTER}_${PROCESS}"
mkdir -p ${OUT} >>${LOG} 2>&1
mv $LASTROOT $OUT
mv ${LOG} ${OUT}

ifdh mkdir ${OUTDIR}/${OUT}
ifdh cp -D ${OUT}/* ${OUTDIR}/${OUT}/
ifdh_exit_code=$?

exit ${ifdh_exit_code}
