#!/bin/bash

# if [ "$0" == "$BASH_SOURCE" ]
# then
#     echo 'Please source this script (it needs to modify your environment)'
#     exit 1
# fi

source srt/srt.sh

echo "SRT setup"

export CVMFS_DISTRO_BASE=/cvmfs/nova.opensciencegrid.org/ || exit 1
source setup/setup_nova.sh -b maxopt -6 $SRT_DIST -e $CVMFS_DISTRO_BASE/externals/ || exit 1

cd releases/development/CAFAna/

# Force SL6 here
export SRT_ARCH=Linux2.6
