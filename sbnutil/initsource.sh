unsetup ifdhc
setup ifdhc              v2_4_1           -q e17:p2714b:prof

echo "Setting up SBNanalysis" 
export SBN_ANALYSIS_DIR=$MRB_INSTALL/sbncode/sbnanalysis/build
export SBN_LIB_DIR=$SBN_ANALYSIS_DIR/lib
export LD_LIBRARY_PATH=$SBN_LIB_DIR:$LD_LIBRARY_PATH
export PATH=$SBN_ANALYSIS_DIR/bin:$PATH
export FHICL_FILE_PATH=.:$MRB_INSTALL/sbncode/sbnanalysis/build/fcl:$FHICL_FILE_PATH
function run_sbn_grid_SBNOsc_NueSelection() {
    sbn-grid -m SBNOsc_NueSelection $@
}
