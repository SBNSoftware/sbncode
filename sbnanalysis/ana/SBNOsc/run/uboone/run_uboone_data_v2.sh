BACK=$PWD
FILELIST=/sbnd/app/users/gputnam/SBNCode/dev/srcs/sbncode/sbnanalysis/ana/SBNOsc/run_analysis/uboone_ubersim_reweight2_initial_files.list
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/uboone
cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/uboone/data_full/v2
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone11.fcl $FILELIST > uboone11.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone12.fcl $FILELIST > uboone12.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone13.fcl $FILELIST > uboone13.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone14.fcl $FILELIST > uboone14.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone15.fcl $FILELIST > uboone15.out &
cd $BACK
