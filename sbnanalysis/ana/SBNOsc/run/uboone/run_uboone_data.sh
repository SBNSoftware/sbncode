BACK=$PWD
FILELIST=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/uboone_files.list
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/uboone
cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/uboone/data_full/v3
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone1.fcl $FILELIST > uboone1.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone2.fcl $FILELIST > uboone2.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone3.fcl $FILELIST > uboone3.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone4.fcl $FILELIST > uboone4.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone5.fcl $FILELIST > uboone5.out &
cd $BACK
