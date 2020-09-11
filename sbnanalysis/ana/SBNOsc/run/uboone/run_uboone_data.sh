BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/uboone
FILELIST=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/uboone_files.list
cd /sbnd/data/users/gputnam/NuMu/workshop-game-0320/private/outputs/uboone/data/211a
sbn 	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone1.fcl $FILELIST > uboone1.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone2.fcl $FILELIST > uboone2.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone3.fcl $FILELIST > uboone3.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone4.fcl $FILELIST > uboone4.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone5.fcl $FILELIST > uboone5.out &
cd $BACK
