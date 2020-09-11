BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd
FILELIST=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/sbnd_files.list
cd /sbnd/data/users/gputnam/NuMu/workshop-game-0320/private/outputs/sbnd/data/211a
sbn 	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND1.fcl $FILELIST  > sbnd1.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND2.fcl $FILELIST  > sbnd2.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND3.fcl $FILELIST  > sbnd3.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND4.fcl $FILELIST  > sbnd4.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND5.fcl $FILELIST  > sbnd5.out &
cd $BACK
