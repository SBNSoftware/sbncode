BACK=$PWD
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd
FILELIST=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/sbnd_files.list
cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/v3
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND1.fcl $FILELIST  > sbnd1.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND2.fcl $FILELIST  > sbnd2.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND3.fcl $FILELIST  > sbnd3.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND4.fcl $FILELIST  > sbnd4.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND5.fcl $FILELIST  > sbnd5.out &
cd $BACK
