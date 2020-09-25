BACK=$PWD
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd
FILELIST=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/sbnd_files.list
cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND1.fcl $FILELIST  > sbnd1.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND2.fcl $FILELIST  > sbnd2.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND3.fcl $FILELIST  > sbnd3.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND4.fcl $FILELIST  > sbnd4.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND5.fcl $FILELIST  > sbnd5.out &
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND11.fcl $FILELIST  > sbnd11.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND12.fcl $FILELIST  > sbnd12.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND13.fcl $FILELIST  > sbnd13.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND14.fcl $FILELIST  > sbnd14.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND15.fcl $FILELIST  > sbnd15.out &
cd $BACK
