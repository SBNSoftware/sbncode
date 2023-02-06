BACK=$PWD
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd
FILELIST=/sbnd/app/users/gputnam/SBNCode/dev/srcs/sbncode/sbnanalysis/ana/SBNOsc/run_analysis/sbnd_ubersim_reweight2_initial_files.list
cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/v2
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND11.fcl $FILELIST  > sbnd11.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND12.fcl $FILELIST  > sbnd12.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND13.fcl $FILELIST  > sbnd13.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND14.fcl $FILELIST  > sbnd14.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND15.fcl $FILELIST  > sbnd15.out &
cd $BACK
