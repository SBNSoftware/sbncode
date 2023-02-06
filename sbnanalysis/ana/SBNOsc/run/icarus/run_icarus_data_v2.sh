BACK=$PWD
FILELIST=/sbnd/app/users/gputnam/SBNCode/dev/srcs/sbncode/sbnanalysis/ana/SBNOsc/run_analysis/icarus_ubersim_reweight2_initial_files.list
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/icarus
cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/icarus/data_full/v2
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus11.fcl $FILELIST  > icarus11.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus12.fcl $FILELIST  > icarus12.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus13.fcl $FILELIST  > icarus13.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus14.fcl $FILELIST  > icarus14.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus15.fcl $FILELIST  > icarus15.out &
cd $BACK
