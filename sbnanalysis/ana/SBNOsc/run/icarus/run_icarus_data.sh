BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/icarus
FILELIST=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/icarus_files.list
cd /sbnd/data/users/gputnam/NuMu/workshop-game-0320/private/outputs/icarus/data/211a
sbn 	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus0.fcl $FILELIST  > icarus0.out &
sbn 	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus1.fcl $FILELIST  > icarus1.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus2.fcl $FILELIST  > icarus2.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus3.fcl $FILELIST  > icarus3.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus4.fcl $FILELIST  > icarus4.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus5.fcl $FILELIST  > icarus5.out &
cd $BACK
