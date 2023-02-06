BACK=$PWD
FILELIST=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/icarus_files.list
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/icarus
cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/icarus/data_full/v3
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus0.fcl $FILELIST  > icarus0.out &
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus1.fcl $FILELIST  > icarus1.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus2.fcl $FILELIST  > icarus2.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus3.fcl $FILELIST  > icarus3.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus4.fcl $FILELIST  > icarus4.out &
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus5.fcl $FILELIST  > icarus5.out &
cd $BACK
