BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/icarus
cd /sbnd/data/users/gputnam/NuMu/workshop-game/outputs/icarus/data
sbn \
	-m SBNOsc_NumuSelection2 -c $CONFIG_DIR/NumuConfigModernIcarus4.fcl \
	-m SBNOsc_NumuSelection2 -c $CONFIG_DIR/NumuConfigModernIcarus5.fcl \
	$CONFIG_DIR/icarus-files.list  > icarus.out &
cd $BACK
