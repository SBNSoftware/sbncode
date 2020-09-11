BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/icarus
cd /sbnd/data/users/gputnam/NuMu/workshop-game-private/outputs/icarus/data
sbn 	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus11.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus12.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus13.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus14.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus15.fcl \
	/sbnd/app/users/gputnam/SBNCode/dev/srcs/sbncode/sbnanalysis/ana/SBNOsc/run_analysis/icarus_ubersim_reweight2_initial_files.list \
	$CONFIG_DIR/icarus-files.list  > icarus.out &
cd $BACK
