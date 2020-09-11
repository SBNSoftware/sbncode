BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd
cd /sbnd/data/users/gputnam/NuMu/workshop-game/outputs/sbnd/data
sbn \
	-m SBNOsc_NumuSelection2 -c $CONFIG_DIR/NumuConfigModernSBND4.fcl \
	-m SBNOsc_NumuSelection2 -c $CONFIG_DIR/NumuConfigModernSBND5.fcl \
	$CONFIG_DIR/sbnd-gen.list  > sbnd.out &
cd $BACK
