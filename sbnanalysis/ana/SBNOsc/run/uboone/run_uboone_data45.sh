BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/uboone
cd /sbnd/data/users/gputnam/NuMu/workshop-game/outputs/uboone/data
sbn \
	-m SBNOsc_NumuSelection2 -c $CONFIG_DIR/NumuConfigModernUboone4.fcl \
	-m SBNOsc_NumuSelection2 -c $CONFIG_DIR/NumuConfigModernUboone5.fcl \
	$CONFIG_DIR/uboone-files.list  > uboone.out &
cd $BACK
