BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/uboone
cd /sbnd/data/users/gputnam/NuMu/workshop-game-private/outputs/uboone/data
sbn 	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone11.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone12.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone13.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone14.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone15.fcl \
	/sbnd/app/users/gputnam/SBNCode/dev/srcs/sbncode/sbnanalysis/ana/SBNOsc/run_analysis/uboone_ubersim_reweight2_initial_files.list \
	$CONFIG_DIR/uboone-files.list  > uboone.out &
cd $BACK
