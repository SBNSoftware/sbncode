BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd
cd /sbnd/data/users/gputnam/NuMu/workshop-game-private/outputs/sbnd/data
sbn 	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND11.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND12.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND13.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND14.fcl \
	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND15.fcl \
	/sbnd/app/users/gputnam/SBNCode/dev/srcs/sbncode/sbnanalysis/ana/SBNOsc/run_analysis/sbnd_ubersim_reweight2_initial_files.list \
	 > sbnd.out &
cd $BACK
