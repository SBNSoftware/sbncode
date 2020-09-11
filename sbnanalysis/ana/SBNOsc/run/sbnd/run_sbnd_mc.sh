BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/config
cd /sbnd/data/users/gputnam/NuMu/workshop-game-0320/211a/sbnd/mc
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND.fcl \
	/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/sbnd_files.list \
	> sbnd_mc.out &
cd $BACK
