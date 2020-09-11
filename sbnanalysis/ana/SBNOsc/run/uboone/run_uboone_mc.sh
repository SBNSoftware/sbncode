BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/uboone
cd /sbnd/data/users/gputnam/NuMu/workshop-game-0320/211a/uboone/mc
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernUboone.fcl \
	/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/uboone_files.list \
	> uboone_mc.out &
cd $BACK
