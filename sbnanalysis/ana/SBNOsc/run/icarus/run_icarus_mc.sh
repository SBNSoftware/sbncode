BACK=$PWD
CONFIG_DIR=/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/icarus
cd /sbnd/data/users/gputnam/NuMu/workshop-game-0320/211a/icarus/mc
sbn	-m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus.fcl \
	/sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/ana/SBNOsc/files/icarus_files.list \
	> icarus_mc.out &
cd $BACK
