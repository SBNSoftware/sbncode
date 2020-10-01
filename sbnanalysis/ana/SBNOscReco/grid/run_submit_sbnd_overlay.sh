NPROCESS=100
jobsub_submit --memory=4GB --group=sbnd --OS=SL6 --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE,DEDICATED -N $NPROCESS file:///sbnd/app/users/gputnam/sbndcode/grid-sbn/submit_combined.sh sbnd-overlay3 $NPROCESS overlay_files.list NumuRecoSBND.fcl 
