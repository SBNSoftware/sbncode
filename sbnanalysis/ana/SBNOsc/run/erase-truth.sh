BASE=/sbnd/data/users/gputnam/NuMu/workshop-game-0320/private/outputs
BASEOUT=/sbnd/data/users/gputnam/NuMu/workshop-game-0320/211a
SBND=sbnd/data/
UBOONE=uboone/data/
ICARUS=icarus/data/

INPUTDIRS="$SBND $ICARUS $UBOONE"
#INPUTDIRS="$SBND $UBOONE"
for exp in $INPUTDIRS
do
  inp=$BASE/$exp/211a
  mkdir -p $BASEOUT/$exp
  for file in $inp/*.root
  do
    out=$BASEOUT/$exp/`basename $file`
 
    echo $file
    echo $out
    echo -e .include /sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/"\n" \
          .L ../macro/erase_truth.cc "\n" \
          .L $SBN_LIB_DIR/libsbnanalysis_Event.so "\n" \
          "read_numu_file("\"${file}\", \"${out}\"");" | root -l
  done
done

#echo -e .include /sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/"\n" \
#      .L ../macro/erase_truth.cc "\n" \
#      .L $SBN_LIB_DIR/libsbnanalysis_Event.so "\n" \
#      "read_numu_file("\"${INPUT_FILE}\", \"${OUTPUT_FILE}\"");" | root -l
