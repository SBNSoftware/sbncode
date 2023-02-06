BASE=/sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a

dets=("SBND" "ICARUS")
locations=("sbnd" "icarus")
fd_scales=(0.832 0.860 0.878)
nd_scales=(0.739 0.758 0.779)

models=("SuSAv2" "SmithMoniz" "LlewellynSmith-Empirical")
model_tags=("susav2" "smith" "lse")

for i in ${!dets[@]}; do
  d=${dets[$i]}
  l=${locations[$i]}

  for j in ${!models[@]}; do
    m=${models[$j]}
    mt=${model_tags[$j]}

    if [ d == "SBND" ]; then
      scale=${nd_scales[$j]}
    elif [ d == "ICARUS" ]; then
      scale=${fd_scales[$j]}
    fi

    inp=${l}/data_full/v3_${mt}_nieves/

    for file in $inp/*.root
    do
      outfile=$out/`basename $file`
   
      echo $file
      echo $outfile
      echo -e .include /sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/"\n" \
            .L ../macro/setpotweight_macro.cc "\n" \
            .L $SBN_LIB_DIR/libsbnanalysis_Event.so "\n" \
            "setpotweight("\"${file}\", ${scale}");" | root -l
  
    done
  done
done

