base=/sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a

detectors=("sbnd" "icarus")
models=("susav2_nieves" "smith_nieves" "lse_nieves")

#models
for i in ${!models[@]}; do
  m=${models[$i]}

  # detectors
  for j in ${!detectors[@]}; do
    d=${detectors[$j]}

    input=$base/${d}/data_full/v3_${m}
    output=$base/${d}/data/v3_${m}

    if [ ! -d $output ]; then
      mkdir -p $output
    fi

    for file in $input/*.root
    do
      outfile=$output/`basename $file`

      echo $file
      echo $outfile
      echo -e .include /sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/"\n" \
        .L ../macro/erase_truth.cc "\n" \
        .L $SBN_LIB_DIR/libsbnanalysis_Event.so "\n" \
        "read_numu_file("\"${file}\", \"${outfile}\"");" | root -l
    done
  done
done

#echo -e .include /sbnd/app/users/gputnam/sbncode-gen/srcs/sbncode/sbnanalysis/"\n" \
#      .L ../macro/erase_truth.cc "\n" \
#      .L $SBN_LIB_DIR/libsbnanalysis_Event.so "\n" \
#      "read_numu_file("\"${INPUT_FILE}\", \"${OUTPUT_FILE}\"");" | root -l
