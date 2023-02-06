BACK=$PWD
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/icarus

model_tags=("6" "7")

for j in ${!model_tags[@]}; do
  mt=${model_tags[$j]}
  FILELIST=/sbnd/app/users/rsjones/LArSoft_v09_10_01/files/geant4_genie_2020a1_icarus.list

  if [ ! -d "/sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/icarus/data_full/v3_2020a1_${mt}" ]; then
    mkdir -p /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/icarus/data_full/v3_2020a1_${mt}
  fi
  cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/icarus/data_full/v3_2020a1_${mt}

  sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernIcarus${mt}.fcl $FILELIST  > icarus${mt}.out &
done
cd $BACK
