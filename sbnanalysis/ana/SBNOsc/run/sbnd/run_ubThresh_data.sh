BACK=$PWD
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd

model_tags=("uBooNEThresh")

for j in ${!model_tags[@]}; do
  mt=${model_tags[$j]}
  FILELIST=/sbnd/app/users/rsjones/LArSoft_v09_10_01/files/geant4_genie_nieves_nieves_active_sbnd.list

  if [ ! -d "/sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/v3_nieves_${mt}" ]; then
    mkdir -p /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/v3_nieves_${mt}
  fi
  cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/v3_nieves_${mt}

  sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND${mt}.fcl $FILELIST  > sbnd${mt}.out &
done
cd $BACK
