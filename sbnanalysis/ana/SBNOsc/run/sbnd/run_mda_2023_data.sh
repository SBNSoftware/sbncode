BACK=$PWD
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd

model_tags=("smith")

for j in ${!model_tags[@]}; do
  mt=${model_tags[$j]}
  FILELIST=/sbnd/app/users/rsjones/LArSoft_v09_66_00/files/${mt}_geant4_files_mda_2023.list

  if [ ! -d "/sbnd/data/users/rsjones/NuMu/mda_2023/sbnd/data_full/v3p2_${mt}" ]; then
    RUN_DIR=/sbnd/data/users/rsjones/NuMu/mda_2023/sbnd/data_full/v3p2_${mt}
    mkdir -p ${RUN_DIR}
  fi
  cd ${RUN_DIR}

  sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND1.fcl $FILELIST  > sbnd1.out &
done
cd $BACK
