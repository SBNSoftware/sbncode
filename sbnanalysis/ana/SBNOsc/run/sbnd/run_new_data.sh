BACK=$PWD
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd

model_tags=("minus50pctBE" "plus50pctBE")

for j in ${!model_tags[@]}; do
  mt=${model_tags[$j]}
  FILELIST=/sbnd/app/users/rsjones/LArSoft_v09_10_01/files/geant4_genie_${mt}_nieves_active_sbnd.list

  if [ ! -d "/sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/v3_${mt}_nieves" ]; then
    mkdir -p /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/v3_${mt}_nieves
  fi
  cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/v3_${mt}_nieves

  sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND1.fcl $FILELIST  > sbnd1.out &
done
cd $BACK
