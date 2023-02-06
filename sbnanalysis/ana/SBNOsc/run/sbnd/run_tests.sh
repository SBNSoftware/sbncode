BACK=$PWD
CONFIG_DIR=/sbnd/app/users/rsjones/SBNCode_Repository/sbncode-v08_47_00/srcs/sbncode/sbnanalysis/ana/SBNOsc/run/sbnd
GEVFILE=/pnfs/sbnd/persistent/users/rsjones/GENIE_G18_10a_02_11a_Driver_1M/g4/39642948_130/cryoTestPOT1M_GEN_G4-20201119T172241.root
LARFILE=/pnfs/sbnd/persistent/users/gputnam/numu_workshop_game_0320/211a/samples/sbnd/prodgenie_bnb_nu_sbn_sbnd_20200722T100022_gen_20200722T101604_g4_20200722T101840_bnbcorrection_20200722T105934_eventweightA_20200722T130322_eventweightB.root

cd /sbnd/data/users/rsjones/NuMu/workshop-game-0320/211a/sbnd/data_full/test

sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND_GEV.fcl $GEVFILE  > sbnd_gev.out &
sbn -m SBNOsc_NumuSelection -c $CONFIG_DIR/NumuConfigModernSBND_LAr.fcl $LARFILE  > sbnd_lar.out &

cd $BACK
