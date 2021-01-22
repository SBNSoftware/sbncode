#ifndef SRENUMS_H
#define SRENUMS_H

namespace caf
{
  /// Which SBN detector?
  enum Det_t
  {
    kUNKNOWN,     ///< Unknown detector
    kSBND,        ///< Near Detector
    kICARUS      ///< Far Detector
  };

  enum Plane_t 
  {
    kUnknown=-1,
    k1stInduction=0,
    k2ndInduction=1,
    kCollection=2
  };

  enum Wall_t 
  {
    kWallNone=0,
    kWallTop=1,
    kWallBottom=2,
    kWallLeft=3,
    kWallRight=4,
    kWallFront=5,
    kWallBack=6
  };

  /// Which type of MC?
  enum MCType_t
  {
    kMCUnknown=0,
    kMCParticleGun=1,
    kMCNeutrino=2,
    kMCCosmic=3,
    kMCOverlay=4
  };

  /// Which generator?
  enum generator_
  {
    kUnknownGenerator = 0,
    kGENIE            = 1
  };

  /// Enum of possible types of truth-matching a TPC slice
  enum interaction_mode_ 
  {
    kCC = 0,           //!< CC neutrino interaction
    kNC = 1,           //!< NC neutrino interaction
    kCosmic = 2,       //!< Cosmic activity
    kIntimeCosmic = 3, //!< Cosmic activity in event triggered by intime-cosmic
    kOther = 4 //!< Release valve value -- if nothing else really fits
  };

  /// Which genie status?
  enum genie_status_ 
  { 
      kIStUndefined                = -1,
      kIStInitialState             =  0,
      kIStStableFinalState         =  1,
      kIStIntermediateState        =  2,
      kIStDecayedState             =  3,
      kIStCorrelatedNucleon        = 10,
      kIStNucleonTarget            = 11,
      kIStDISPreFragmHadronicState = 12,
      kIStPreDecayResonantState    = 13,
      kIStHadronInTheNucleus       = 14,
      kIStFinalStateNuclearRemnant = 15,
      kIStNucleonClusterTarget     = 16,
      kNotGenie                    = 17 //!< Not a genie particle
  };//genie_status

  /// Which G4 process ?
  enum g4_process_ 
  {
    kG4primary=0,
    kG4CoupledTransportation=1,
    kG4FastScintillation=2,
    kG4Decay=3,
    kG4anti_neutronInelastic=4,
    kG4neutronInelastic=5,
    kG4anti_protonInelastic=6,
    kG4protonInelastic=7,
    kG4hadInelastic=8,
    kG4pipInelastic=9,
    kG4pimInelastic=10,
    kG4xipInelastic=11,
    kG4ximInelastic=12,
    kG4kaonpInelastic=13,
    kG4kaonmInelastic=14,
    kG4sigmapInelastic=15,
    kG4sigmamInelastic=16,
    kG4kaon0LInelastic=17,
    kG4kaon0SInelastic=18,
    kG4lambdaInelastic=19,
    kG4anti_lambdaInelastic=20,
    kG4He3Inelastic=21,
    kG4ionInelastic=22,
    kG4xi0Inelastic=23,
    kG4alphaInelastic=24,
    kG4tInelastic=25,
    kG4dInelastic=26,
    kG4anti_neutronElastic=27,
    kG4neutronElastic=28,
    kG4anti_protonElastic=29,
    kG4protonElastic=30,
    kG4hadElastic=31,
    kG4pipElastic=32,
    kG4pimElastic=33,
    kG4kaonpElastic=34,
    kG4kaonmElastic=35,
    kG4conv=36,
    kG4phot=37,
    kG4annihil=38,
    kG4nCapture=39,
    kG4nKiller=40,
    kG4muMinusCaptureAtRest=41,
    kG4muIoni=42,
    kG4eBrem=43,
    kG4CoulombScat=44,
    kG4hBertiniCaptureAtRest=45,
    kG4hFritiofCaptureAtRest=46,
    kG4photonNuclear=47,
    kG4muonNuclear=48,
    kG4electronNuclear=49,
    kG4positronNuclear=50,
    kG4compt=51,
    kG4eIoni=52,
    kG4muBrems=53,
    kG4hIoni=54,
    kG4muPairProd=55,
    kG4hPairProd=56,
    kG4LArVoxelReadoutScoringProcess=57
  };// g4_process_

}

#endif
