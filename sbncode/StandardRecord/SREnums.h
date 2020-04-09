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
    kG4primary,
    kG4CoupledTransportation,
    kG4FastScintillation,
    kG4Decay,
    kG4anti_neutronInelastic,
    kG4neutronInelastic,
    kG4anti_protonInelastic,
    kG4protonInelastic,
    kG4hadInelastic,
    kG4pipInelastic,
    kG4pimInelastic,
    kG4xipInelastic,
    kG4ximInelastic,
    kG4kaonpInelastic,
    kG4kaonmInelastic,
    kG4sigmapInelastic,
    kG4sigmamInelastic,
    kG4kaon0LInelastic,
    kG4kaon0SInelastic,
    kG4lambdaInelastic,
    kG4anti_lambdaInelastic,
    kG4He3Inelastic,
    kG4ionInelastic,
    kG4xi0Inelastic,
    kG4alphaInelastic,
    kG4tInelastic,
    kG4dInelastic,
    kG4anti_neutronElastic,
    kG4neutronElastic,
    kG4anti_protonElastic,
    kG4protonElastic,
    kG4hadElastic,
    kG4pipElastic,
    kG4pimElastic,
    kG4kaonpElastic,
    kG4kaonmElastic,
    kG4conv,
    kG4phot,
    kG4annihil,
    kG4nCapture,
    kG4nKiller,
    kG4muMinusCaptureAtRest,
    kG4muIoni,
    kG4eBrem,
    kG4CoulombScat,
    kG4hBertiniCaptureAtRest,
    kG4hFritiofCaptureAtRest,
    kG4photonNuclear,
    kG4muonNuclear,
    kG4electronNuclear,
    kG4positronNuclear,
    kG4compt,
    kG4eIoni,
    kG4muBrems,
    kG4hIoni,
    kG4muPairProd,
    kG4hPairProd
  };// g4_process_

}

#endif
