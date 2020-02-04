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

  enum Wall_t {
    kWallNone=0,
    kWallTop=1,
    kWallBottom=2,
    kWallLeft=3,
    kWallRight=4,
    kWallFront=5,
    kWallBack=6
  };

  /// Which generator?
  enum generator_{
    kUnknownGenerator = 0,
    kGENIE            = 1
  };

  /// Enum of possible types of truth-matching a TPC slice
  enum interaction_mode_ {
    kCC = 0,           //!< CC neutrino interaction
    kCCNonPrimary = 1, //!< A non-primary part of a CC neutrino interaction
    kNC = 2,           //!< NC neutrino interaction
    kNCNonPrimary = 3, //!< A non-primary part of an NC neutrino interaction
    kCosmic = 4,       //!< Cosmic activity
    kIntimeCosmic = 5, //!< Cosmic activity in event triggered by intime-cosmic
    mOther = 6 //!< Release valve value -- if nothing else really fits
  };

  /// Which genie status?
  enum genie_status_ : int32_t { 
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
  enum g4_process_ : int32_t {
      kPrimary,
      kCoupledTransportation,
      kFastScintillation,
      kDecay,
      kAntiNeutronInelastic,
      kNeutronInelastic,
      kAntiProtonInelastic,
      kProtonInelastic,
      kHadInelastic,
      kPipInelastic,
      kPimInelastic,
      kXipInelastic,
      kXimInelastic,
      kKaonpInelastic,
      kKaonmInelastic,
      kSigmapInelastic,
      kSigmamInelastic,
      kkaon0LInelastic,
      kKaon0SInelastic,
      kLambdaInelastic,
      kAntiLambdaInelastic,
      kHe3Inelastic,
      kIonInelastic,
      kXi0Inelastic,
      kAlphaInelastic,
      kTInelastic,
      kDInelastic,
      kAntiNeutronElastic,
      kNeutronElastic,
      kAntiProtonElastic,
      kProtonElastic,
      kHadElastic,
      kPipElastic,
      kPimElastic,
      kKaonpElastic,
      kKaonmElastic,
      kKonv,
      kPhot,
      kAnnihil,
      kNCapture,
      kNKiller,
      kMuMinusCaptureAtRest,
      kMuIoni,
      kEBrem,
      kCoulombScat,
      kHBertiniCaptureAtRest,
      kHFritiofCaptureAtRest,
      kPhotonNuclear,
      kPuonNuclear,
      kElectronNuclear,
      kPositronNuclear,
      kCompt,
      kEIoni,
      kMuBrems,
      kHIoni,
      kMuPairProd,
      kHPairProd
  };// g4_process_

}

#endif
