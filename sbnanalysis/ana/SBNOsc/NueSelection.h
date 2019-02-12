#ifndef __sbnanalysis_ana_SBNOsc_NueSelection__
#define __sbnanalysis_ana_SBNOsc_NueSelection__

/**
 * \file NueSelection.h
 *
 * SBN nue selection.
 *
 * Author:
 */

//C++ Includes 
#include <iostream>
#include <vector> 

//Framework Includes 
#include "canvas/Utilities/InputTag.h"

//SBN Includes 
#include "core/SelectionBase.hh"
#include "core/Event.hh"

//Larsoft Includes 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"

// take the geobox stuff from uboonecode                          
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"

//Root Includes
#include "TRandom3.h"


class TH2D;

namespace ana {
  namespace SBNOsc {

/**
 * \class NueSelection
 * \brief Electron neutrino event selection
 */
class NueSelection : public core::SelectionBase {
public:
  /** Constructor. */
  NueSelection();

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL);

  /** Finalize and write objects to the output file. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco);

protected:

  unsigned EventCounter;  //!< Count processed events
  unsigned NuCount;  //!< Count selected events

  /** Configuration parameters */
  struct Config {
    //    art::InputTag fTruthTag;  //!< art tag for MCTruth information
    std::vector<geoalgo::AABox> fiducial_volumes; //!< List of FV containers -- set by "fiducial_volumes"
    std::vector<geoalgo::AABox> active_volumes; //!< List of active volumes
    double minLengthExitingTrack; //!< Minimum length [cm] of exiting tracks.  Will not apply cut if value is negative.
    double trackVisibleEnergyThreshold; //!< Energy threshold for track to be acounted in visible energy calculation [GeV].
    double ShowerVisibleEnergyThreshold; //!< Energy at which the showers are visible in the detector.
    double showerEnergyDistortion; //!< Energy distortion of showers for visible energy calculation (%).
    double trackEnergyDistortion; //!< Energy distortion of tracks for visible energy calculation (%).
    double leptonEnergyDistortionContained; //<! Energy distortion of lepton (primary track) for visible energy calculation (%).
    double leptonEnergyDistortionLeavingA; //!< Parameter to be used in the energy distortion of primary lepton for visible energy calculation. 
    // (%) = -A * Log(B * L)  where L is the lepton contained length
    double leptonEnergyDistortionLeavingB; //!< Parameter to be used in the energy distortion of primary lepton for visible energy calculation. 
    double nueEfficencyWeight; 
    double photonEnergyCut; //!< Energy cut at which photons can be identified. 
    double VtxEnergyCut; //!< Hadronic Energy cut on the vertext to identify photons from there containment length.
    double PhotonConvLenghCut; //!< Photon conversion length cut.  
    double dEdxPhotonCut; //!< dEdx cut from photons and electrons.

    bool ApplyElectronEnergyCut;
    bool NueIDEfficiency;
    bool doFVCut;
    bool verbose; 
  };

    /** Additional information used by the selection per neutrino interaction */
  struct NueInteraction {
    double hadronic_energy;
    double shower_energy;
    double leptonic_energy;
    double nue_energy; 

    NueInteraction(){
      hadronic_energy = -99999;
      shower_energy   = -99999;
      leptonic_energy = -99999;
      nue_energy      = -99999;
    }
  };
  


  Config fConfig; //!< The config
  
  TRandom *randomnum_gen = new TRandom3();

  bool Select(const gallery::Event& ev, const simb::MCTruth& mctruth, std::vector<const simb::MCTrack>& mctracks,unsigned truth_ind,const Config fConfig,const std::vector<double> visible_energy, std::map<int, const simb::MCParticle*>& mcparticles, NueInteraction& intInfo);

  //Selection Functions 
  bool passFV(const TVector3 &v) {return !fConfig.doFVCut || containedInFV(v);}
  bool passNueIDEfficiencyCut(const simb::MCNeutrino& nu);
  bool passeEnergyCut(const simb::MCNeutrino& nu){return !fConfig.ApplyElectronEnergyCut || eEnergyCut(nu);}
  bool passPhotonEnergyCut(std::vector<int> photons, std::map<int, const simb::MCParticle*>& mcparticles, int photonTrackID);
  bool passConversionGapCut(int photonTrackID, std::map<int, const simb::MCParticle*>& mcparticles);
  bool passMuLengthCut(std::vector<const sim::MCTrack>& mctrack_list);

  bool containedInFV(const TVector3 &v);
  bool eEnergyCut(const simb::MCNeutrino& nu);
  std::vector<int> findNeutralPions(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth);
  std::vector<int> findPhotons(std::vector<int>& pi_zeros,std::map<int, const simb::MCParticle*>& mcparticles);

};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelection__

