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
#include <ctime> 

//Framework Includes 
#include "canvas/Utilities/InputTag.h"

//SBN Includes 
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/PostProcessorBase.hh"

//Larsoft Includes 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"

// take the geobox stuff from uboonecode                          
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"

//Root Includes
#include "TRandom3.h"
#include "THStack.h"

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
    double showerVisibleEnergyThreshold; //!< Energy at which the showers are visible in the detector.
    double showerEnergyDistortion; //!< Energy distortion of showers for visible energy calculation (%).
    double trackEnergyDistortion; //!< Energy distortion of tracks for visible energy calculation (%).
    double leptonEnergyDistortionContained; //<! Energy distortion of lepton (primary track) for visible e
    double nueEfficency; 
    double vtxEnergyCut; //!< Hadronic Energy cut on the vertext to identify photons from there containment length.
    double photonConvLenghCut; //!< Photon conversion length cut.  
    double dEdxPhotonCut; //!< dEdx cut from photons and electrons.
    double GlobalWeight; //Global weight to accoutn for e.g. change in size for SBND weight. 
    double POTWeight; 
    double photonVisibleEnergyThreshold; //Secondary photons can be seen at 100 MeV
    double CosmicWeight;
    double CosmicVolumeRadius;

    std::string Detector;  

    std::vector<std::string> UniformWeights;

    bool ApplyKMECCut; //Wasn't in the proposal remove?  
    bool ApplyNueEfficiency; 
    bool ApplyElectronEnergyCut;
    bool ApplyFVCut;
    bool ApplyAVCut;
    bool ApplyPhotonEnergyCut;
    bool ApplyConversionGapCut;
    bool ApplydEdxCut;
    bool ApplyMuonLenghtCut;
    bool ApplyCosmicCylinderCut;
    bool ApplyCosmicFVCut;

    bool IncludeCosmics;
    bool IncludeDirt;
    bool CosmicsOnly;
    bool DirtOnly;
    bool Verbose; 


  };

  //Additional information used by the selection per neutrino interaction/ 
  //C++ is smart and doesn't need the constructor.
  struct NueInteraction {
    double hadronic_energy;
    double shower_energy;
    double leptonic_energy;
    double weight;
    int initnu;
    int leptontrackID;
    int leptonpdgID;
    int fnd;
    double GetNueEnergy(){return hadronic_energy + shower_energy + leptonic_energy;}
  };
  
  //Output Root Histogram map for background selections
  struct RootHistograms{

    double emin = 0.0, emax = 3.0;
    int ebins = 120;

    std::map<std::string,TH1D*> TrueEnergyAll_Hist;
    std::map<std::string,TH1D*> TrueEnergy_Hist;
    std::map<std::string,TH1D*> CCQEEnergy_Hist;

    std::map<std::string,TH1D*> VisibleEnergy_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_AVCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_FVCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_EnergyCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PhotonEnergyCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ConversionGapCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MuLenghtCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_NCCut_Hist;
   
    std::map<std::string,TH1D*> VisibleEnergy_LeptonPlusPhotonCut_Hist;

    TH1D* VisibleEnergy_CosmicFVCut_Hist;
    TH1D* VisibleEnergy_CosmicClyinderCut_Hist;
    TH1D* VisibleEnergy_CosmicdEdxCut_Hist;
    TH1D* VisibleEnergy_CosmicWeightCut_Hist;
    TH1D* VisibleEnergy_CosmicEnergyCut_Hist;

    TH1D* VisibleEnergy_FNKPDecay_Hist;
    TH1D* VisibleEnergy_FNKMDecay_Hist;
    TH1D* VisibleEnergy_MuDecays_Hist;

    //Final selection histograms
    TH1D* VisibleEnergy_CosmicSelection_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_Selection_Hist;
    std::vector<std::string> HistTypes = {"NuMu","InNuE","OscNuE","NCInNuE","NCOscNuE","NCNuMu","DirtNuMu","DirtInNuE","DirtOscNuE","DirtNCInNuE","DirtNCOscNuE","DirtNCNuMu"};

  };

  Config fConfig; //!< The config
  RootHistograms fRootHists;

  TRandom *randomnum_gen = new TRandom3();

  bool Select(const gallery::Event& ev, const simb::MCTruth& mctruth, const std::vector<sim::MCTrack>& mctracks, unsigned truth_ind, std::map<int, const simb::MCParticle*>& mcparticles, NueInteraction& intInfo);
  bool SelectCosmic(const sim::MCShower& mcs, std::map<int, const simb::MCParticle*>& mcparticles, NueSelection::NueInteraction& intInfo);

  //Selection Functions 
  bool passFV(const TVector3 &v) {return fConfig.ApplyFVCut && containedInFV(v);}
  bool passAV(const TVector3 &v) {return fConfig.ApplyAVCut && containedInAV(v);}
  bool passNueIDEfficiencyCut(const simb::MCNeutrino& nu){return fConfig.ApplyNueEfficiency && NueIDEfficiencyCut(nu);}
  bool passeEnergyCut(double& lepton_energy){return fConfig.ApplyElectronEnergyCut && eEnergyCut(lepton_energy);}
  bool passPhotonEnergyCut(std::vector<int>& photons, std::map<int, const simb::MCParticle*>& mcparticles, int& photonTrackID){return fConfig.ApplyPhotonEnergyCut && PhotonEnergyCut(photons,mcparticles,photonTrackID);}
  bool passMuLengthCut(const std::vector<sim::MCTrack>& mctrack_list, const simb::MCTruth& mctruth){return fConfig.ApplyMuonLenghtCut && MuLengthCut(mctrack_list, mctruth);}
  bool passdEdxCut(int pdgcode){return fConfig.ApplydEdxCut && dEdxCut(pdgcode);}
  bool passConversionGapCut(std::map<int, const simb::MCParticle*>& mcparticles, int photonTrackID, const float& hadronic_energy, const simb::MCNeutrino& nu){return fConfig.ApplyConversionGapCut && ConversionGapCut(mcparticles,photonTrackID,hadronic_energy,nu);}
  bool PassCosmicInFV(const sim::MCShower& mcs, TVector3& vertex){return fConfig.ApplyCosmicFVCut || CosmicInFV(mcs, vertex);}
  bool PassCosmicCylinderCut(const sim::MCShower& mcs, TVector3& vertex, std::map<int, const simb::MCParticle*>& mcparticles){return fConfig.ApplyCosmicCylinderCut || CosmicCylinderCut(mcs, vertex, mcparticles);}


  //Cut Functions
  bool containedInFV(const TVector3 &v);
  bool containedInAV(const TVector3 &v);
  bool NueIDEfficiencyCut(const simb::MCNeutrino& nu);
  bool eEnergyCut(double leptonenergy);
  bool PhotonEnergyCut(std::vector<int>& photons, std::map<int, const simb::MCParticle*>& mcparticles, int &photonTrackID);
  bool ConversionGapCut(std::map<int, const simb::MCParticle*>& mcparticles, int photonTrackID, const float& hadronicE, const simb::MCNeutrino& nu);
  bool dEdxCut(int pdgcode);
  bool MuLengthCut(const std::vector<sim::MCTrack>& mctrack_list, const simb::MCTruth& mctruth);
  bool CosmicInFV(const sim::MCShower& mcs, TVector3& vertex);
  bool CosmicCylinderCut(const sim::MCShower& mcs, TVector3& vertex, std::map<int, const simb::MCParticle*>& mcparticles);

  std::vector<int> findNeutralPions(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth);
  std::vector<int> findPhotons(std::vector<int>& pi_zeros,std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth);

  void InitialiseHistograms();
  void FillHistograms(std::map<std::string,TH1D*>& HistMap, const simb::MCNeutrino& nu, NueSelection::NueInteraction& intInfo, bool& booldirtevent);
  void FillHistograms(std::map<std::string,TH1D*>& HistMap, const simb::MCNeutrino& nu, NueSelection::NueInteraction& intInfo,double Energy, bool &booldirtevent);

  void PrintInformation(const simb::MCTruth& mctruth, NueSelection::NueInteraction& intInfo);

  std::vector<int> ElectronNuDecays = {1,2,6,9};
  std::vector<int> MuonNuDecays = {3,4,5,7,8,10,11,12,13,14};

  //Measure CPU times 
  std::clock_t c_start;
  std::clock_t c_end;

  TRandom rand;
  double fPOT;
  TH1I*  fEventHist;
  bool dirtevent;
  int MCTruthCounter; 
  
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelection__

