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
#include <map>

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
#include "TH2D.h"

#include "Utilities.h"

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
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco);

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
    double CosmicGlobalWeight;
    double CosmicVolumeRadius;
    double ElectronWeight; 
    double PhotonWeight;

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
    bool ApplyCosmicInSpillWindowCut;
    bool ApplyGENIEVersionWeight;


    bool IgnoreNeutrinoDepsInCosmic;

    bool IncludeCosmics;
    bool IncludeDirt;
    bool CosmicsOnly;
    bool DirtOnly;
    bool Verbose;
    bool IntrinsicOnly;
    bool NuMuOnly;
    bool OscOnly;  
    bool FillHistograms;
    bool FillModeHistograms;
    bool FillIntTypeHistograms;
    bool UseAllCosmics;
    bool DontUseSimChannels;
    bool UseGenieHists;

    float TrackEnergyCorrection;
    float ShowerEnergyCorrection;
    float NutrinoEnergyCorrection;

    float GlobalTimeOffset;
    float SpillTime;
    float ReadoutWindowSize;
  };

  //Additional information used by the selection per neutrino interaction/ 
  //C++ is smart and doesn't need the constructor.
  struct NueInteraction {
    double hadronic_energy = 0;
    double shower_energy   = 0;
    double leptonic_energy = 0;
    double other_energy    = 0; 
    double weight;
    int initnu;
    int leptontrackID;
    int leptonpdgID;
    int fnd;
    double GetNueEnergy(){return hadronic_energy + shower_energy + leptonic_energy + other_energy;}

    double GetRecoEnergy(){return (hadronic_energy+shower_energy);}
  };
  
  //Output Root Histogram map for background selections
  struct RootHistograms{

    double emin = 0.0, emax = 10;
    int ebins = 400;

    std::map<std::string,TH1D*> TrueNumber_Hist;
    std::map<std::string,TH1D*> TrueEnergyAll_Hist;
    std::map<std::string,TH1D*> TrueEnergy_Hist;
    std::map<std::string,TH1D*> CCQEEnergy_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_AVCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_FVCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_FVBefore_Hist;

    std::map<std::string,TH1D*> VisibleEnergy_EnergyCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PhotonEnergyCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ConversionGapCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MuLenghtCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_NCCut_Hist;
    std::map<std::string,TH1D*> Weights_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_LeptonPlusPhotonCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PiZero_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_Photon_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PhotonSmall_Hist;
    
     std::map<std::string,TH1D*> VisibleEnergy_NoShowerCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_TwoPhotonCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_TwoShowerCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PhotonFVCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_TwoLeptonsAfter_Hist;

     std::map<std::string,TH1D*> VisibleEnergy_NoShowerCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_TwoPhotonCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_TwoShowerCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PhotonFVCut_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_TwoLeptons_Hist;


    std::map<std::string,TH1D*> VisibleEnergyAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_AVCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_FVCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_FVBeforeAfter_Hist;

    std::map<std::string,TH1D*> VisibleEnergy_EnergyCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PhotonEnergyCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ConversionGapCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MuLenghtCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_NCCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_LeptonPlusPhotonCutAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PiZeroAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PhotonAfter_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_PhotonSmallAfter_Hist;


    std::map<std::string,TH1D*> LowNCEnergy_Hist;

    std::map<std::string,TH1D*> ProtonE_Hist;
    std::map<std::string,TH1D*> PionE_Hist;
    std::map<std::string,TH1D*> KaonE_Hist;
    std::map<std::string,TH1D*> ProtonN_Hist;
    std::map<std::string,TH1D*> PionN_Hist;
    std::map<std::string,TH1D*> KaonN_Hist;
    std::map<std::string,TH1D*> HadronE_Hist;
    std::map<std::string,TH1D*> PhotonCon_Hist;
    

    std::map<std::string,std::map<int,TH1D*> > TrueNumber_HistMode;
    std::map<std::string,std::map<int,TH1D*> > TrueEnergyAll_HistMode;
    std::map<std::string,std::map<int,TH1D*> > TrueEnergy_HistMode;
    std::map<std::string,std::map<int,TH1D*> > CCQEEnergy_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_AVCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_FVCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_FVBefore_HistMode;

    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_EnergyCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PhotonEnergyCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_ConversionGapCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_MuLenghtCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_NCCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_LeptonPlusPhotonCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> > Weights_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PiZero_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_Photon_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PhotonSmall_HistMode;
    std::map<std::string,std::map<int,TH1D*> > LowNCEnergy_HistMode;
    std::map<std::string,std::map<int,TH1D*> > ProtonE_HistMode;
    std::map<std::string,std::map<int,TH1D*> > PionE_HistMode;
    std::map<std::string,std::map<int,TH1D*> > KaonE_HistMode;
    std::map<std::string,std::map<int,TH1D*> > ProtonN_HistMode;
    std::map<std::string,std::map<int,TH1D*> > PionN_HistMode;
    std::map<std::string,std::map<int,TH1D*> > KaonN_HistMode;
    std::map<std::string,std::map<int,TH1D*> > HadronE_HistMode;
    std::map<std::string,std::map<int,TH1D*> > PhotonCon_HistMode;

    std::map<std::string,std::map<int,TH1D*> > VisibleEnergyAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_AVCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_FVCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_FVBeforeAfter_HistMode;

    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_EnergyCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PhotonEnergyCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_ConversionGapCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_MuLenghtCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_NCCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_LeptonPlusPhotonCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PiZeroAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PhotonAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PhotonSmallAfter_HistMode;
    
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_NoShowerCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoPhotonCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoShowerCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_PhotonFVCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoLeptonsAfter_HistMode;

    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_NoShowerCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoPhotonCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoShowerCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_PhotonFVCut_HistMode;
     std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoLeptons_HistMode;



    std::map<std::string,std::map<int,TH1D*> > TrueNumber_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > TrueEnergyAll_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > TrueEnergy_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > CCQEEnergy_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_AVCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_FVCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_FVBefore_HistIntType;

    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_EnergyCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PhotonEnergyCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_ConversionGapCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_MuLenghtCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_NCCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_LeptonPlusPhotonCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > Weights_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PiZero_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_Photon_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PhotonSmall_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > LowNCEnergy_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > ProtonE_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > PionE_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > KaonE_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > ProtonN_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > PionN_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > KaonN_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > HadronE_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > PhotonCon_HistIntType;

    std::map<std::string,std::map<int,TH1D*> > VisibleEnergyAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_AVCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_FVCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_FVBeforeAfter_HistIntType;

    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_EnergyCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_PhotonEnergyCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_ConversionGapCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_MuLenghtCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_NCCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_LeptonPlusPhotonCutAfter_HistIntType;

    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_NoShowerCutAfter_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoPhotonCutAfter_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoShowerCutAfter_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_PhotonFVCutAfter_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoLeptonsAfter_HistIntType;

    std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_NoShowerCut_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoPhotonCut_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoShowerCut_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_PhotonFVCut_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergy_TwoLeptons_HistIntType;



    std::map<std::string,TH1D* >VisibleEnergyLep_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_AVCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_FVBefore_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_FVCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_EnergyCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_PhotonEnergyCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_ConversionGapCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_MuLenghtCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_NCCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_LeptonPlusPhotonCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_Selection_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_PiZero_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_Photon_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_PhotonSmall_Hist;

    std::map<std::string,TH1D* >VisibleEnergyLep_NoShowerCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_TwoPhotonCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_TwoShowerCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_PhotonFVCut_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_TwoLeptons_Hist;

    std::map<std::string,TH1D* >VisibleEnergyLepAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_AVCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_FVBeforeAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_FVCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_EnergyCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_PhotonEnergyCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_ConversionGapCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_MuLenghtCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_NCCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_LeptonPlusPhotonCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_SelectionAfter_Hist;

    std::map<std::string,TH1D* >VisibleEnergyLep_NoShowerCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_TwoPhotonCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_TwoShowerCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_PhotonFVCutAfter_Hist;
    std::map<std::string,TH1D* >VisibleEnergyLep_TwoLeptonsAfter_Hist;

    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_AVCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_FVCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_FVBefore_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_EnergyCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonEnergyCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_ConversionGapCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_MuLenghtCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_NCCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_LeptonPlusPhotonCut_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_Selection_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PiZero_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_Photon_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonSmall_HistMode;

    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLepAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_AVCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_FVCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_FVBeforeAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_EnergyCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonEnergyCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_ConversionGapCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_MuLenghtCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_NCCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_LeptonPlusPhotonCutAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_SelectionAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PiZeroAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonAfter_HistMode;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonSmallAfter_HistMode;
    

   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_NoShowerCutAfter_HistMode;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoPhotonCutAfter_HistMode;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoShowerCutAfter_HistMode;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonFVCutAfter_HistMode;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoLeptonsAfter_HistMode;

   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_NoShowerCut_HistMode;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoPhotonCut_HistMode;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoShowerCut_HistMode;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonFVCut_HistMode;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoLeptons_HistMode;



    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_AVCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_FVCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_FVBefore_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_EnergyCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonEnergyCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_ConversionGapCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_MuLenghtCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_NCCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_LeptonPlusPhotonCut_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_Selection_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PiZero_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_Photon_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonSmall_HistIntType;
    
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLepAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_AVCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_FVCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_FVBeforeAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_EnergyCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonEnergyCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_ConversionGapCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_MuLenghtCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_NCCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_LeptonPlusPhotonCutAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_SelectionAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PiZeroAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonSmallAfter_HistIntType;
	      
    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_NoShowerCutAfter_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoPhotonCutAfter_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoShowerCutAfter_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonFVCutAfter_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoLeptonsAfter_HistIntType;

    std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_NoShowerCut_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoPhotonCut_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoShowerCut_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_PhotonFVCut_HistIntType;
   std::map<std::string,std::map<int,TH1D*> >VisibleEnergyLep_TwoLeptons_HistIntType;

    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_SelectionAfter_HistIntType;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_SelectionAfter_HistMode;

    TH1D* VisibleEnergy_CosmicAVCut_Hist;
    TH1D* VisibleEnergy_CosmicFVCut_Hist;
    TH1D* VisibleEnergy_CosmicClyinderCut_Hist;
    TH1D* VisibleEnergy_CosmicdEdxCut_Hist;
    TH1D* VisibleEnergy_CosmicWeightCut_Hist;
    TH1D* VisibleEnergy_CosmicEnergyCut_Hist;
    TH1D* VisibleEnergy_SpillWindow_Hist;

    TH1D* VisibleEnergy_FNKPDecay_Hist;
    TH1D* VisibleEnergy_FNKMDecay_Hist;
    TH1D* VisibleEnergy_MuDecays_Hist;

    TH1D* NeutrinoT0;
    TH1D* CosmicShowerT0;

    TH2D* XZCosmic = new TH2D("XZCosmic","XZCosmic",100,-200,200,100,0,500);
    TH2D* YZCosmic = new TH2D("YZCosmic","YZCosmic",100,-200,200,100,0,500);
    TH2D* XYCosmic = new TH2D("XYCosmic","XYCosmic",100,-200,200,100,-200,200);
    TH1D* DCACosmic = new TH1D("DCACosmic","DCACosmic",100,0,200);

    TH2D* StartXZCosmic = new TH2D("StartXZCosmic","StartXZCosmic",100,-400,400,100,-100,600);
    TH2D* StartYZCosmic = new TH2D("StartYZCosmic","StartYZCosmic",100,-400,400,100,-100,600);
    TH2D* StartXYCosmic = new TH2D("StartXYCosmic","StartXYCosmic",100,-400,400,400,-400,200);
    TH2D* StartMCPXZCosmic = new TH2D("StartMCPXZCosmic","StartMCPXZCosmic",100,-400,400,100,-100,600);
    TH2D* StartMCPYZCosmic = new TH2D("StartMCPYZCosmic","StartMCPYZCosmic",100,-400,400,100,-100,600);
    TH2D* StartMCPXYCosmic = new TH2D("StartMCPXYCosmic","StartMCPXYCosmic",100,-400,400,400,-400,200);
    TH2D* EndMCPXZCosmic = new TH2D("EndMCPXZCosmic","EndMCPXZCosmic",100,-400,400,100,-100,600);
    TH2D* EndMCPYZCosmic = new TH2D("EndMCPYZCosmic","EndMCPYZCosmic",100,-400,400,100,-100,600);
    TH2D* EndMCPXYCosmic = new TH2D("EndMCPXYCosmic","EndMCPXYCosmic",100,-400,400,400,-400,200);

    


    TH2D* XZPhot = new TH2D("XZPhot","XZPhot",100,-200,200,100,0,500);
    TH2D* YZPhot = new TH2D("YZPhot","YZPhot",100,-200,200,100,0,500);
    TH2D* XYPhot = new TH2D("XYPhot","XYPhot",100,-200,200,100,-200,200);
    TH1D* DCAPhot = new TH1D("DCAPhot","DCAPhot",100,0,200);

    TH2D* XZNu = new TH2D("XZNu","XZNu",200,-400,400,200,-100,600);
    TH2D* YZNu = new TH2D("YZNu","YZNu",200,-400,400,200,-100,600);
    TH2D* XYNu = new TH2D("XYNu","XYNu",200,-400,400,200,-400,400);
    TH1D* DCANu = new TH1D("DCANu","DCANu",100,0,200);
  
   
    //Final selection histograms
    TH1D* VisibleEnergy_CosmicSelection_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_Selection_Hist;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_Selection_HistMode;
    std::map<std::string,std::map<int,TH1D*> > VisibleEnergy_Selection_HistIntType;
    std::vector<std::string> HistTypes = {"NuMu","InNuE","OscNuE","NCInNuE","NCOscNuE","NCNuMu","DirtNuMu","DirtInNuE","DirtOscNuE","DirtNCInNuE","DirtNCOscNuE","DirtNCNuMu","AllBackground","AllSignal","AllDirt"};
  };

  Config fConfig; //!< The config
  RootHistograms fRootHists;

  TRandom *randomnum_gen = new TRandom3();

  
  bool Select(const gallery::Event& ev, const simb::MCTruth& mctruth, unsigned truth_ind, std::map<int, const simb::MCParticle*>& mcparticles, std::map<int,double>& visible_mcparticles, NueInteraction& intInfo, std::vector<int>& leptontrackIDs);
  bool SelectCosmic(const sim::MCShower& mcs, std::map<int, const simb::MCParticle*>& mcparticles, NueSelection::NueInteraction& intInfo, std::vector<simb::MCTruth>& mcneutrinotruth);

  //Selection Functions 
  bool passFV(const TVector3 &v) {if(fConfig.ApplyFVCut){return containedInFV(v);} else return true;}
  bool passAV(const TVector3 &v) {if(fConfig.ApplyAVCut){return containedInAV(v);} else return true;}
  bool passNueIDEfficiencyCut(const simb::MCNeutrino& nu){if(fConfig.ApplyNueEfficiency){return NueIDEfficiencyCut(nu);} else return true;}
  bool passeEnergyCut(double& lepton_energy){if(fConfig.ApplyElectronEnergyCut){return eEnergyCut(lepton_energy);} else return true;}
  bool passPhotonEnergyCut(std::vector<int>& photons, std::map<int, const simb::MCParticle*>& mcparticles, int& photonTrackID,std::map<int,double>& visible_mcparticles){if(fConfig.ApplyPhotonEnergyCut){return PhotonEnergyCut(photons,mcparticles,photonTrackID, visible_mcparticles);} else return true;}
  bool passMuLengthCut(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth){if(fConfig.ApplyMuonLenghtCut){return MuLengthCut(mcparticles, mctruth);} else return true;}
  bool passdEdxCut(int pdgcode){if(fConfig.ApplydEdxCut){return dEdxCut(pdgcode);} else return true;}
  bool passConversionGapCut(std::map<int, const simb::MCParticle*>& mcparticles, int photonTrackID, const float& hadronic_energy, const simb::MCNeutrino& nu){if(fConfig.ApplyConversionGapCut){return ConversionGapCut(mcparticles,photonTrackID,hadronic_energy,nu);} else return true;}
  bool PassCosmicInFV(const sim::MCShower& mcs, TVector3& vertex, const simb::MCParticle*& mcparticle, double& time ){if(fConfig.ApplyCosmicFVCut){return CosmicInFV(mcs, vertex,mcparticle,time);} else return true;}
  bool PassCosmicCylinderCut(const sim::MCShower& mcs, TVector3& vertex, std::map<int, const simb::MCParticle*>& mcparticles){if(fConfig.ApplyCosmicCylinderCut){return CosmicCylinderCut(mcs, vertex, mcparticles);} else return true;}
  bool PassCosmicInSpillWindow(const sim::MCShower& mcs, std::map<int, const simb::MCParticle*>& mcparticles, NueSelection::NueInteraction& intInfo,std::vector<simb::MCTruth>& mcneutrinotruth, double& time){if(fConfig.ApplyCosmicInSpillWindowCut){return CosmicInSpillWindow(mcs,mcparticles,intInfo, mcneutrinotruth,time);} else return true;}

  //Cut Functions
  bool containedInFV(const TVector3 &v);
  bool containedInAV(const TVector3 &v);
  bool NueIDEfficiencyCut(const simb::MCNeutrino& nu);
  bool eEnergyCut(double leptonenergy);
  bool PhotonEnergyCut(std::vector<int>& photons, std::map<int, const simb::MCParticle*>& mcparticles, int &photonTrackID, std::map<int,double>& visible_mcparticles);
  bool ConversionGapCut(std::map<int, const simb::MCParticle*>& mcparticles, int photonTrackID, const float& hadronicE, const simb::MCNeutrino& nu);
  bool dEdxCut(int pdgcode);
  bool MuLengthCut(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth);
  bool CosmicInFV(const sim::MCShower& mcs, TVector3& vertex,  const simb::MCParticle*& mcparticle, double& time);
  bool CosmicInAV(const sim::MCShower& mcs, const simb::MCParticle*& mcparticle);

  bool CosmicCylinderCut(const sim::MCShower& mcs, TVector3& vertex, std::map<int, const simb::MCParticle*>& mcparticles);
  bool CosmicInSpillWindow(const sim::MCShower& mcs, std::map<int, const simb::MCParticle*>& mcparticles, NueSelection::NueInteraction& intInfo,std::vector<simb::MCTruth>& mcneutrinotruth, double& time);

  std::vector<int> findNeutralPions(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth);
  std::vector<int> findPhotons(std::vector<int>& pi_zeros,std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth, std::map<int,double>& visible_mcparticles);
  double PhotonVisibleEnergy(std::map<int, const simb::MCParticle*>& mcparticles, std::vector<int>& photons);
  double MECWeight(const simb::MCNeutrino& nu,std::string& Detector, NueSelection::NueInteraction& intInfo);
  double GENIEWeight(const simb::MCNeutrino& nu);

  void InitialiseHistograms();
  void FillHistograms(std::map<std::string,TH1D*>& HistMap, std::map<std::string,std::map<int,TH1D*> >& HistMapMode, std::map<std::string,std::map<int,TH1D*> >& HistMapIntType, const simb::MCNeutrino& nu, NueSelection::NueInteraction& intInfo, bool& booldirtevent);
  void FillHistograms(std::map<std::string,TH1D*>& HistMap, std::map<std::string,std::map<int,TH1D*> >& HistMapMode, std::map<std::string,std::map<int,TH1D*> >& HistMapIntType, const simb::MCNeutrino& nu, NueSelection::NueInteraction& intInfo,double Energy, bool &booldirtevent);

  void PrintInformation(const simb::MCTruth& mctruth, NueSelection::NueInteraction& intInfo);

  double DistanceToClosestSurface(const TVector3 &v);


  std::vector<int> ElectronNuDecays = {1,2,6,9};
  std::vector<int> MuonNuDecays = {3,4,5,7,8,10,11,12,13,14};

  //Measure CPU times 
  std::clock_t c_start;
  std::clock_t c_end;

  TRandom3 rand;
  double fPOT;
  TH1I*  fEventHist;
  bool dirtevent;
  int MCTruthCounter; 
  TFile* GenieDiffFile;

  std::map<int,TH1D*> GenieHists;
  std::map<std::string,std::map<std::string,std::map<int, std::vector<double> > > > MECWeight_map;
  std::map<int,std::vector<double> > GENIEWeight_map;
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelection__

