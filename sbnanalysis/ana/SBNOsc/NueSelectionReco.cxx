//C++ Includes 
#include <iostream>
#include <vector>

//Root Includes 
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/CrossValidation.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

//Framework Includes 
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "cetlib/search_path.h"

//Larsoft Includes 
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "RecoUtils/RecoUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

//SBN Includes 
#include "core/Event.hh"
#include "NueSelectionReco.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

    NueSelectionReco::NueSelectionReco() : SelectionBase(), EventCounter(0), NuCount(0) {}

    
    void NueSelectionReco::Initialize(fhicl::ParameterSet* config) {

      hello();
      fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("NueSelectionReco");

      // setup active volume bounding boxes
      std::vector<fhicl::ParameterSet> AVs =				\
        pconfig.get<std::vector<fhicl::ParameterSet> >("active_volumes");
      for (auto const& AV : AVs) {
        double xmin = AV.get<double>("xmin");
        double ymin = AV.get<double>("ymin");
        double zmin = AV.get<double>("zmin");
        double xmax = AV.get<double>("xmax");
        double ymax = AV.get<double>("ymax");
        double zmax = AV.get<double>("zmax");
        fConfig.active_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
      }

      std::vector<fhicl::ParameterSet> FVs =				\
        pconfig.get<std::vector<fhicl::ParameterSet> >("fiducial_volumes");
      for (auto const& FV : FVs) {
        double xmin = FV.get<double>("xmin");
        double ymin = FV.get<double>("ymin");
        double zmin = FV.get<double>("zmin");
        double xmax = FV.get<double>("xmax");
        double ymax = FV.get<double>("ymax");
        double zmax = FV.get<double>("zmax");
        fConfig.fiducial_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
      }

      //Set up the fcl parameters.
      fConfig.GlobalWeight            = pconfig.get<double>("GlobalWeight",1);
      fConfig.POTWeight               = pconfig.get<double>("POTWeight",1);
      fConfig.ElectronWeight          = pconfig.get<double>("ElectronWeight",1);
      fConfig.PhotonWeight            = pconfig.get<double>("PhotonWeight",1);
      fConfig.GlobalTimeOffset        = pconfig.get<double>("GlobalTimeOffset",1);
      fConfig.SpillTime               = pconfig.get<double>("SpillTime",1);
      fConfig.ReadoutWindowSize       = pconfig.get<double>("ReadoutWindowSize",1);
      fConfig.NSegments               = pconfig.get<double>("NSegments",1);
      fConfig.TrackEnergyCorrection   = pconfig.get<double>("TrackEnergyCorrection",1);
      fConfig.ShowerEnergyCorrection  = pconfig.get<double>("ShowerEnergyCorrection",1);
      fConfig.NutrinoEnergyCorrection = pconfig.get<double>("NutrinoEnergyCorrection",1);
      fConfig.EnergyConversion = pconfig.get<std::vector<float>>("EnergyConversion",{});
      fConfig.XOffset = pconfig.get<double>("XOffset",1);

      fConfig.GeneratorTag            = pconfig.get<std::string>("GeneratorTag","generator");
      fConfig.PandoraTag              = pconfig.get<std::string>("PandoraTag","pandora");
      fConfig.HitTag                  = pconfig.get<std::string>("HitTag","linecluster");
      fConfig.ShowerTag               = pconfig.get<std::string>("ShowerTag","tracs");
      fConfig.FluxTag                 = pconfig.get<std::string>("FluxTag","generator");
      fConfig.TrackTag                = pconfig.get<std::string>("TrackTag","pandoraTracktune");
      fConfig.TrackPIDTag             = pconfig.get<std::string>("TrackPIDTag","pandoraTrackPidtune");
      fConfig.CalorimetryTag          = pconfig.get<std::string>("CalorimetryTag","TrackRedoCaloUboonetune");

      fConfig.SetupMVA                = pconfig.get<bool>("SetupMVA",false);
      fConfig.Weightfile              = pconfig.get<std::string>("Weightfile","");
      fConfig.MVAMethod               = pconfig.get<std::string>("MVAMethod","");
      fConfig.MVAValuesNames          = pconfig.get<std::vector<std::string> >("MVAValuesNames",{});

      fConfig.ApplyFVCut                    = pconfig.get<bool>("ApplyFVCut",false);
      fConfig.ApplyOneShowerEnergyCut       = pconfig.get<bool>("ApplyOneShowerEnergyCut",false);
      fConfig.ApplyShowerResidualCut        = pconfig.get<bool>("ApplyShowerResidualCut",false);
      fConfig.ApplyConversionGapCut         = pconfig.get<bool>("ApplyConversionGapCut",false);
      fConfig.ApplydEdxCut                  = pconfig.get<bool>("ApplydEdxCut",false);
      fConfig.ApplyLengthCut                = pconfig.get<bool>("ApplyLengthCut",false);
      fConfig.ApplyOpeningAngleCut          = pconfig.get<bool>("ApplyOpeningAngleCut",false);
      fConfig.ApplyShowerEnergyCut          = pconfig.get<bool>("ApplyShowerEnergyCut",false);
      fConfig.ApplyShowerDensityGradientCut = pconfig.get<bool>("ApplyShowerDensityGradientCut",false);
      fConfig.ApplyShowerDensityPowerCut    = pconfig.get<bool>("ApplyShowerTrackLengthCut",false);
      fConfig.ApplyShowerTrackLengthCut     = pconfig.get<bool>("ApplyShowerTrackLengthCut",false);
      fConfig.ApplyShowerTrackWidthCut      = pconfig.get<bool>("ApplyShowerTrackWidthCut",false);
      fConfig.ApplyMaxTrackLengthCut        = pconfig.get<bool>("ApplyMaxTrackLengthCut",false);
      fConfig.ApplyMaxTrackPIDACut          = pconfig.get<bool>("ApplyMaxTrackPIDACut",false);
      fConfig.ApplyMaxTrackLengthPIDACut    = pconfig.get<bool>("ApplyMaxTrackLengthPIDACut",false);
      fConfig.ApplyMVACut                   = pconfig.get<bool>("ApplyMVACut",false);
      fConfig.ApplyNeutrinoPdgCodeCut       = pconfig.get<bool>("ApplyNeutrinoPdgCodeCut",false);
      fConfig.ApplyNumneutrinosCut          = pconfig.get<bool>("ApplyNumneutrinosCut",false);

      fConfig.IncludeCosmics                = pconfig.get<bool>("IncludeCosmics",false);
      fConfig.IncludeDirt                   = pconfig.get<bool>("IncludeDirt",false);
      fConfig.CosmicsOnly                   = pconfig.get<bool>("CosmicsOnly",false);
      fConfig.DirtOnly                      = pconfig.get<bool>("DirtOnly",false);
      fConfig.Verbose                       = pconfig.get<bool>("Verbose",false);
      fConfig.IntrinsicOnly                 = pconfig.get<bool>("IntrinsicOnly",false);
      fConfig.NuMuOnly                      = pconfig.get<bool>("NuMuOnly",false);
      fConfig.OscOnly                       = pconfig.get<bool>("OscOnly",false); 
      fConfig.NoDirt                        = pconfig.get<bool>("NoDirt",false);
      fConfig.FillHistograms                = pconfig.get<bool>("FillHistograms",false);
      fConfig.FillModeHistograms            = pconfig.get<bool>("FillModeHistograms",false);
      fConfig.FillIntTypeHistograms         = pconfig.get<bool>("FillIntTypeHistograms",false);
      fConfig.UseAllCosmics                 = pconfig.get<bool>("UseAllCosmics",false);

      fConfig.SecondaryShowerEnergyCut = pconfig.get<float>("SecondaryShowerEnergyCut",1);
      fConfig.ShowerResidualEnergyCut  = pconfig.get<float>("ShowerResidualEnergyCut",1);
      fConfig.ResiudalCut              = pconfig.get<float>("ResiudalCut",1);
      fConfig.NumShowersCut            = pconfig.get<float>("NumShowersCut",1);
      fConfig.dEdxCut                  = pconfig.get<float>("dEdxCut",1);
      fConfig.ConversionGapCut         = pconfig.get<float>("ConversionGapCut",1);
      fConfig.LengthCut                = pconfig.get<float>("LengthCut",1);
      fConfig.OpeningAngleCut          = pconfig.get<float>("OpeningAngleCut",1);
      fConfig.ShowerEnergyCut          = pconfig.get<float>("ShowerEnergyCut",1);
      fConfig.ShowerDensityGradientCut = pconfig.get<float>("ShowerDensityGradientCut",1);
      fConfig.ShowerDensityPowerCut    = pconfig.get<float>("ShowerDensityPowerCut",1);
      fConfig.ShowerTrackLengthCut     = pconfig.get<float>("ShowerTrackLengthCut",1);
      fConfig.ShowerTrackWidthCut      = pconfig.get<float>("ShowerTrackWidthCut",1);
      fConfig.MaxTrackLengthCut        = pconfig.get<float>("MaxTrackLengthCut",1);
      fConfig.MaxTrackPIDACut          = pconfig.get<float>("MaxTrackPIDACut",1);
      fConfig.MaxTrackLengthPIDACut    = pconfig.get<float>("MaxTrackLengthPIDACut",1);
      fConfig.MVACut                   = pconfig.get<float>("MVACut",1);
      fConfig.NeutrinoPdgCodeCut       = pconfig.get<float>("NeutrinoPdgCodeCut",1);
      fConfig.NumneutrinosCut          = pconfig.get<float>("NumneutrinosCut",1);

      //Set up the selection histograms 
      fOutputFile->cd();

      //Initilise the histograms 
      InitialiseHistograms();


      //Time the CPU.
      //c_start = std::clock();

      MCTruthCounter =0;

      //Setup the services
      fDetProp    = fProviderManager->GetDetectorPropertiesProvider();
      fpi_service = fProviderManager->GetParticleInventoryProvider();

      //Setup the mva 
      if(fConfig.SetupMVA){
	for(auto const& MVAValuesName: fConfig.MVAValuesNames){
	  std::cout << "adding variable: " << MVAValuesName << std::endl;
	  MVAValues[MVAValuesName] = -999;
	}
	SetupMVA();
      }
      std::cout << "gestting to the end of the set up" << std::endl;

    }
    
    
    void NueSelectionReco::Finalize(){

      std::cout << "Neutrinos selected: " <<  NuCount << " in " << EventCounter << " readout events" << std::endl;
      std::cout << "MCTruth Events count was: " <<  MCTruthCounter << std::endl;

      //Set up the selection histograms                                                  
      fOutputFile->cd();    
      gDirectory->mkdir("Histograms");
      gDirectory->mkdir("ModeHistograms");
      gDirectory->mkdir("InteractionTypeHistograms");
      fOutputFile->cd("Histograms");

      if(fConfig.FillHistograms){
	
	fRootHists.XDiff->Write();
	fRootHists.YDiff->Write();
	fRootHists.ZDiff->Write();

	for(int i=0; i<fRootHists.HistTypes.size(); ++i){

	  fOutputFile->cd("Histograms");
	  fRootHists.VisibleEnergy_MCTruth_Hist[fRootHists.HistTypes[i]]->Write(); 

	  fRootHists.TrueNumber_Hist[fRootHists.HistTypes[i]]->Write();
	  fRootHists.VisibleEnergy_BeforeSel_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_AfterSel_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_AfterSelOne_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_AfterSelExtra_Hist[fRootHists.HistTypes[i]]->Write(); 

	  fRootHists.VisibleEnergy_FVRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_FVPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_GtrOneShowerRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_GtrOneShowerPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_OneShowerECutRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_OneShowerECutPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_OneShowerResidualRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_OneShowerResidualPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ConversionGapRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ConversionGapPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_dEdxRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_dEdxPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_LengthRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_LengthPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_OpeningAngleRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_OpeningAnglePassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerEnergyRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerEnergyPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerDensityGradientPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerDensityPowerPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerTrackLengthPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_ShowerTrackWidthPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_MaxTrackLengthRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_MaxTrackLengthPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_MaxTrackPIDARemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_MaxTrackPIDAPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_MVARemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_MVAPassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_NumneutrinosRemoved_Hist[fRootHists.HistTypes[i]]->Write(); 
	  fRootHists.VisibleEnergy_NumneutrinosPassed_Hist[fRootHists.HistTypes[i]]->Write(); 

	  if(fConfig.FillModeHistograms){
	    for(int mode=-1; mode<6; ++mode){
	      fOutputFile->cd("ModeHistograms");
	      fRootHists.VisibleEnergy_MCTruth_HistMode[fRootHists.HistTypes[i]][mode]->Write();

	      fRootHists.TrueNumber_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_BeforeSel_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_AfterSel_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_AfterSelOne_HistMode[fRootHists.HistTypes[i]][mode]->Write();

	      fRootHists.VisibleEnergy_AfterSelExtra_HistMode[fRootHists.HistTypes[i]][mode]->Write();

	      fRootHists.VisibleEnergy_FVRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_FVPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_GtrOneShowerRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_GtrOneShowerPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_OneShowerECutRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_OneShowerECutPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_OneShowerResidualRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_OneShowerResidualPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ConversionGapRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ConversionGapPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_dEdxRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_dEdxPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_LengthRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_LengthPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_OpeningAngleRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_OpeningAnglePassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerEnergyRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerEnergyPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerDensityGradientPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerDensityPowerPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerTrackLengthPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_ShowerTrackWidthPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MaxTrackLengthRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MaxTrackLengthPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MaxTrackPIDARemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MaxTrackPIDAPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MVARemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_MVAPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_NumneutrinosRemoved_HistMode[fRootHists.HistTypes[i]][mode]->Write();
	      fRootHists.VisibleEnergy_NumneutrinosPassed_HistMode[fRootHists.HistTypes[i]][mode]->Write();

	    }
	  }

	  if(fConfig.FillIntTypeHistograms){
	    for(int inttype=1000; inttype<1101; ++inttype){
	      fOutputFile->cd("InteractionTypeHistograms");

	      if(fRootHists.TrueNumber_HistIntType[fRootHists.HistTypes[i]][inttype]->GetEntries() == 0){continue;} 
	      fRootHists.VisibleEnergy_MCTruth_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();


	      fRootHists.TrueNumber_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_BeforeSel_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_AfterSel_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_AfterSelOne_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
 	      fRootHists.VisibleEnergy_AfterSelExtra_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();

	      fRootHists.VisibleEnergy_FVRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_FVPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_GtrOneShowerRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_GtrOneShowerPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_OneShowerECutRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_OneShowerECutPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_OneShowerResidualRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_OneShowerResidualPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ConversionGapRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ConversionGapPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_dEdxRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_dEdxPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_LengthRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_LengthPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_OpeningAngleRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_OpeningAnglePassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerEnergyRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerEnergyPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerDensityGradientPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerDensityPowerPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerTrackLengthPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_ShowerTrackWidthPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MaxTrackLengthRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MaxTrackLengthPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MaxTrackPIDARemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MaxTrackPIDAPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MVARemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_MVAPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();

	      fRootHists.VisibleEnergy_NumneutrinosRemoved_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();
	      fRootHists.VisibleEnergy_NumneutrinosPassed_HistIntType[fRootHists.HistTypes[i]][inttype]->Write();


	    }
	  }
	}
	
	fOutputFile->cd("Histograms");

      }
    }


    void NueSelectionReco::SetupMVA(){
      
      reader = new TMVA::Reader( "!Color:!Silent" );

      //Define the variables.
      for(auto& MVAValue: MVAValues){
	std::cout << "Adding Variable: " << MVAValue.first << std::endl;
	reader->AddVariable(MVAValue.first,&MVAValue.second);
      }

      std::string fname;
      cet::search_path sp("FW_SEARCH_PATH");
      if (!sp.find_file((std::string) fConfig.Weightfile, fname)) {
	throw cet::exception("NueSelectionReco") << "Could not find the weight file: " << fname << std::endl; ;
      }
      reader->BookMVA(fConfig.MVAMethod, fname);
      return; 
    }


    bool NueSelectionReco::ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco){

      fDetProp    = fProviderManager->GetDetectorPropertiesProvider();
      fpi_service = fProviderManager->GetParticleInventoryProvider();


      if(fConfig.Verbose){
	if (EventCounter % 10 == 0) {
	  std::cout << "NueSelectionReco: Processing event " << EventCounter << " "
		    << "(" << NuCount << " neutrinos selected)"
		    << std::endl;
	}

	std::cout << "New Event" << std::endl;
      }

      bool selected = false;
      int truth_int = 0; 
      EventCounter++;

      //Collect all the truths in a map where the key is the origin and energy and is th value is the interger used in the interaction
      std::map<simb::Origin_t,std::map<int,int> > mctruthkeys;
      int mc_iter=0;
      for(int fTruthTag=0; fTruthTag<fTruthTags.size(); ++fTruthTag){
	auto const& mctruthsHandle =  *ev.getValidHandle< std::vector<simb::MCTruth> >(fTruthTags[fTruthTag]);
	for(auto const& mctruth: mctruthsHandle){
	  if(fConfig.Verbose){
	    std::cout << "mctruth: " << std::endl;
	    std::cout << mctruth << std::endl;
	  }

	  double Energy = 0;
	  if(mctruth.NeutrinoSet()){
	    Energy = mctruth.GetNeutrino().Nu().E();
	  }
	  mctruthkeys[mctruth.Origin()][Energy] = ++mc_iter;
	}
      }
      
      //Association truth and flux
      auto const& mctruthsHandle =  ev.getValidHandle< std::vector<simb::MCTruth> >(fConfig.GeneratorTag);

      art::FindManyP<simb::MCFlux> fof(mctruthsHandle, ev, fConfig.GeneratorTag);
      if(!fof.isValid()){
	throw cet::exception("NueSelectionReco") << "Flux MCTruth association is somehow not valid. Stopping";
	return false;
      }

      //Save a histogram of the values before reconstruction.
      int iter = -1;
      std::vector<art::Ptr<simb::MCTruth> > mctruth_ptrs;
      art::fill_ptr_vector(mctruth_ptrs,mctruthsHandle);

      for(auto const& mctruth: mctruth_ptrs){   
	iter++;
	bool passAV = containedInAV(mctruth->GetNeutrino().Nu().Position().Vect());
	if(passAV){
	  std::vector<art::Ptr<simb::MCFlux> > mcflux = fof.at(mctruth.key());
	  int initnu = mcflux.at(0)->fntype;
	  int fnd = mcflux.at(0)->fndecay;
	  int pdg = mctruth->GetNeutrino().Nu().PdgCode();
	  float True_energy = mctruth->GetNeutrino().Nu().E();
	  NueSelectionReco::NueInteraction intInfo({0, 0, 0, 0, 1,initnu,fnd, True_energy,mctruth,false}); 
	  FillHistograms(fRootHists.VisibleEnergy_MCTruth_Hist,fRootHists.VisibleEnergy_MCTruth_HistMode,fRootHists.VisibleEnergy_MCTruth_HistIntType,intInfo);
	}
      }
      
	
      //Identify the neutrino particles 
      auto const &pfpHandle =						\
      	ev.getValidHandle<std::vector<recob::PFParticle> >(fConfig.PandoraTag);
      std::vector<art::Ptr<recob::PFParticle> > pfps;
      art::fill_ptr_vector(pfps,pfpHandle);

      std::vector<art::Ptr<recob::PFParticle> > neutrinos;
      for(auto const& pfp: pfps){
	if(TMath::Abs(pfp->PdgCode()) == 12 || TMath::Abs(pfp->PdgCode()) == 14){
	  neutrinos.push_back(pfp);
	}
      }
      
      //Identify the cluster data products 
      auto const &clusterHandle =  \
	ev.getValidHandle<std::vector<recob::Cluster> >(fConfig.PandoraTag);

      //Create the assoications
      art::FindManyP<recob::Cluster> fmpfc(pfpHandle, ev, fConfig.PandoraTag);
      art::FindManyP<recob::Hit> fmch(clusterHandle, ev, fConfig.PandoraTag);

      //Make a pfp map
      std::map<int, art::Ptr<recob::PFParticle> > pfp_map;
      for(auto const& pfp: pfps){
	pfp_map[pfp->Self()] = pfp;
      }

      //Keep track if a neutrino is already selected.
      std::vector<int> SelectedKeys;

      //Match the neutrinos to the truth 
      for(auto const& neutrino: neutrinos){

	if(fConfig.Verbose){
	  std::cout << "Particle list: " << std::endl;
	  const sim::ParticleList& particles = fpi_service->ParticleList();
	  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
	    const simb::MCParticle *particle = particleIt->second;
	    std::cout << "Track ID: " << particle->TrackId() << " pdg: " << particle->PdgCode() << " E: " <<  particle->E() << " mother: " << particle->Mother() << std::endl;
	  }
	}


	//Get the hits
	std::vector<art::Ptr<recob::Hit> > pfpNeutrinoHits;
	GetSliceHits(neutrino, pfpNeutrinoHits, pfp_map, fmpfc, fmch);

	if(fConfig.Verbose){
	  std::cout << " number of hits: " << pfpNeutrinoHits.size() << std::endl;
	}
	
	//Get the mc track id which has depoisted the most energy
	int MostHitsTrackID = TMath::Abs(SBNRecoUtils::TrueParticleIDFromTotalRecoHits(*fProviderManager,pfpNeutrinoHits, false));

	//Backtracker no set up I presume lets try to match with geomtry
	if(MostHitsTrackID == 99999){
	  if(fConfig.Verbose){
	    std::cout << "failed ot get info from backtracker trying geo" << std::endl;
	  }

	  //Assocations between pfparticle and vertex.
	  art::FindManyP<recob::Vertex> fmv(pfpHandle, ev, fConfig.PandoraTag);
	  if(!fmv.isValid()){
	    throw cet::exception("NueSelectionReco") << "Vertex and PF particle association is somehow not valid. Stopping";
	    return false;
	  }
	  
	  //Get the neutrino vertex 
	  std::vector<art::Ptr<recob::Vertex> > vertex = fmv.at(neutrino.key());
	  if(vertex.size() > 1){
	    throw cet::exception("NueSelectionReco") << "we have too many recob vertex for pfparticles";
	    return false;
	  }
      
	  double vtx_xyz[3] = {-999,-999,-999};
	  vertex[0]->XYZ(vtx_xyz);
	  TVector3 vertex_position = {vtx_xyz[0],vtx_xyz[1],vtx_xyz[2]};

	  if(fConfig.Verbose){ 
	    std::cout << "vertex position X: " << vtx_xyz[0] << " Y: " << vtx_xyz[1] << " Z: " << vtx_xyz[2] << std::endl;
	  } 

	  MostHitsTrackID = TMath::Abs(SBNRecoUtils::FindClosestParticle(*fProviderManager,vertex_position));

	  if(MostHitsTrackID == 99999){
	    throw cet::exception("NueSelectionReco") << "I've tried to match the reco neutrino using calo and geomtry and failed. I guess something is wrong";
	    return false;
	  }
	}

	if(fConfig.Verbose){ 
	  std::cout << "MostHitsTrackID: " << MostHitsTrackID << std::endl;
	}

	//Get the MCTruth associated to this particle
	const art::Ptr< simb::MCTruth > mctruth = fpi_service->TrackIdToMCTruth_P(MostHitsTrackID);

	if(fConfig.Verbose){
	  std::cout << "Matched MC Truth is" << std::endl;
	  std::cout << mctruth << std::endl;
	}



	//Get the truth iter for the interaction
	int truth_iter = 999;
	if(mctruth->NeutrinoSet()){
	  truth_iter = mctruthkeys[mctruth->Origin()][mctruth->GetNeutrino().Nu().E()];
	}
	else{
	  truth_iter = mctruthkeys[mctruth->Origin()][0];
	}

	//Calculate Energies
	double Hadronic_energy = 0;
	double Shower_energy   = 0;
	double Leptonic_energy = 0;
	double Other_energy    = 0;
	double True_energy     = -999;

	double  weight = 1;
	int initnu = -999;
	int fnd    = -999;
	int pdg    = -999;
	

	//Get the initial interaction infromation
	if(mctruth->NeutrinoSet()){
	  std::vector<art::Ptr<simb::MCFlux> > mcflux = fof.at(mctruth.key());
	  initnu = mcflux.at(0)->fntype;
	  fnd = mcflux.at(0)->fndecay;
	  pdg = mctruth->GetNeutrino().Nu().PdgCode();
	  True_energy = mctruth->GetNeutrino().Nu().E();
	}

	//Identify if the event was a dirt event or not.
	bool booldirtevent = false;
	bool passAV = containedInAV(mctruth->GetNeutrino().Nu().Position().Vect());
	if(!passAV){
	  if(fConfig.Verbose){std::cout << "is dirt" << std::endl;}
	  booldirtevent = true;
	  if(fConfig.NoDirt){
	    if(fConfig.Verbose){std::cout << "removing due to dirt" << std::endl;}
	    return false;
	  }
	}
       

	//Fill the interaction information (a mix between truth and reco)
	NueSelectionReco::NueInteraction intInfo({Hadronic_energy, Shower_energy, Leptonic_energy, Other_energy, weight,initnu,fnd, True_energy,mctruth,booldirtevent}); 
	

	///Weight with the POT
	intInfo.weight *= fConfig.POTWeight;
	
	//Add the Global weightings if any.
	intInfo.weight *= fConfig.GlobalWeight;
	
	//Check the graphs before 
	FillHistograms(fRootHists.VisibleEnergy_BeforeSel_Hist,fRootHists.VisibleEnergy_BeforeSel_HistMode,fRootHists.VisibleEnergy_BeforeSel_HistIntType,intInfo);

	//Perform the selection 
	bool selection = Select(ev, neutrino, pfp_map,intInfo);

	//Calculate the energy in the neutrino now.
	CalculateEnergy(ev, neutrino, pfp_map,intInfo);

	if(selection){

	  if(fConfig.Verbose){
	    std::cout << "Even selected Energy is" << intInfo.GetRecoEnergy()*fConfig.NutrinoEnergyCorrection << " weight: " << intInfo.weight;
	    std::cout << "Matched MC Info: " << mctruth << std::endl;
	  }
	  
	  //Calculate the energy in the neutrino now.
	  //  CalculateEnergy(ev, neutrino, pfp_map,intInfo);

	  if(std::find(SelectedKeys.begin(),SelectedKeys.end(),intInfo.mctruth.key()) == SelectedKeys.end()){
	    
	    SelectedKeys.push_back(intInfo.mctruth.key());
	    FillHistograms(fRootHists.VisibleEnergy_AfterSelOne_Hist,fRootHists.VisibleEnergy_AfterSelOne_HistMode,fRootHists.VisibleEnergy_AfterSelOne_HistIntType,intInfo);
	  }
	  else{
	    FillHistograms(fRootHists.VisibleEnergy_AfterSelExtra_Hist,fRootHists.VisibleEnergy_AfterSelExtra_HistMode,fRootHists.VisibleEnergy_AfterSelExtra_HistIntType,intInfo);
	  }

	  //check the graphs after
	  FillHistograms(fRootHists.VisibleEnergy_AfterSel_Hist,fRootHists.VisibleEnergy_AfterSel_HistMode,fRootHists.VisibleEnergy_AfterSel_HistIntType,intInfo);

	  
	  //Make the reco interaction event. 
	  event::RecoInteraction reco_interaction(truth_int);
	  reco_interaction.weight = intInfo.weight; 
	  reco_interaction.reco_energy= (intInfo.GetRecoEnergy()*fConfig.NutrinoEnergyCorrection)/1000.;
	  reco.push_back(reco_interaction);
	  selected =  selection;
	  
	  NuCount++;
	}
      }
    

      return selected;
    }

    // void NueSelectionReco::CalculateEnergy(const art::Ptr<recob::PFParticle>& pfp, 
    // 					std::vector< art::Ptr< recob::Hit> >& pfpHits,
    // 					std::map<int, art::Ptr<recob::PFParticle> >& pfpMap,
    // 					art::FindManyP<recob::Cluster>& fmpfc,
    // 					art::FindManyP<recob::Hit>& fmch){
      
    //   // Get the hits from the PFParticle
    //   const std::vector< art::Ptr< recob::Cluster> >& clusters = fmpfc.at(pfp.key());
    //   for (const auto& cluster: clusters){
    // 	const std::vector< art::Ptr< recob::Hit> >& hits = fmch.at(cluster.key());
    // 	pfpHits.insert(pfpHits.end(), hits.begin(), hits.end());
    //   }
      
    //   // Get the daughters
    //   const std::vector<long unsigned int> daughters = pfp->Daughters();
    //   for (const auto daughterIter: daughters){
    // 	art::Ptr<recob::PFParticle> daughter =  pfpMap.at(daughterIter);
    // 	// Get the hits from the daughters
    // 	GetSliceHits(daughter, pfpHits, pfpMap, fmpfc,fmch);
    //   }
    //   return;
    // }

    void NueSelectionReco::CalculateEnergy(const gallery::Event& ev,
					   const art::Ptr<recob::PFParticle>& neutrino,
					   std::map<int, art::Ptr<recob::PFParticle> >& pfp_map,
					   NueSelectionReco::NueInteraction& intInfo){
      
      //Ceebs to pass the handle 
      auto const &pfpHandle =						\
      	ev.getValidHandle<std::vector<recob::PFParticle> >(fConfig.PandoraTag);

      auto const &calo =						\
      	ev.getValidHandle<std::vector<anab::Calorimetry> >(fConfig.CalorimetryTag);
      std::vector<art::Ptr<anab::Calorimetry> > calo_vec;
      art::fill_ptr_vector(calo_vec,calo);

      //Association between Showers and pfParticle
      art::FindManyP<recob::Shower> fmsh(pfpHandle, ev, fConfig.ShowerTag);
      if(!fmsh.isValid()){
	throw cet::exception("NueSelectionReco") << "Shower and PF particle association is somehow not valid. Stopping";
	return;
      }
      auto const &trackListHandle =					\
      	ev.getValidHandle<std::vector<recob::Track> >(fConfig.TrackTag);
      art::FindManyP<recob::Track> fmt(pfpHandle, ev, fConfig.TrackTag);
      if(!fmt.isValid()){
	throw cet::exception("NueSelectionReco") << "Track and PF particle association is somehow not valid. Stopping";
	return;
      }
      art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, ev, fConfig.CalorimetryTag);
      if(!fmcal.isValid()){
	throw cet::exception("NueSelectionReco") << "Track and Calo association is somehow not valid. Stopping";
	return;
      }
     
      float HadronicEnergy = 0;
      float ShowerEnergy   = 0; 

      //Loop over the particles
      for(auto const& pfp: pfp_map){
	art::Ptr<recob::PFParticle> parent   = pfp.second;
	art::Ptr<recob::PFParticle> temp_pfp = pfp.second;
      
	//Get the primary particle
	while(!parent->IsPrimary()){
	  temp_pfp = parent;
	  if(pfp_map.find(parent->Self()) == pfp_map.end()){parent = temp_pfp; break;}
	  parent = pfp_map[parent->Parent()];
	}

	//Check we come from the neutrino
	if(parent->Self() != neutrino->Self()){continue;}


	//Calculate the track energy
	if(pfp.second->PdgCode() == 13){

	  //Get the track info 
	  std::vector<art::Ptr<recob::Track> > track = fmt.at(pfp.second.key());
	  if(track.size() == 0){continue;}
	  if(track.size() > 1){
	    throw cet::exception("NueSelectionReco") << "we have too many recob tracks for pfparticles";
	    return;
	  }
	  //Get the calorimetry info
	  std::vector<art::Ptr<anab::Calorimetry> > calo = fmcal.at(track[0].key());

	  float TrackEnergy = -999;
	  if(calo.size() > 0){
	    TrackEnergy = calo.at(calo.size()-1)->KineticEnergy();
	  }
	  if(TrackEnergy < 0){
	    float max_tack_energy = 0;
	    int trackit =0;
	    for(auto const& cal: calo){
	      if(cal->KineticEnergy() < 0){continue;}
	      max_tack_energy += cal->KineticEnergy();
	      ++trackit;
	    }
	    if(trackit == 0){
	      max_tack_energy = -999;
	    }
	    else{
	      max_tack_energy /= trackit;
	      TrackEnergy = max_tack_energy;

	    }
	  }

	  if(TrackEnergy > 0){
	    HadronicEnergy += TrackEnergy*fConfig.TrackEnergyCorrection;
	  }
	}


	//Calculate the shower energy 
	if(pfp.second->PdgCode() == 11){
	  
	  //Do we have a corresponding shower particle.
	  std::vector<art::Ptr<recob::Shower> > shower = fmsh.at(pfp.second.key());
	  
	  //Did we succeed at characterising the shower particle?
	  if(shower.size() == 0){std::cout << "no reco shower" << std::endl; continue;}
	  
	  //If we have two then our charactisation did a silly.
	  if(shower.size() != 1){
	    throw cet::exception("NueSelectionReco") << "we have too many recob showers for pfparticles";
	    return;
	  }

	  //Apply reconstruction energy cut.                 
	  const int ShowerBest_Plane = shower[0]->best_plane();
	  const std::vector<double> ShowerEnergyPlanes = shower[0]->Energy(); 
	  float max_shower_energy = ShowerEnergyPlanes[0];
	  if(ShowerBest_Plane != -999){
	    if((int) ShowerEnergyPlanes.size() >= ShowerBest_Plane+1){
	      max_shower_energy = ShowerEnergyPlanes[ShowerBest_Plane];
	    }
	  }
	  if(max_shower_energy < 0){
	    max_shower_energy = 0;
	    int showerit =0;
	    for(auto const& showerenergy: ShowerEnergyPlanes){
	      if(showerenergy < 0){continue;}
	      max_shower_energy += showerenergy;
	      ++showerit;
	    }
	    if(showerit == 0){
	      max_shower_energy = -999;
	    }
	    else{
	      max_shower_energy /= showerit;
	    }
	  }
	  ShowerEnergy += max_shower_energy*fConfig.ShowerEnergyCorrection;
	}
      }

      if(fConfig.Verbose){
	std::cout << "HadronicEnergy: " << HadronicEnergy << std::endl;
	std::cout << "ShowerEnergy:   " << ShowerEnergy   << std::endl;
      }

      intInfo.hadronic_energy = HadronicEnergy;
      intInfo.shower_energy   = ShowerEnergy;
      return;
    }

      

    
    bool NueSelectionReco::Select(const gallery::Event& ev,
				  const art::Ptr<recob::PFParticle>& neutrino,
				  std::map<int, art::Ptr<recob::PFParticle> >& pfp_map,
				  NueSelectionReco::NueInteraction& intInfo){

      //Ceebs to pass the handle 
      auto const &pfpHandle =						\
      	ev.getValidHandle<std::vector<recob::PFParticle> >(fConfig.PandoraTag);


      //Assocations between pfparticle and vertex.
      art::FindManyP<recob::Vertex> fmv(pfpHandle, ev, fConfig.PandoraTag);
      if(!fmv.isValid()){
	throw cet::exception("NueSelectionReco") << "Vertex and PF particle association is somehow not valid. Stopping";
	return false;
      }

      //Get the neutrino vertex 
      std::vector<art::Ptr<recob::Vertex> > vertex = fmv.at(neutrino.key());
      if(vertex.size() > 1){
	throw cet::exception("NueSelectionReco") << "we have too many recob vertex for pfparticles";
	return false;
      }
      
      double vtx_xyz[3] = {-999,-999,-999};
      vertex[0]->XYZ(vtx_xyz);
      TVector3 vertex_position = {vtx_xyz[0],vtx_xyz[1],vtx_xyz[2]};
      
      //#################
      //##### FV Cut ####
      //#################

      //Check If the event was in the fiducal volume.
      TVector3 vertex_position_corrected = {vtx_xyz[0]-fConfig.XOffset,vtx_xyz[1],vtx_xyz[2]};

      fRootHists.XDiff->Fill(vertex_position_corrected.X()-intInfo.mctruth->GetNeutrino().Nu().Position().Vect().X());
      fRootHists.YDiff->Fill(vertex_position_corrected.Y()-intInfo.mctruth->GetNeutrino().Nu().Position().Vect().Y());
      fRootHists.ZDiff->Fill(vertex_position_corrected.Z()-intInfo.mctruth->GetNeutrino().Nu().Position().Vect().Z());


      bool pass_FV = passFV(vertex_position_corrected);
      if(!pass_FV){
	if(fConfig.Verbose){
	  std::cout << "Failed the FV cut. Event not selected with true position X:" << intInfo.mctruth->GetNeutrino().Nu().Position().Vect().X() << " Y: " << intInfo.mctruth->GetNeutrino().Nu().Position().Vect().Y() << " Z: " << intInfo.mctruth->GetNeutrino().Nu().Position().Vect().Z() << std::endl;
	  std::cout << "Reco Position X: " << vertex_position.X() << " Y: " << vertex_position.Y() << " Z: " << vertex_position.Z() << std::endl;
	}
	FillHistograms(fRootHists.VisibleEnergy_FVRemoved_Hist,fRootHists.VisibleEnergy_FVRemoved_HistMode,fRootHists.VisibleEnergy_FVRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_FVPassed_Hist,fRootHists.VisibleEnergy_FVPassed_HistMode,fRootHists.VisibleEnergy_FVPassed_HistIntType,intInfo);

      //Association between Showers and pfParticle
      art::FindManyP<recob::Shower> fmsh(pfpHandle, ev, fConfig.ShowerTag);
      if(!fmsh.isValid()){
	throw cet::exception("NueSelectionReco") << "Shower and PF particle association is somehow not valid. Stopping";
	return false ;
      }


      //Identify any primary showers 
      const std::vector<long unsigned int> daughters = neutrino->Daughters();
      
      std::vector<art::Ptr<recob::Shower> > neutrino_showers; 
      for(auto const& daughter: daughters){
	
	//is the daughter a shower.
	if(pfp_map[daughter]->PdgCode() != 11){continue;}

	//Do we have a corresponding shower particle.
	std::vector<art::Ptr<recob::Shower> > shower = fmsh.at(pfp_map[daughter].key());
	
	//Did we succeed at characterising the shower particle?
	if(shower.size() == 0){
	  if(fConfig.Verbose){
	    std::cout << " this event has a pfp electron but no reco shower" << std::endl; 
	  }
	    continue;
	}
	
	//If we have two then our charactisation did a silly.
	if(shower.size() != 1){
	  throw cet::exception("NueSelectionReco") << "we have too many recob showers for pfparticles";
	  return false;
	}

	//Remove shower with energy less than 10 MeV. what is a a shower when it is just a few hits.
	if(shower[0]->Energy().at(shower[0]->best_plane()) < 10){continue;}

	//Then we have a shower
	neutrino_showers.push_back(shower[0]);
      }

      if(fConfig.Verbose){
	std::cout << "number of reco showers in event: " <<  neutrino_showers.size() << std::endl;
      }
      
      //########################################
      //########## Global TPC Cuts #############
      //########################################

      //Remove any event without a shower
      if(neutrino_showers.size() == 0){
	if(fConfig.Verbose){std::cout << "No showers in the event side" << std::endl;}
	FillHistograms(fRootHists.VisibleEnergy_GtrOneShowerRemoved_Hist,fRootHists.VisibleEnergy_GtrOneShowerRemoved_HistMode,fRootHists.VisibleEnergy_GtrOneShowerRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_GtrOneShowerPassed_Hist,fRootHists.VisibleEnergy_GtrOneShowerPassed_HistMode,fRootHists.VisibleEnergy_GtrOneShowerPassed_HistIntType,intInfo);


      //Order the showers with regards to their energy
      std::sort(neutrino_showers.begin(), neutrino_showers.end(),[](const art::Ptr<recob::Shower> & a, const art::Ptr<recob::Shower> & b){ 
	  return a->Energy().at(a->best_plane()) > b->Energy().at(b->best_plane()); 
	});
      

      //Remove event if there is more than x showers above a certain energy threshold
      bool pass_osec = passOneShowerEnergyCut(neutrino_showers);
      if(!pass_osec){
	if(fConfig.Verbose){
	  std::cout << "too many showers one shower cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_OneShowerECutRemoved_Hist,fRootHists.VisibleEnergy_OneShowerECutRemoved_HistMode,fRootHists.VisibleEnergy_OneShowerECutRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_OneShowerECutPassed_Hist,fRootHists.VisibleEnergy_OneShowerECutPassed_HistMode,fRootHists.VisibleEnergy_OneShowerECutPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "passed the one shower cut" << std::endl;
      }


      float neutrinopdg = neutrino->PdgCode();
      //Remove event if there is more than x showers above a certain energy threshold
      bool pass_neutrinopdg = passNeutrinoPdgCodeCut(neutrinopdg);
      if(!pass_neutrinopdg){
	if(fConfig.Verbose){
	  std::cout << "failed neutrinopdgcode cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_Hist,fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_HistMode,fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_Hist,fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_HistMode,fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "passed the neutrinopdgcode cut" << std::endl;
      }

      //Count number of neutrinos
      float numneutrinos = 0;
      for(auto const& pfp: pfp_map){ 
	if(TMath::Abs(pfp.second->PdgCode()) == 12 || TMath::Abs(pfp.second->PdgCode()) == 14){
	  ++numneutrinos;
	}
      }
      bool pass_numneutrinos = passNumneutrinosCut(numneutrinos);
      if(!pass_numneutrinos){
	if(fConfig.Verbose){
	  std::cout << "failed numneutrinoscode cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_NumneutrinosRemoved_Hist,fRootHists.VisibleEnergy_NumneutrinosRemoved_HistMode,fRootHists.VisibleEnergy_NumneutrinosRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_NumneutrinosPassed_Hist,fRootHists.VisibleEnergy_NumneutrinosPassed_HistMode,fRootHists.VisibleEnergy_NumneutrinosPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "passed the numneutrinoscode cut" << std::endl;
      }



      //Remove events via the one shower residual cut
      //Removal event if there is more than x shower above a certain energy threhold 
      //that are within y cm of the a best fit line from the cloest shower and vertex. 
      double numshowers = -1;
      bool pass_res = passShowerResiudalCut(neutrino_showers,numshowers,vertex_position);
      if(!pass_res && fConfig.ApplyShowerResidualCut){
	if(fConfig.Verbose){
	  std::cout << "too many showers one shower res cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_OneShowerResidualRemoved_Hist,fRootHists.VisibleEnergy_OneShowerResidualRemoved_HistMode,fRootHists.VisibleEnergy_OneShowerResidualRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_OneShowerResidualPassed_Hist,fRootHists.VisibleEnergy_OneShowerResidualPassed_HistMode,fRootHists.VisibleEnergy_OneShowerResidualPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "one shower res cut continuing" << std::endl;
      }     

      //#######################################
      //##### Calculate Shower Properties #####
      //#######################################

      //Now the shower metrics are ont he most energetic shower. 
      const art::Ptr<recob::Shower>  BiggestShower = neutrino_showers.at(0);

      //Get the shower properties 
      const TVector3            ShowerStart        = BiggestShower->ShowerStart();
      const TVector3            ShowerDirection    = BiggestShower->Direction();
      const int                 ShowerBestPlane    = BiggestShower->best_plane();
      const double              ShowerLength       = BiggestShower->Length();
      const double              ShowerOpeningAngle = BiggestShower->OpenAngle();
      const std::vector<double> ShowerEnergyVec    = BiggestShower->Energy();
      const std::vector<double> ShowerdEdxVec      = BiggestShower->dEdx();
      double ShowerEnergy = -999;
      double ShowerdEdx   = -999;

      if(ShowerEnergyVec.size() >= (ShowerBestPlane+1)){
	ShowerEnergy = ShowerEnergyVec.at(ShowerBestPlane);
      }
      if(ShowerdEdxVec.size() >= (ShowerBestPlane+1)){  
	ShowerdEdx = ShowerdEdxVec.at(ShowerBestPlane);     
      }

      //Get the showers 
      auto const &showerListHandle =						\
      	ev.getValidHandle<std::vector<recob::Shower> >(fConfig.ShowerTag);
      
      //Get The shower spacepoints 
      art::FindManyP<recob::Hit> fmh(showerListHandle, ev, fConfig.ShowerTag);
      if(!fmh.isValid()){
	throw cet::exception("NueSelectionReco") << "Hit-Recob::Shower association is somehow not valid. Stopping";
	return false;
      }

      //Association between Showers and 2d Hits
      art::FindManyP<recob::SpacePoint> fmsp(showerListHandle, ev, fConfig.ShowerTag);
      if(!fmsp.isValid()){
	throw cet::exception("NueSelectionReco") << "Spacepoint-Recob::Shower association is somehow not valid. Stopping";
	return false;
      }

      //Get the shower track stubs and hits
      art::FindManyP<recob::Track> fmst(showerListHandle, ev, fConfig.ShowerTag);
      if(!fmst.isValid()){
	throw cet::exception("NueSelectionReco") << "Track and Shower association is somehow not valid. Stopping";
	return false;
      }
      //and the hits
      auto const &showertrackHandle =					\
      	ev.getValidHandle<std::vector<recob::Track> >(fConfig.ShowerTag);
      art::FindManyP<recob::Hit> fmsth(showertrackHandle, ev, fConfig.ShowerTag);
      if(!fmsth.isValid()){
	throw cet::exception("NueSelectionReco") << "Track  and hit association is somehow not valid. Stopping";
	return false;
      }

      //Get the track hits to spacepoint assn
      auto const &HitHandle =					\
      	ev.getValidHandle<std::vector<recob::Hit> >(fConfig.HitTag);
      art::FindManyP<recob::SpacePoint> fmsphsp(HitHandle, ev,fConfig.PandoraTag);
      if(!fmsphsp.isValid()){
	throw cet::exception("NueSelectionReco") << "this one Spacepoint and hit association not valid. Stopping.";
	return false;
      }

      //Get the spacepoints handle and the hit assoication              
      auto const &spHandle = ev.getValidHandle<std::vector<recob::SpacePoint> >(fConfig.PandoraTag);
      art::FindManyP<recob::Hit> fmsph(spHandle, ev, fConfig.PandoraTag);
      if(!fmsph.isValid()){
	throw cet::exception("RecoEfficencyFinder") << "Spacepoint to hit association not valid. Stopping.";
	return false;
     }




      //Get this shower spacepoints
      std::vector<art::Ptr<recob::SpacePoint> > sps = fmsp.at(BiggestShower.key());
      
      //Calculate the shower density gradient. Scale for correct energy
      double pw = -1;

      double showerdensitygradient = ShowerDensityGradient(sps,ShowerStart,ShowerDirection,ShowerLength,ShowerOpeningAngle,fmsph,ShowerEnergy*1.145,pw);

      //Get the shower track 
      float ShowerTrackLength = -999;
      float ShowerTrackWidth   = -999; 
      std::vector<art::Ptr<recob::Track> > showertrack = fmst.at(BiggestShower.key());
      //Check the track was sucessful.
      if(showertrack.size()>0){
	
	//Get the Length
	ShowerTrackLength = showertrack.at(0)->Length();

	//Get the hits 
	std::vector<art::Ptr<recob::Hit> > showertrackhits = fmsth.at(showertrack.at(0).key());

	//Get the width
	ShowerTrackWidth = ShowerTrackWidthCal(showertrackhits,fmsphsp,ShowerStart,ShowerDirection);
      }

      //##################################
      //########## Shower Cuts ###########
      //##################################


      //Perform the conversion gap cut.
      double conversion_gap = -1;
      bool pass_conv = passConversionGapCut(ShowerStart,conversion_gap,vertex_position);
      if(!pass_conv && fConfig.ApplyConversionGapCut){
	if(fConfig.Verbose){
	  std::cout << "Failed the conversion gap cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_ConversionGapRemoved_Hist,fRootHists.VisibleEnergy_ConversionGapRemoved_HistMode,fRootHists.VisibleEnergy_ConversionGapRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_ConversionGapPassed_Hist,fRootHists.VisibleEnergy_ConversionGapPassed_HistMode,fRootHists.VisibleEnergy_ConversionGapPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the conversion gap cut" << std::endl;
      }

      //Perform the dE/dx Cut 
      bool pass_dEdx = passdEdxCut(ShowerdEdx);
      if(!pass_dEdx){
	if(fConfig.Verbose){
	  std::cout << "Failed the dEdx gap cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_dEdxRemoved_Hist,fRootHists.VisibleEnergy_dEdxRemoved_HistMode,fRootHists.VisibleEnergy_dEdxRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_dEdxPassed_Hist,fRootHists.VisibleEnergy_dEdxPassed_HistMode,fRootHists.VisibleEnergy_dEdxPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the dEdx cut" << std::endl;
      }

      //Peform the shower length cut which we never perform.
      bool pass_Length = passLengthCut(ShowerLength);
      if(!pass_Length){
	if(fConfig.Verbose){
	  std::cout << "Failed the Length  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_LengthRemoved_Hist,fRootHists.VisibleEnergy_LengthRemoved_HistMode,fRootHists.VisibleEnergy_LengthRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_LengthPassed_Hist,fRootHists.VisibleEnergy_LengthPassed_HistMode,fRootHists.VisibleEnergy_LengthPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the Length  cut" << std::endl;
      }
      

      //Peform the shower length cut which we never perform.
      bool pass_OpeningAngle = passOpeningAngleCut(ShowerOpeningAngle);
      if(!pass_OpeningAngle){
	if(fConfig.Verbose){
	  std::cout << "Failed the OpeningAngle  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_OpeningAngleRemoved_Hist,fRootHists.VisibleEnergy_OpeningAngleRemoved_HistMode,fRootHists.VisibleEnergy_OpeningAngleRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_OpeningAnglePassed_Hist,fRootHists.VisibleEnergy_OpeningAnglePassed_HistMode,fRootHists.VisibleEnergy_OpeningAnglePassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the OpeningAngle  cut" << std::endl;
      }

      //Perform the shower energy cut 
      bool pass_ShowerEnergy = passShowerEnergyCut(ShowerEnergy);
      if(!pass_ShowerEnergy){
	if(fConfig.Verbose){
	  std::cout << "Failed the ShowerEnergy  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_ShowerEnergyRemoved_Hist,fRootHists.VisibleEnergy_ShowerEnergyRemoved_HistMode,fRootHists.VisibleEnergy_ShowerEnergyRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_ShowerEnergyPassed_Hist,fRootHists.VisibleEnergy_ShowerEnergyPassed_HistMode,fRootHists.VisibleEnergy_ShowerEnergyPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the ShowerEnergy  cut" << std::endl;
      }

      //Perform the shower energy gradient cut
      bool pass_ShowerDensityGradient = passShowerDensityGradientCut(showerdensitygradient);
      if(!pass_ShowerDensityGradient){
	if(fConfig.Verbose){
	  std::cout << "Failed the ShowerDensityGradient  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_Hist,fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_HistMode,fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_ShowerDensityGradientPassed_Hist,fRootHists.VisibleEnergy_ShowerDensityGradientPassed_HistMode,fRootHists.VisibleEnergy_ShowerDensityGradientPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the ShowerDensityGradient  cut" << std::endl;
      }
      bool pass_ShowerDensityPower = passShowerDensityPowerCut(pw);
      if(!pass_ShowerDensityPower){
	if(fConfig.Verbose){
	  std::cout << "Failed the ShowerDensityPower  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_Hist,fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_HistMode,fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_ShowerDensityPowerPassed_Hist,fRootHists.VisibleEnergy_ShowerDensityPowerPassed_HistMode,fRootHists.VisibleEnergy_ShowerDensityPowerPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the ShowerDensityPower  cut" << std::endl;
      }
      bool pass_ShowerTrackLength = passShowerTrackLengthCut(ShowerTrackLength);
      if(!pass_ShowerTrackLength){
	if(fConfig.Verbose){
	  std::cout << "Failed the ShowerTrackLength  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_Hist,fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_HistMode,fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_ShowerTrackLengthPassed_Hist,fRootHists.VisibleEnergy_ShowerTrackLengthPassed_HistMode,fRootHists.VisibleEnergy_ShowerTrackLengthPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the ShowerTrackLength  cut" << std::endl;
      }
      bool pass_ShowerTrackWidth = passShowerTrackWidthCut(ShowerTrackWidth);
      if(!pass_ShowerTrackWidth){
	if(fConfig.Verbose){
	  std::cout << "Failed the ShowerTrackWidth  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_Hist,fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_HistMode,fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_ShowerTrackWidthPassed_Hist,fRootHists.VisibleEnergy_ShowerTrackWidthPassed_HistMode,fRootHists.VisibleEnergy_ShowerTrackWidthPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the ShowerTrackWidth  cut" << std::endl;
      }

      
      //##############################
      //######## Track Info ##########
      //##############################

      //Get the tracks
      auto const &trackListHandle =					\
      	ev.getValidHandle<std::vector<recob::Track> >(fConfig.TrackTag);
      art::FindManyP<recob::Track> fmt(pfpHandle, ev, fConfig.TrackTag);
      if(!fmt.isValid()){
	throw cet::exception(" NueSelectionReco") << "Track and PF particle association is somehow not valid. Stopping";
	return false;
      }

      //Get the pida information
      art::FindManyP<anab::ParticleID> fmpid(trackListHandle, ev, fConfig.TrackPIDTag);
      if(!fmpid.isValid()){
	throw cet::exception("RecoEfficencyFinder") << "Track and PID association is somehow not valid. Stopping";
	return false;
      }


      double MaxTrackLength = -999;
      double MaxTrackPIDA   = -999;

      //Identify the largest track in the event
      for(auto const& daughter: daughters){
	
      	//is the daughter a shower.
      	if(pfp_map[daughter]->PdgCode() != 13){continue;}

	//Get the track info 
	std::vector<art::Ptr<recob::Track> > daughter_track = fmt.at(pfp_map[daughter].key());
	if(daughter_track.size() == 0){
	  std::cout << "No Track" << std::endl;
	  continue;
	}
	if(daughter_track.size() != 1){
	  throw cet::exception("NueSelectionReco") << "we have too many recob tracks for pfparticles";
	return false;
	}
	
	if(daughter_track.at(0)->Length() > MaxTrackLength){
	  MaxTrackLength = daughter_track.at(0)->Length();

	  //Get the the PID information 
	  std::vector<art::Ptr<anab::ParticleID > > pids = fmpid.at(daughter_track[0].key());
	  float PIDA = -999;
	  if(pids.size() > 0){
	    PIDA = pids.at(pids.size()-1)->PIDA();
	  }
	  MaxTrackPIDA   = PIDA;
	}
      }

      //##############################
      //######## Track Cuts ##########
      //##############################

      bool pass_MaxTrackLength = passMaxTrackLengthCut(MaxTrackLength);
      if(!pass_MaxTrackLength){
	if(fConfig.Verbose){
	  std::cout << "Failed the MaxTrackLength  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_MaxTrackLengthRemoved_Hist,fRootHists.VisibleEnergy_MaxTrackLengthRemoved_HistMode,fRootHists.VisibleEnergy_MaxTrackLengthRemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_MaxTrackLengthPassed_Hist,fRootHists.VisibleEnergy_MaxTrackLengthPassed_HistMode,fRootHists.VisibleEnergy_MaxTrackLengthPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the MaxTrackLength  cut" << std::endl;
      }
      bool pass_MaxTrackPIDA = passMaxTrackPIDACut(MaxTrackPIDA);
      if(!pass_MaxTrackPIDA){
	if(fConfig.Verbose){
	  std::cout << "Failed the MaxTrackPIDA  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_MaxTrackPIDARemoved_Hist,fRootHists.VisibleEnergy_MaxTrackPIDARemoved_HistMode,fRootHists.VisibleEnergy_MaxTrackPIDARemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_MaxTrackPIDAPassed_Hist,fRootHists.VisibleEnergy_MaxTrackPIDAPassed_HistMode,fRootHists.VisibleEnergy_MaxTrackPIDAPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the MaxTrackPIDA  cut" << std::endl;
      }
      bool pass_MaxTrackLengthPIDA = passMaxTrackLengthPIDACut(MaxTrackLength,MaxTrackPIDA);
      if(!pass_MaxTrackLengthPIDA){
	if(fConfig.Verbose){
	  std::cout << "Failed the MaxTrackLengthPIDA  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_Hist,fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_HistMode,fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_HistIntType,intInfo);
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_Hist,fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_HistMode,fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_HistIntType,intInfo);
      if(fConfig.Verbose){
	std::cout << "Passed the MaxTrackLengthPIDA  cut" << std::endl;
      }

      
      //#########################
      //######## MVA Cut ########
      //#########################

      //Using values above pass to the BDT 
      if(fConfig.SetupMVA){
	if(MVAValues.find("dEdxMVA") != MVAValues.end()){
	  MVAValues["dEdxMVA"] = ShowerdEdx;
	}
	if(MVAValues.find("CoversionGapMVA") != MVAValues.end()){
	  MVAValues["CoversionGapMVA"] = conversion_gap;
	}
	if(MVAValues.find("MaxTrackLenghMVA") != MVAValues.end()){
	  MVAValues["MaxTrackLenghMVA"] = MaxTrackLength;
	}
	if(MVAValues.find("MaxTrackPIDAMVA") != MVAValues.end()){
	  MVAValues["MaxTrackPIDAMVA"] = MaxTrackPIDA;
	}
	if(MVAValues.find("ShowerEnergyMVA") != MVAValues.end()){
	  MVAValues["ShowerEnergyMVA"] = ShowerEnergy;
	}
	if(MVAValues.find("ShowerLengthEMVA") != MVAValues.end()){
	  MVAValues["ShowerLengthEMVA"] = ShowerLength;
	}
	if(MVAValues.find("ShowerDensityOpeningAngleMVA") != MVAValues.end()){
	  MVAValues["ShowerDensityOpeningAngleMVA"] = ShowerOpeningAngle;
	}
	if(MVAValues.find("NumberOfNeutrinosMVA") != MVAValues.end()){
	  MVAValues["NumberOfNeutrinosMVA"] = numneutrinos;
	}
	if(MVAValues.find("ShowerResidualNumShowersMVA") != MVAValues.end()){
	  MVAValues["ShowerResidualNumShowersMVA"] = numshowers;
	}
	if(MVAValues.find("ShowerDensityPWMVA") != MVAValues.end()){
	  MVAValues["ShowerDensityPWMVA"] = pw;
	}
	if(MVAValues.find("ShowerDensityGradNewMVA") != MVAValues.end()){
	  MVAValues["ShowerDensityGradNewMVA"] = showerdensitygradient;
	}
	if(MVAValues.find("NuetrinoPdGMVA") != MVAValues.end()){
	  MVAValues["NuetrinoPdGMVA"] = neutrinopdg;
	}
	if(MVAValues.find("ShowerTrackWidthMVA") != MVAValues.end()){
	  MVAValues["ShowerTrackWidthMVA"] = ShowerTrackWidth;
	}
	if(MVAValues.find("ShowerTrackLengthMVA") !=  MVAValues.end()){
	  MVAValues["ShowerTrackLengthMVA"] = ShowerTrackLength;
	}
	
	bool pass_MVA = passMVACut();
	if(!pass_MVA){
	if(fConfig.Verbose){
	  std::cout << "Failed the MVA  cut" << std::endl;
	}	
	FillHistograms(fRootHists.VisibleEnergy_MVARemoved_Hist,fRootHists.VisibleEnergy_MVARemoved_HistMode,fRootHists.VisibleEnergy_MVARemoved_HistIntType,intInfo);
	return false;
	}
	FillHistograms(fRootHists.VisibleEnergy_MVAPassed_Hist,fRootHists.VisibleEnergy_MVAPassed_HistMode,fRootHists.VisibleEnergy_MVAPassed_HistIntType,intInfo);
	if(fConfig.Verbose){
	  std::cout << "Passed the MVA  cut" << std::endl;
	}
      }
      
      std::cout << "Event has been selected" << std::endl;

      //You can passed all the cuts 
      return true;
    }
    
    //Check if the point is the fiducal volume.
    bool NueSelectionReco::containedInFV(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& FV: fConfig.fiducial_volumes) {
	if (FV.Contain(p)) return true;
      }
      return false;
    }

    //Check if the point is in the Active volume.
    bool NueSelectionReco::containedInAV(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& AV: fConfig.active_volumes) {
	if (AV.Contain(p)) return true;
      }
      return false;
    }

    void NueSelectionReco::GetSliceHits(const art::Ptr<recob::PFParticle>& pfp, 
					std::vector< art::Ptr< recob::Hit> >& pfpHits,
					std::map<int, art::Ptr<recob::PFParticle> >& pfpMap,
					art::FindManyP<recob::Cluster>& fmpfc,
					art::FindManyP<recob::Hit>& fmch){


      // Get the hits from the PFParticle
      const std::vector< art::Ptr< recob::Cluster> >& clusters = fmpfc.at(pfp.key());
      for (const auto& cluster: clusters){
	const std::vector< art::Ptr< recob::Hit> >& hits = fmch.at(cluster.key());
	pfpHits.insert(pfpHits.end(), hits.begin(), hits.end());
      }
      
      // Get the daughters
      const std::vector<long unsigned int> daughters = pfp->Daughters();
      for (const auto daughterIter: daughters){
	art::Ptr<recob::PFParticle> daughter =  pfpMap.at(daughterIter);
	// Get the hits from the daughters
	GetSliceHits(daughter, pfpHits, pfpMap, fmpfc,fmch);
      }
  

      return;
    }

    double NueSelectionReco::ShowerDensityGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, 
						   const TVector3& ShowerStartPosition, 
						   const TVector3& ShowerDirection,
						   const double ShowerLength, 
						   const double& OpenAngle, 
						   art::FindManyP<recob::Hit> const& fmh, 
						   const double ShowerEnergy, 
						   double& pw){

      std::map<int, std::vector<art::Ptr<recob::SpacePoint> > > len_segment_map;
      double segmentsize = ShowerLength/fConfig.NSegments;
      
      //Split the the spacepoints into segments.
      for(auto const& sp: sps){
    
	//Get the position of the spacepoint
	TVector3 pos = SpacePointPosition(sp) - ShowerStartPosition;
    
	//Get the the projected length
	double len = pos.Dot(ShowerDirection);
	
	//Get the length to the projection
	TVector3 perp = pos - (len*ShowerDirection);
	double perpLen = perp.Mag();

	//Get where the sp should be place.
	int sg_len = round(len/segmentsize);

	//Only add if the hit within the cone
	if(perpLen > TMath::Abs(TMath::Tan(OpenAngle)*len)){
          continue;
	}
	
	len_segment_map[sg_len].push_back(sp);
      }

      TGraph* graph = new TGraph();

      //Calculate the density gradent.
      for(auto& segment: len_segment_map){
	double sg_len = segment.first;
    
	if(segment.second.size() < 10){continue;}

	//Calculate the charge in the segement
	double SegmentEnergy = TotalEnergy(segment.second,fmh);
    
	//Calculate the voume
	double lower_dist = sg_len*segmentsize - segmentsize/2;
	double upper_dist = sg_len*segmentsize + segmentsize/2;
	
	//if(fRemoveStartFin){if(sg_len==0 || sg_len==fNSegments){continue;}}

	if(sg_len==0)         {lower_dist = 0;}
	if(sg_len==fConfig.NSegments){upper_dist = sg_len*segmentsize;}
    
	double littlevolume = lower_dist*TMath::Power((TMath::Tan(0.5*OpenAngle)*lower_dist),2)*TMath::Pi()/3;
	double bigvolume    = upper_dist*TMath::Power((TMath::Tan(0.5*OpenAngle)*upper_dist),2)*TMath::Pi()/3;
	double volume       = bigvolume - littlevolume;
    
	double SegmentDensity = SegmentEnergy/volume;

	double LengthToSegment = (lower_dist+upper_dist)/2;

	graph->SetPoint(graph->GetN(),LengthToSegment,SegmentDensity/ShowerEnergy);
      }
  
      if(graph->GetN() < 3 ){return -999;}
      
      TF1 *fit = new TF1("fit", "[0]/x^[1]");
      fit->SetParLimits(1,1,2);
      fit->SetParLimits(0,0,1);

      graph->Fit(fit,"Q");
      
      double grad = fit->GetParameter(0);
      pw   = fit->GetParameter(1);
  

      delete fit;
      delete graph;
      return grad;
    }
    
    double  NueSelectionReco::TotalEnergy(const std::vector<art::Ptr<recob::SpacePoint> >&sps, 
					  art::FindManyP<recob::Hit> const& fmh){
  
      std::map<int,double> TotalCharge_plane;
      std::map<int,int> TotalCharge_int;
      for(auto const& sp: sps){
	double Charge =  SpacePointCharge(sp,fmh);
	double Time   =  SpacePointTime(sp,fmh);
	int    Plane  =  SpacePointPlane(sp,fmh);
	Charge *= TMath::Exp((fDetProp->SamplingRate() * Time ) / (fDetProp->ElectronLifetime()*1e3));
	TotalCharge_plane[Plane] += Charge;
	++TotalCharge_int[Plane];
      }

      int max_plane = -999;
      for(auto const& plane: TotalCharge_int){
	if(plane.second > max_plane){
	  max_plane = plane.first;
	}
      }
      if(max_plane == -999){return 0;}
    
      double TotalEnergy = TotalCharge_plane[max_plane]*fConfig.EnergyConversion[max_plane];

      return TotalEnergy; 
    }
    
    int  NueSelectionReco::SpacePointPlane(art::Ptr<recob::SpacePoint> const& sp,
					   art::FindManyP<recob::Hit> const& fmh) const {

      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());

      if(hits.size() != 1){
	throw cet::exception("RecoEfficencyFinder") << "Spacepoint is not matched to 1 hit. This was unexpected";
      }

      return hits[0]->WireID().Plane;

    }

    double  NueSelectionReco::SpacePointCharge(art::Ptr<recob::SpacePoint> const& sp,
					       art::FindManyP<recob::Hit> const& fmh) const {
  
      double Charge = 0;
  
      //Average over the charge even though there is only one
      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());
      for(auto const& hit: hits){
	Charge += hit->Integral();
      }
  
      Charge /= (float) hits.size();
  
      return Charge;
    }

    //Return the spacepoint time.
    double  NueSelectionReco::SpacePointTime(art::Ptr<recob::SpacePoint> const& sp,
					     art::FindManyP<recob::Hit> const& fmh) const {
  
      double Time = 0;

      //Avergae over the hits
      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());
      for(auto const& hit: hits){
	Time += hit->PeakTime();
      }

      Time /= (float) hits.size();
      return Time;
    }


    double  NueSelectionReco::SpacePointProjection(const art::Ptr<recob::SpacePoint>&sp,
						   TVector3 const& vertex, 
						   TVector3 const& direction) const {

      // Get the position of the spacepoint
      TVector3 pos = SpacePointPosition(sp) - vertex;

      // Get the the projected length
      double projLen = pos.Dot(direction);

      return projLen;
    }

    TVector3 NueSelectionReco::SpacePointPosition(art::Ptr<recob::SpacePoint> const& sp) const {

      const Double32_t* sp_xyz = sp->XYZ();
      TVector3 sp_postiion = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
      return sp_postiion;
    }

    double NueSelectionReco::ShowerTrackWidthCal(std::vector<art::Ptr<recob::Hit> > const& showertrackhits
						 ,art::FindManyP<recob::SpacePoint> const& fmsphsp,
						 const TVector3& ShowerStart,
						 const TVector3& ShowerDirection) const{

      //Calculate the average spread from the shower direction
      float Perp = 0;
      for(auto const& showertrackhit: showertrackhits){
	  
	//Get the spacepoint. 
	std::vector<art::Ptr<recob::SpacePoint> > sps = fmsphsp.at(showertrackhit.key());
	for(auto const& sp: sps){
	    
	  //Calulate the perpendicular distance 
	  const Double32_t* sp_xyz = sp->XYZ();
	  TVector3 sp_position = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
	  TVector3 pos = sp_position - ShowerStart;
	  double  proj = pos.Dot(ShowerDirection);
	  pos = pos - proj * ShowerDirection;
	  Perp += pos.Mag(); 
	}
      }
      Perp= Perp/(int)showertrackhits.size();
      return Perp;
    }
  



    bool NueSelectionReco::OneShowerEnergyCut(std::vector<art::Ptr<recob::Shower> > const& neutrino_showers) const{
      //Remove event if thee is more than x showers above a certain energy threshold
      if(neutrino_showers.size() > fConfig.NumShowersCut){
	//Check that secondary and further showers have less energy
	const art::Ptr<recob::Shower> SecondaryShower = neutrino_showers.at(fConfig.NumShowersCut);
	if(SecondaryShower->Energy().at(SecondaryShower->best_plane()) > fConfig.SecondaryShowerEnergyCut){
	  if(fConfig.Verbose){std::cout << "too many showers one shower cut. The next biggest shower has energy: " << SecondaryShower->Energy().at(SecondaryShower->best_plane()) << std::endl;}
	return false;
	}
      }
      return true;
    }


    bool NueSelectionReco::OneShowerResiudalCut(std::vector<art::Ptr<recob::Shower> > const& neutrino_showers,
						double& numshowers, 
						const TVector3& vertex_position){
      
      //Find the distance to the closest shower
      if(neutrino_showers.size() < 0){return true;}

      const TVector3 biggestshower_start     = neutrino_showers[0]->ShowerStart(); 
      const TVector3 biggestshower_direction = neutrino_showers[0]->Direction(); 
      const double shower_openingangle       = neutrino_showers[0]->OpenAngle();

      //Identify the resiudal distance and check if shower is within cone
      int shower_num = 0;  
      for(auto const& neutrino_shower: neutrino_showers){
	
	const TVector3 ShowerStart     = neutrino_shower->ShowerStart();//cm 
	TVector3 conversion_vec  = (ShowerStart - biggestshower_start);
	
		//Calculate the projected length and perpendicualr length
 	double len  = conversion_vec.Dot(biggestshower_direction);
	double perp = (conversion_vec - len*biggestshower_direction).Mag();

	//Calculate the cone length of the biggest shower for that shower 
	double Shower_dist = len*TMath::Tan(0.5*shower_openingangle);

	//Get the Best Shower energy pass the error value if these things have not been correctly made.
	double ShowerEnergy = -999;
	const std::vector<double> ShowerEnergyPlanes = neutrino_shower->Energy(); //MeV

	const int ShowerBest_Plane = neutrino_shower->best_plane();
	if(ShowerBest_Plane != -999){
	  if((int) ShowerEnergyPlanes.size() >= ShowerBest_Plane+1){
	    ShowerEnergy = ShowerEnergyPlanes[ShowerBest_Plane];
	  }
	}
	else{
	  std::cerr << "Dom Has not considered the energy not being reconstructed in the resiudal cut please check."<< std::endl;
	}

	if(fConfig.Verbose){
	  std::cout << " shower has energy: " << ShowerEnergy << " and residual: " << perp;
	}

	//Count the number of showers above the resiudal cut allow for 3 cm off
	if(ShowerEnergy > fConfig.ShowerResidualEnergyCut && perp > TMath::Abs(fConfig.ResiudalCut*Shower_dist)+2){
	  ++shower_num;
	}
      }

      if(fConfig.Verbose){
	std::cout << "number of showers above the cuts: " << shower_num << std::endl;
      }

      numshowers = shower_num;
      if(shower_num > 0){ 
	return false;
      }
      
      return true;
    }	  

    bool NueSelectionReco::ConversionGapCut(const TVector3 ShowerStart, 
					    double& conversion_gap,
					    const TVector3& vertex_position){

      if(fConfig.Verbose){
	std::cout << "Conversion gap is: " << (ShowerStart-vertex_position).Mag()  << " for a Cut of: " << fConfig.ConversionGapCut << std::endl;
      }
      
      
      if((ShowerStart-vertex_position).Mag() > fConfig.ConversionGapCut){
	return false;
      }
      return true;
    }

    bool NueSelectionReco::NeutrinoPdgCodeCut(const double& NeutrinoPdgCode) const{

      if(fConfig.Verbose){std::cout << "NeutrinoPdgCode is: " << NeutrinoPdgCode << std::endl;}
      if(NeutrinoPdgCode > fConfig.NeutrinoPdgCodeCut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::NumneutrinosCut(const double& Numneutrinos) const{

      if(fConfig.Verbose){std::cout << "Numneutrinos is: " << Numneutrinos << std::endl;}
      if(Numneutrinos > fConfig.NumneutrinosCut){
	return false;
      }
      return true;
    }

    bool NueSelectionReco::dEdxCut(const double& dEdx) const{

      if(fConfig.Verbose){std::cout << "dEdx is: " << dEdx << std::endl;}
      if(dEdx > fConfig.dEdxCut){
	return false;
      }
      return true;
    }

    bool NueSelectionReco::LengthCut(const double& Length) const{

      if(fConfig.Verbose){std::cout << "Length is: " << Length << std::endl;}
      if(Length < fConfig.LengthCut){
	return false;
      }
      return true;
    }

    bool NueSelectionReco::OpeningAngleCut(const double& OpeningAngle) const{

      if(fConfig.Verbose){std::cout << "OpeningAngle is: " << OpeningAngle << std::endl;}
      if(OpeningAngle < fConfig.OpeningAngleCut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::ShowerEnergyCut(const double& ShowerEnergy) const{

      if(fConfig.Verbose){std::cout << "ShowerEnergy is: " << ShowerEnergy << std::endl;}
      if(ShowerEnergy < fConfig.ShowerEnergyCut){
	return false;
      }
      return true;
    }

    bool NueSelectionReco::ShowerDensityGradientCut(const double& ShowerDensityGradient) const{

      if(fConfig.Verbose){std::cout << "ShowerDensityGradient is: " << ShowerDensityGradient << std::endl;}
      if(ShowerDensityGradient > fConfig.ShowerDensityGradientCut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::ShowerDensityPowerCut(const double& ShowerDensityPower) const{

      if(fConfig.Verbose){std::cout << "ShowerDensityPower is: " << ShowerDensityPower << std::endl;}
      if(ShowerDensityPower > fConfig.ShowerDensityPowerCut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::ShowerTrackLengthCut(const double& ShowerTrackLength) const{

      if(fConfig.Verbose){std::cout << "ShowerTrackLength is: " << ShowerTrackLength << std::endl;}
      if(ShowerTrackLength < fConfig.ShowerTrackLengthCut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::ShowerTrackWidthCut(const double& ShowerTrackWidth) const{

      if(fConfig.Verbose){std::cout << "ShowerTrackWidth is: " << ShowerTrackWidth << std::endl;}
      if(ShowerTrackWidth < fConfig.ShowerTrackWidthCut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::MaxTrackLengthCut(const double& MaxTrackLength) const{

      if(fConfig.Verbose){std::cout << "MaxTrackLength is: " << MaxTrackLength << std::endl;}
      if(MaxTrackLength > fConfig.MaxTrackLengthCut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::MaxTrackPIDACut(const double& MaxTrackPIDA) const{

      if(fConfig.Verbose){std::cout << "MaxTrackPIDA is: " << MaxTrackPIDA << std::endl;}
      if(MaxTrackPIDA > fConfig.MaxTrackPIDACut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::MaxTrackLengthPIDACut(const double& MaxTrackLength,const double& MaxTrackPIDA) const{

      if(fConfig.Verbose){std::cout << "MaxTrackLength is: " << MaxTrackLength << " and PIDA is: " << MaxTrackPIDA << std::endl;}
      if(MaxTrackLength > fConfig.MaxTrackLengthPIDACut && MaxTrackPIDA < fConfig.MaxTrackPIDACut){
	return false;
      }
      return true;
    }
    bool NueSelectionReco::MVACut() const{

      if(fConfig.Verbose){std::cout << "MVA Value is: " << reader->EvaluateMVA(fConfig.MVAMethod) << std::endl;}

      //Perform the BDT Cut
      if(reader->EvaluateMVA(fConfig.MVAMethod) < fConfig.MVACut){
	return false;
      }
      return true;
    }


	
  
    //Fill the histograms.... can't think of a better comment.
    void NueSelectionReco::FillHistograms(std::map<std::string,TH1D*>& HistMap, std::map<std::string,std::map<int,TH1D*> >& HistMapMode, std::map<std::string,std::map<int,TH1D*> >& HistMapIntType, NueSelectionReco::NueInteraction& intInfo){
      
      double Energy = intInfo.True_energy;
      
      FillHistograms(HistMap,HistMapMode,HistMapIntType,intInfo.mctruth->GetNeutrino(),intInfo,Energy,intInfo.booldirtevent);

    }

    void NueSelectionReco::FillHistograms(std::map<std::string,TH1D*>& HistMap, std::map<std::string,std::map<int,TH1D*> >& HistMapMode, std::map<std::string,std::map<int,TH1D*> >& HistMapIntType, const simb::MCNeutrino& nu, NueSelectionReco::NueInteraction& intInfo,double Energy, bool& booldirtevent){
      
      if(!fConfig.FillHistograms){return;}

      if(nu.Nu().PdgCode() != intInfo.initnu && nu.Nu().PdgCode() == 11 && nu.CCNC() == 0  && dirtevent == false){
	HistMap["AllSignal"]->Fill(Energy,intInfo.weight);
      }
      else{
	HistMap["AllBackground"]->Fill(Energy,intInfo.weight);
      }


      int mode = -1;
      if(nu.Mode() > -1 && nu.Mode() < 6){mode = nu.Mode();}

      if(booldirtevent == true){
	if(nu.CCNC() == simb::kCC){

	  //Check to see if we are Numu background
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["DirtNuMu"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtNuMu"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)    HistMapMode["DirtNuMu"][mode]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["DirtInNuE"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtInNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)    HistMapMode["DirtInNuE"][mode]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["DirtOscNuE"]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtOscNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillModeHistograms)    HistMapMode["DirtOscNuE"][mode]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
	else{ 
	  //Neutral current histogtam fill
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["DirtNCNuMu"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtNCNuMu"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)    HistMapMode["DirtNCNuMu"][mode]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["DirtNCInNuE"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtNCInNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)    HistMapMode["DirtNCInNuE"][mode]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["DirtNCOscNuE"]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillIntTypeHistograms) HistMapIntType["DirtNCOscNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillModeHistograms)     HistMapMode["DirtNCOscNuE"][mode]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
      }
      else{
	//Check if we are charged current
	if(nu.CCNC() == simb::kCC){
	  //Check to see if we are Numu background
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["NuMu"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["NuMu"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)     HistMapMode["NuMu"][mode]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["InNuE"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["InNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)     HistMapMode["InNuE"][mode]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["OscNuE"]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillIntTypeHistograms) HistMapIntType["OscNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillModeHistograms)     HistMapMode["OscNuE"][mode]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
	else{ 
	  //Neutral current histogtam fill
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["NCNuMu"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["NCNuMu"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)     HistMapMode["NCNuMu"][mode]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["NCInNuE"]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillIntTypeHistograms) HistMapIntType["NCInNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	    if(fConfig.FillModeHistograms)     HistMapMode["NCInNuE"][mode]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["NCOscNuE"]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillIntTypeHistograms) HistMapIntType["NCOscNuE"][nu.InteractionType()]->Fill(Energy,intInfo.weight);
	      if(fConfig.FillModeHistograms)     HistMapMode["NCOscNuE"][mode]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
      }
      return;
    }


    void NueSelectionReco::InitialiseHistograms(){

      fRootHists.XDiff = new TH1D("XDiff","XDiff",200,-200,200);
      fRootHists.YDiff = new TH1D("YDiff","YDiff",200,-200,200);
      fRootHists.ZDiff = new TH1D("ZDiff","ZDiff",200,-200,200);

      //Loop over the types of interactions and set the bins
      if(fConfig.FillHistograms){
	for(auto& Type: fRootHists.HistTypes){

	  std::string  TrueNumber_String                     = Type + " TrueNumber";
	  const char* TrueNumber_Name                     = TrueNumber_String.c_str();
	  fRootHists.TrueNumber_Hist[Type] = new TH1D(TrueNumber_Name,TrueNumber_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_BeforeSel_String                     = Type + " VisibleEnergy_BeforeSel";
	  const char* VisibleEnergy_BeforeSel_Name                     = VisibleEnergy_BeforeSel_String.c_str();
	  fRootHists.VisibleEnergy_BeforeSel_Hist[Type] = new TH1D(VisibleEnergy_BeforeSel_Name,VisibleEnergy_BeforeSel_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_AfterSel_String                     = Type + " VisibleEnergy_AfterSel";
	  const char* VisibleEnergy_AfterSel_Name                     = VisibleEnergy_AfterSel_String.c_str();
	  fRootHists.VisibleEnergy_AfterSel_Hist[Type] = new TH1D(VisibleEnergy_AfterSel_Name,VisibleEnergy_AfterSel_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_AfterSelOne_String                     = Type + " VisibleEnergy_AfterSelOne";
	  const char* VisibleEnergy_AfterSelOne_Name                     = VisibleEnergy_AfterSelOne_String.c_str();
	  fRootHists.VisibleEnergy_AfterSelOne_Hist[Type] = new TH1D(VisibleEnergy_AfterSelOne_Name,VisibleEnergy_AfterSelOne_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_AfterSelExtra_String                     = Type + " VisibleEnergy_AfterSelExtra";
	  const char* VisibleEnergy_AfterSelExtra_Name                     = VisibleEnergy_AfterSelExtra_String.c_str();
	  fRootHists.VisibleEnergy_AfterSelExtra_Hist[Type] = new TH1D(VisibleEnergy_AfterSelExtra_Name,VisibleEnergy_AfterSelExtra_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	  std::string  VisibleEnergy_MCTruth_String                     = Type + " VisibleEnergy_MCTruth";
	  const char* VisibleEnergy_MCTruth_Name                     = VisibleEnergy_MCTruth_String.c_str();
	  fRootHists.VisibleEnergy_MCTruth_Hist[Type] = new TH1D(VisibleEnergy_MCTruth_Name,VisibleEnergy_MCTruth_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);



	  std::string  VisibleEnergy_FVRemoved_String                     = Type + " VisibleEnergy_FVRemoved";
	  const char* VisibleEnergy_FVRemoved_Name                     = VisibleEnergy_FVRemoved_String.c_str();
	  fRootHists.VisibleEnergy_FVRemoved_Hist[Type] = new TH1D(VisibleEnergy_FVRemoved_Name,VisibleEnergy_FVRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_FVPassed_String                     = Type + " VisibleEnergy_FVPassed";
	  const char* VisibleEnergy_FVPassed_Name                     = VisibleEnergy_FVPassed_String.c_str();
	  fRootHists.VisibleEnergy_FVPassed_Hist[Type] = new TH1D(VisibleEnergy_FVPassed_Name,VisibleEnergy_FVPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_GtrOneShowerRemoved_String                     = Type + " VisibleEnergy_GtrOneShowerRemoved";
	  const char* VisibleEnergy_GtrOneShowerRemoved_Name                     = VisibleEnergy_GtrOneShowerRemoved_String.c_str();
	  fRootHists.VisibleEnergy_GtrOneShowerRemoved_Hist[Type] = new TH1D(VisibleEnergy_GtrOneShowerRemoved_Name,VisibleEnergy_GtrOneShowerRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_GtrOneShowerPassed_String                     = Type + " VisibleEnergy_GtrOneShowerPassed";
	  const char* VisibleEnergy_GtrOneShowerPassed_Name                     = VisibleEnergy_GtrOneShowerPassed_String.c_str();
	  fRootHists.VisibleEnergy_GtrOneShowerPassed_Hist[Type] = new TH1D(VisibleEnergy_GtrOneShowerPassed_Name,VisibleEnergy_GtrOneShowerPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_OneShowerECutRemoved_String                     = Type + " VisibleEnergy_OneShowerECutRemoved";
	  const char* VisibleEnergy_OneShowerECutRemoved_Name                     = VisibleEnergy_OneShowerECutRemoved_String.c_str();
	  fRootHists.VisibleEnergy_OneShowerECutRemoved_Hist[Type] = new TH1D(VisibleEnergy_OneShowerECutRemoved_Name,VisibleEnergy_OneShowerECutRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_OneShowerECutPassed_String                     = Type + " VisibleEnergy_OneShowerECutPassed";
	  const char* VisibleEnergy_OneShowerECutPassed_Name                     = VisibleEnergy_OneShowerECutPassed_String.c_str();
	  fRootHists.VisibleEnergy_OneShowerECutPassed_Hist[Type] = new TH1D(VisibleEnergy_OneShowerECutPassed_Name,VisibleEnergy_OneShowerECutPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_OneShowerResidualRemoved_String                     = Type + " VisibleEnergy_OneShowerResidualRemoved";
	  const char* VisibleEnergy_OneShowerResidualRemoved_Name                     = VisibleEnergy_OneShowerResidualRemoved_String.c_str();
	  fRootHists.VisibleEnergy_OneShowerResidualRemoved_Hist[Type] = new TH1D(VisibleEnergy_OneShowerResidualRemoved_Name,VisibleEnergy_OneShowerResidualRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_OneShowerResidualPassed_String                     = Type + " VisibleEnergy_OneShowerResidualPassed";
	  const char* VisibleEnergy_OneShowerResidualPassed_Name                     = VisibleEnergy_OneShowerResidualPassed_String.c_str();
	  fRootHists.VisibleEnergy_OneShowerResidualPassed_Hist[Type] = new TH1D(VisibleEnergy_OneShowerResidualPassed_Name,VisibleEnergy_OneShowerResidualPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_ConversionGapRemoved_String                     = Type + " VisibleEnergy_ConversionGapRemoved";
	  const char* VisibleEnergy_ConversionGapRemoved_Name                     = VisibleEnergy_ConversionGapRemoved_String.c_str();
	  fRootHists.VisibleEnergy_ConversionGapRemoved_Hist[Type] = new TH1D(VisibleEnergy_ConversionGapRemoved_Name,VisibleEnergy_ConversionGapRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_ConversionGapPassed_String                     = Type + " VisibleEnergy_ConversionGapPassed";
	  const char* VisibleEnergy_ConversionGapPassed_Name                     = VisibleEnergy_ConversionGapPassed_String.c_str();
	  fRootHists.VisibleEnergy_ConversionGapPassed_Hist[Type] = new TH1D(VisibleEnergy_ConversionGapPassed_Name,VisibleEnergy_ConversionGapPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_dEdxRemoved_String                     = Type + " VisibleEnergy_dEdxRemoved";
	  const char* VisibleEnergy_dEdxRemoved_Name                     = VisibleEnergy_dEdxRemoved_String.c_str();
	  fRootHists.VisibleEnergy_dEdxRemoved_Hist[Type] = new TH1D(VisibleEnergy_dEdxRemoved_Name,VisibleEnergy_dEdxRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_dEdxPassed_String                     = Type + " VisibleEnergy_dEdxPassed";
	  const char* VisibleEnergy_dEdxPassed_Name                     = VisibleEnergy_dEdxPassed_String.c_str();
	  fRootHists.VisibleEnergy_dEdxPassed_Hist[Type] = new TH1D(VisibleEnergy_dEdxPassed_Name,VisibleEnergy_dEdxPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_LengthRemoved_String                     = Type + " VisibleEnergy_LengthRemoved";
	  const char* VisibleEnergy_LengthRemoved_Name                     = VisibleEnergy_LengthRemoved_String.c_str();
	  fRootHists.VisibleEnergy_LengthRemoved_Hist[Type] = new TH1D(VisibleEnergy_LengthRemoved_Name,VisibleEnergy_LengthRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_LengthPassed_String                     = Type + " VisibleEnergy_LengthPassed";
	  const char* VisibleEnergy_LengthPassed_Name                     = VisibleEnergy_LengthPassed_String.c_str();
	  fRootHists.VisibleEnergy_LengthPassed_Hist[Type] = new TH1D(VisibleEnergy_LengthPassed_Name,VisibleEnergy_LengthPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_OpeningAngleRemoved_String                     = Type + " VisibleEnergy_OpeningAngleRemoved";
	  const char* VisibleEnergy_OpeningAngleRemoved_Name                     = VisibleEnergy_OpeningAngleRemoved_String.c_str();
	  fRootHists.VisibleEnergy_OpeningAngleRemoved_Hist[Type] = new TH1D(VisibleEnergy_OpeningAngleRemoved_Name,VisibleEnergy_OpeningAngleRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_OpeningAnglePassed_String                     = Type + " VisibleEnergy_OpeningAnglePassed";
	  const char* VisibleEnergy_OpeningAnglePassed_Name                     = VisibleEnergy_OpeningAnglePassed_String.c_str();
	  fRootHists.VisibleEnergy_OpeningAnglePassed_Hist[Type] = new TH1D(VisibleEnergy_OpeningAnglePassed_Name,VisibleEnergy_OpeningAnglePassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_ShowerEnergyRemoved_String                     = Type + " VisibleEnergy_ShowerEnergyRemoved";
	  const char* VisibleEnergy_ShowerEnergyRemoved_Name                     = VisibleEnergy_ShowerEnergyRemoved_String.c_str();
	  fRootHists.VisibleEnergy_ShowerEnergyRemoved_Hist[Type] = new TH1D(VisibleEnergy_ShowerEnergyRemoved_Name,VisibleEnergy_ShowerEnergyRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_ShowerEnergyPassed_String                     = Type + " VisibleEnergy_ShowerEnergyPassed";
	  const char* VisibleEnergy_ShowerEnergyPassed_Name                     = VisibleEnergy_ShowerEnergyPassed_String.c_str();
	  fRootHists.VisibleEnergy_ShowerEnergyPassed_Hist[Type] = new TH1D(VisibleEnergy_ShowerEnergyPassed_Name,VisibleEnergy_ShowerEnergyPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_ShowerDensityGradientRemoved_String                     = Type + " VisibleEnergy_ShowerDensityGradientRemoved";
	  const char* VisibleEnergy_ShowerDensityGradientRemoved_Name                     = VisibleEnergy_ShowerDensityGradientRemoved_String.c_str();
	  fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_Hist[Type] = new TH1D(VisibleEnergy_ShowerDensityGradientRemoved_Name,VisibleEnergy_ShowerDensityGradientRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_ShowerDensityGradientPassed_String                     = Type + " VisibleEnergy_ShowerDensityGradientPassed";
	  const char* VisibleEnergy_ShowerDensityGradientPassed_Name                     = VisibleEnergy_ShowerDensityGradientPassed_String.c_str();
	  fRootHists.VisibleEnergy_ShowerDensityGradientPassed_Hist[Type] = new TH1D(VisibleEnergy_ShowerDensityGradientPassed_Name,VisibleEnergy_ShowerDensityGradientPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_ShowerDensityPowerRemoved_String                     = Type + " VisibleEnergy_ShowerDensityPowerRemoved";
	  const char* VisibleEnergy_ShowerDensityPowerRemoved_Name                     = VisibleEnergy_ShowerDensityPowerRemoved_String.c_str();
	  fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_Hist[Type] = new TH1D(VisibleEnergy_ShowerDensityPowerRemoved_Name,VisibleEnergy_ShowerDensityPowerRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_ShowerDensityPowerPassed_String                     = Type + " VisibleEnergy_ShowerDensityPowerPassed";
	  const char* VisibleEnergy_ShowerDensityPowerPassed_Name                     = VisibleEnergy_ShowerDensityPowerPassed_String.c_str();
	  fRootHists.VisibleEnergy_ShowerDensityPowerPassed_Hist[Type] = new TH1D(VisibleEnergy_ShowerDensityPowerPassed_Name,VisibleEnergy_ShowerDensityPowerPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_ShowerTrackLengthRemoved_String                     = Type + " VisibleEnergy_ShowerTrackLengthRemoved";
	  const char* VisibleEnergy_ShowerTrackLengthRemoved_Name                     = VisibleEnergy_ShowerTrackLengthRemoved_String.c_str();
	  fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_Hist[Type] = new TH1D(VisibleEnergy_ShowerTrackLengthRemoved_Name,VisibleEnergy_ShowerTrackLengthRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_ShowerTrackLengthPassed_String                     = Type + " VisibleEnergy_ShowerTrackLengthPassed";
	  const char* VisibleEnergy_ShowerTrackLengthPassed_Name                     = VisibleEnergy_ShowerTrackLengthPassed_String.c_str();
	  fRootHists.VisibleEnergy_ShowerTrackLengthPassed_Hist[Type] = new TH1D(VisibleEnergy_ShowerTrackLengthPassed_Name,VisibleEnergy_ShowerTrackLengthPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_ShowerTrackWidthRemoved_String                     = Type + " VisibleEnergy_ShowerTrackWidthRemoved";
	  const char* VisibleEnergy_ShowerTrackWidthRemoved_Name                     = VisibleEnergy_ShowerTrackWidthRemoved_String.c_str();
	  fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_Hist[Type] = new TH1D(VisibleEnergy_ShowerTrackWidthRemoved_Name,VisibleEnergy_ShowerTrackWidthRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_ShowerTrackWidthPassed_String                     = Type + " VisibleEnergy_ShowerTrackWidthPassed";
	  const char* VisibleEnergy_ShowerTrackWidthPassed_Name                     = VisibleEnergy_ShowerTrackWidthPassed_String.c_str();
	  fRootHists.VisibleEnergy_ShowerTrackWidthPassed_Hist[Type] = new TH1D(VisibleEnergy_ShowerTrackWidthPassed_Name,VisibleEnergy_ShowerTrackWidthPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	  std::string  VisibleEnergy_MaxTrackLengthRemoved_String                     = Type + " VisibleEnergy_MaxTrackLengthRemoved";
	  const char* VisibleEnergy_MaxTrackLengthRemoved_Name                     = VisibleEnergy_MaxTrackLengthRemoved_String.c_str();
	  fRootHists.VisibleEnergy_MaxTrackLengthRemoved_Hist[Type] = new TH1D(VisibleEnergy_MaxTrackLengthRemoved_Name,VisibleEnergy_MaxTrackLengthRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_MaxTrackLengthPassed_String                     = Type + " VisibleEnergy_MaxTrackLengthPassed";
	  const char* VisibleEnergy_MaxTrackLengthPassed_Name                     = VisibleEnergy_MaxTrackLengthPassed_String.c_str();
	  fRootHists.VisibleEnergy_MaxTrackLengthPassed_Hist[Type] = new TH1D(VisibleEnergy_MaxTrackLengthPassed_Name,VisibleEnergy_MaxTrackLengthPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	  std::string  VisibleEnergy_MaxTrackPIDARemoved_String                     = Type + " VisibleEnergy_MaxTrackPIDARemoved";
	  const char* VisibleEnergy_MaxTrackPIDARemoved_Name                     = VisibleEnergy_MaxTrackPIDARemoved_String.c_str();
	  fRootHists.VisibleEnergy_MaxTrackPIDARemoved_Hist[Type] = new TH1D(VisibleEnergy_MaxTrackPIDARemoved_Name,VisibleEnergy_MaxTrackPIDARemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_MaxTrackPIDAPassed_String                     = Type + " VisibleEnergy_MaxTrackPIDAPassed";
	  const char* VisibleEnergy_MaxTrackPIDAPassed_Name                     = VisibleEnergy_MaxTrackPIDAPassed_String.c_str();
	  fRootHists.VisibleEnergy_MaxTrackPIDAPassed_Hist[Type] = new TH1D(VisibleEnergy_MaxTrackPIDAPassed_Name,VisibleEnergy_MaxTrackPIDAPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_MaxTrackLengthPIDARemoved_String                     = Type + " VisibleEnergy_MaxTrackLengthPIDARemoved";
	  const char* VisibleEnergy_MaxTrackLengthPIDARemoved_Name                     = VisibleEnergy_MaxTrackLengthPIDARemoved_String.c_str();
	  fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_Hist[Type] = new TH1D(VisibleEnergy_MaxTrackLengthPIDARemoved_Name,VisibleEnergy_MaxTrackLengthPIDARemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_MaxTrackLengthPIDAPassed_String                     = Type + " VisibleEnergy_MaxTrackLengthPIDAPassed";
	  const char* VisibleEnergy_MaxTrackLengthPIDAPassed_Name                     = VisibleEnergy_MaxTrackLengthPIDAPassed_String.c_str();
	  fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_Hist[Type] = new TH1D(VisibleEnergy_MaxTrackLengthPIDAPassed_Name,VisibleEnergy_MaxTrackLengthPIDAPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_MVARemoved_String                     = Type + " VisibleEnergy_MVARemoved";
	  const char* VisibleEnergy_MVARemoved_Name                     = VisibleEnergy_MVARemoved_String.c_str();
	  fRootHists.VisibleEnergy_MVARemoved_Hist[Type] = new TH1D(VisibleEnergy_MVARemoved_Name,VisibleEnergy_MVARemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_MVAPassed_String                     = Type + " VisibleEnergy_MVAPassed";
	  const char* VisibleEnergy_MVAPassed_Name                     = VisibleEnergy_MVAPassed_String.c_str();
	  fRootHists.VisibleEnergy_MVAPassed_Hist[Type] = new TH1D(VisibleEnergy_MVAPassed_Name,VisibleEnergy_MVAPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_NeutrinoPdgCodeRemoved_String                     = Type + " VisibleEnergy_NeutrinoPdgCodeRemoved";
	  const char* VisibleEnergy_NeutrinoPdgCodeRemoved_Name                     = VisibleEnergy_NeutrinoPdgCodeRemoved_String.c_str();
	  fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_Hist[Type] = new TH1D(VisibleEnergy_NeutrinoPdgCodeRemoved_Name,VisibleEnergy_NeutrinoPdgCodeRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_NeutrinoPdgCodePassed_String                     = Type + " VisibleEnergy_NeutrinoPdgCodePassed";
	  const char* VisibleEnergy_NeutrinoPdgCodePassed_Name                     = VisibleEnergy_NeutrinoPdgCodePassed_String.c_str();
	  fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_Hist[Type] = new TH1D(VisibleEnergy_NeutrinoPdgCodePassed_Name,VisibleEnergy_NeutrinoPdgCodePassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	  std::string  VisibleEnergy_NumneutrinosRemoved_String                     = Type + " VisibleEnergy_NumneutrinosRemoved";
	  const char* VisibleEnergy_NumneutrinosRemoved_Name                     = VisibleEnergy_NumneutrinosRemoved_String.c_str();
	  fRootHists.VisibleEnergy_NumneutrinosRemoved_Hist[Type] = new TH1D(VisibleEnergy_NumneutrinosRemoved_Name,VisibleEnergy_NumneutrinosRemoved_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	  std::string  VisibleEnergy_NumneutrinosPassed_String                     = Type + " VisibleEnergy_NumneutrinosPassed";
	  const char* VisibleEnergy_NumneutrinosPassed_Name                     = VisibleEnergy_NumneutrinosPassed_String.c_str();
	  fRootHists.VisibleEnergy_NumneutrinosPassed_Hist[Type] = new TH1D(VisibleEnergy_NumneutrinosPassed_Name,VisibleEnergy_NumneutrinosPassed_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);




	  if(fConfig.FillModeHistograms){
	    for(int mode=-1; mode<6; ++mode){
	      std::string  TrueNumber_StringMode                     = TrueNumber_String + " " +  std::to_string(mode);
	      const char* TrueNumber_NameMode                     = TrueNumber_StringMode.c_str();
	      fRootHists.TrueNumber_HistMode[Type][mode] = new TH1D(TrueNumber_NameMode,TrueNumber_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      std::string  VisibleEnergy_BeforeSel_StringMode                     = VisibleEnergy_BeforeSel_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_BeforeSel_NameMode                     = VisibleEnergy_BeforeSel_StringMode.c_str();
	      fRootHists.VisibleEnergy_BeforeSel_HistMode[Type][mode] = new TH1D(VisibleEnergy_BeforeSel_NameMode,VisibleEnergy_BeforeSel_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_AfterSel_StringMode                     = VisibleEnergy_AfterSel_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_AfterSel_NameMode                     = VisibleEnergy_AfterSel_StringMode.c_str();
	      fRootHists.VisibleEnergy_AfterSel_HistMode[Type][mode] = new TH1D(VisibleEnergy_AfterSel_NameMode,VisibleEnergy_AfterSel_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_AfterSelOne_StringMode                     = VisibleEnergy_AfterSelOne_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_AfterSelOne_NameMode                     = VisibleEnergy_AfterSelOne_StringMode.c_str();
	      fRootHists.VisibleEnergy_AfterSelOne_HistMode[Type][mode] = new TH1D(VisibleEnergy_AfterSelOne_NameMode,VisibleEnergy_AfterSelOne_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_AfterSelExtra_StringMode                     = VisibleEnergy_AfterSelExtra_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_AfterSelExtra_NameMode                     = VisibleEnergy_AfterSelExtra_StringMode.c_str();
	      fRootHists.VisibleEnergy_AfterSelExtra_HistMode[Type][mode] = new TH1D(VisibleEnergy_AfterSelExtra_NameMode,VisibleEnergy_AfterSelExtra_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      
	      std::string  VisibleEnergy_MCTruth_StringMode                     = VisibleEnergy_MCTruth_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MCTruth_NameMode                     = VisibleEnergy_MCTruth_StringMode.c_str();
	      fRootHists.VisibleEnergy_MCTruth_HistMode[Type][mode] = new TH1D(VisibleEnergy_MCTruth_NameMode,VisibleEnergy_MCTruth_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);



	      
	      std::string  VisibleEnergy_FVRemoved_StringMode                     = VisibleEnergy_FVRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_FVRemoved_NameMode                     = VisibleEnergy_FVRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_FVRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_FVRemoved_NameMode,VisibleEnergy_FVRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_FVPassed_StringMode                     = VisibleEnergy_FVPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_FVPassed_NameMode                     = VisibleEnergy_FVPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_FVPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_FVPassed_NameMode,VisibleEnergy_FVPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_GtrOneShowerRemoved_StringMode                     = VisibleEnergy_GtrOneShowerRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_GtrOneShowerRemoved_NameMode                     = VisibleEnergy_GtrOneShowerRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_GtrOneShowerRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_GtrOneShowerRemoved_NameMode,VisibleEnergy_GtrOneShowerRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_GtrOneShowerPassed_StringMode                     = VisibleEnergy_GtrOneShowerPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_GtrOneShowerPassed_NameMode                     = VisibleEnergy_GtrOneShowerPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_GtrOneShowerPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_GtrOneShowerPassed_NameMode,VisibleEnergy_GtrOneShowerPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      std::string  VisibleEnergy_OneShowerECutRemoved_StringMode                     = VisibleEnergy_OneShowerECutRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_OneShowerECutRemoved_NameMode                     = VisibleEnergy_OneShowerECutRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_OneShowerECutRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_OneShowerECutRemoved_NameMode,VisibleEnergy_OneShowerECutRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_OneShowerECutPassed_StringMode                     = VisibleEnergy_OneShowerECutPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_OneShowerECutPassed_NameMode                     = VisibleEnergy_OneShowerECutPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_OneShowerECutPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_OneShowerECutPassed_NameMode,VisibleEnergy_OneShowerECutPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_OneShowerResidualRemoved_StringMode                     = VisibleEnergy_OneShowerResidualRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_OneShowerResidualRemoved_NameMode                     = VisibleEnergy_OneShowerResidualRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_OneShowerResidualRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_OneShowerResidualRemoved_NameMode,VisibleEnergy_OneShowerResidualRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_OneShowerResidualPassed_StringMode                     = VisibleEnergy_OneShowerResidualPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_OneShowerResidualPassed_NameMode                     = VisibleEnergy_OneShowerResidualPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_OneShowerResidualPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_OneShowerResidualPassed_NameMode,VisibleEnergy_OneShowerResidualPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      std::string  VisibleEnergy_ConversionGapRemoved_StringMode                     = VisibleEnergy_ConversionGapRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ConversionGapRemoved_NameMode                     = VisibleEnergy_ConversionGapRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_ConversionGapRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_ConversionGapRemoved_NameMode,VisibleEnergy_ConversionGapRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ConversionGapPassed_StringMode                     = VisibleEnergy_ConversionGapPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ConversionGapPassed_NameMode                     = VisibleEnergy_ConversionGapPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_ConversionGapPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_ConversionGapPassed_NameMode,VisibleEnergy_ConversionGapPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_dEdxRemoved_StringMode                     = VisibleEnergy_dEdxRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_dEdxRemoved_NameMode                     = VisibleEnergy_dEdxRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_dEdxRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_dEdxRemoved_NameMode,VisibleEnergy_dEdxRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_dEdxPassed_StringMode                     = VisibleEnergy_dEdxPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_dEdxPassed_NameMode                     = VisibleEnergy_dEdxPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_dEdxPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_dEdxPassed_NameMode,VisibleEnergy_dEdxPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_LengthRemoved_StringMode                     = VisibleEnergy_LengthRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_LengthRemoved_NameMode                     = VisibleEnergy_LengthRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_LengthRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_LengthRemoved_NameMode,VisibleEnergy_LengthRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_LengthPassed_StringMode                     = VisibleEnergy_LengthPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_LengthPassed_NameMode                     = VisibleEnergy_LengthPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_LengthPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_LengthPassed_NameMode,VisibleEnergy_LengthPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_OpeningAngleRemoved_StringMode                     = VisibleEnergy_OpeningAngleRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_OpeningAngleRemoved_NameMode                     = VisibleEnergy_OpeningAngleRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_OpeningAngleRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_OpeningAngleRemoved_NameMode,VisibleEnergy_OpeningAngleRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_OpeningAnglePassed_StringMode                     = VisibleEnergy_OpeningAnglePassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_OpeningAnglePassed_NameMode                     = VisibleEnergy_OpeningAnglePassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_OpeningAnglePassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_OpeningAnglePassed_NameMode,VisibleEnergy_OpeningAnglePassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerEnergyRemoved_StringMode                     = VisibleEnergy_ShowerEnergyRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerEnergyRemoved_NameMode                     = VisibleEnergy_ShowerEnergyRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerEnergyRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerEnergyRemoved_NameMode,VisibleEnergy_ShowerEnergyRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerEnergyPassed_StringMode                     = VisibleEnergy_ShowerEnergyPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerEnergyPassed_NameMode                     = VisibleEnergy_ShowerEnergyPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerEnergyPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerEnergyPassed_NameMode,VisibleEnergy_ShowerEnergyPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerDensityGradientRemoved_StringMode                     = VisibleEnergy_ShowerDensityGradientRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerDensityGradientRemoved_NameMode                     = VisibleEnergy_ShowerDensityGradientRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerDensityGradientRemoved_NameMode,VisibleEnergy_ShowerDensityGradientRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerDensityGradientPassed_StringMode                     = VisibleEnergy_ShowerDensityGradientPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerDensityGradientPassed_NameMode                     = VisibleEnergy_ShowerDensityGradientPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerDensityGradientPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerDensityGradientPassed_NameMode,VisibleEnergy_ShowerDensityGradientPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerDensityPowerRemoved_StringMode                     = VisibleEnergy_ShowerDensityPowerRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerDensityPowerRemoved_NameMode                     = VisibleEnergy_ShowerDensityPowerRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerDensityPowerRemoved_NameMode,VisibleEnergy_ShowerDensityPowerRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerDensityPowerPassed_StringMode                     = VisibleEnergy_ShowerDensityPowerPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerDensityPowerPassed_NameMode                     = VisibleEnergy_ShowerDensityPowerPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerDensityPowerPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerDensityPowerPassed_NameMode,VisibleEnergy_ShowerDensityPowerPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerTrackLengthRemoved_StringMode                     = VisibleEnergy_ShowerTrackLengthRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerTrackLengthRemoved_NameMode                     = VisibleEnergy_ShowerTrackLengthRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerTrackLengthRemoved_NameMode,VisibleEnergy_ShowerTrackLengthRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerTrackLengthPassed_StringMode                     = VisibleEnergy_ShowerTrackLengthPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerTrackLengthPassed_NameMode                     = VisibleEnergy_ShowerTrackLengthPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerTrackLengthPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerTrackLengthPassed_NameMode,VisibleEnergy_ShowerTrackLengthPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerTrackWidthRemoved_StringMode                     = VisibleEnergy_ShowerTrackWidthRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerTrackWidthRemoved_NameMode                     = VisibleEnergy_ShowerTrackWidthRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerTrackWidthRemoved_NameMode,VisibleEnergy_ShowerTrackWidthRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerTrackWidthPassed_StringMode                     = VisibleEnergy_ShowerTrackWidthPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_ShowerTrackWidthPassed_NameMode                     = VisibleEnergy_ShowerTrackWidthPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_ShowerTrackWidthPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_ShowerTrackWidthPassed_NameMode,VisibleEnergy_ShowerTrackWidthPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_MaxTrackLengthRemoved_StringMode                     = VisibleEnergy_MaxTrackLengthRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MaxTrackLengthRemoved_NameMode                     = VisibleEnergy_MaxTrackLengthRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_MaxTrackLengthRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_MaxTrackLengthRemoved_NameMode,VisibleEnergy_MaxTrackLengthRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_MaxTrackLengthPassed_StringMode                     = VisibleEnergy_MaxTrackLengthPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MaxTrackLengthPassed_NameMode                     = VisibleEnergy_MaxTrackLengthPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_MaxTrackLengthPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_MaxTrackLengthPassed_NameMode,VisibleEnergy_MaxTrackLengthPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_MaxTrackPIDARemoved_StringMode                     = VisibleEnergy_MaxTrackPIDARemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MaxTrackPIDARemoved_NameMode                     = VisibleEnergy_MaxTrackPIDARemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_MaxTrackPIDARemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_MaxTrackPIDARemoved_NameMode,VisibleEnergy_MaxTrackPIDARemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_MaxTrackPIDAPassed_StringMode                     = VisibleEnergy_MaxTrackPIDAPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MaxTrackPIDAPassed_NameMode                     = VisibleEnergy_MaxTrackPIDAPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_MaxTrackPIDAPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_MaxTrackPIDAPassed_NameMode,VisibleEnergy_MaxTrackPIDAPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_MaxTrackLengthPIDARemoved_StringMode                     = VisibleEnergy_MaxTrackLengthPIDARemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MaxTrackLengthPIDARemoved_NameMode                     = VisibleEnergy_MaxTrackLengthPIDARemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_MaxTrackLengthPIDARemoved_NameMode,VisibleEnergy_MaxTrackLengthPIDARemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_MaxTrackLengthPIDAPassed_StringMode                     = VisibleEnergy_MaxTrackLengthPIDAPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MaxTrackLengthPIDAPassed_NameMode                     = VisibleEnergy_MaxTrackLengthPIDAPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_MaxTrackLengthPIDAPassed_NameMode,VisibleEnergy_MaxTrackLengthPIDAPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_MVARemoved_StringMode                     = VisibleEnergy_MVARemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MVARemoved_NameMode                     = VisibleEnergy_MVARemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_MVARemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_MVARemoved_NameMode,VisibleEnergy_MVARemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_MVAPassed_StringMode                     = VisibleEnergy_MVAPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_MVAPassed_NameMode                     = VisibleEnergy_MVAPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_MVAPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_MVAPassed_NameMode,VisibleEnergy_MVAPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_NeutrinoPdgCodeRemoved_StringMode                     = VisibleEnergy_NeutrinoPdgCodeRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_NeutrinoPdgCodeRemoved_NameMode                     = VisibleEnergy_NeutrinoPdgCodeRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_NeutrinoPdgCodeRemoved_NameMode,VisibleEnergy_NeutrinoPdgCodeRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_NeutrinoPdgCodePassed_StringMode                     = VisibleEnergy_NeutrinoPdgCodePassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_NeutrinoPdgCodePassed_NameMode                     = VisibleEnergy_NeutrinoPdgCodePassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_NeutrinoPdgCodePassed_NameMode,VisibleEnergy_NeutrinoPdgCodePassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_NumneutrinosRemoved_StringMode                     = VisibleEnergy_NumneutrinosRemoved_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_NumneutrinosRemoved_NameMode                     = VisibleEnergy_NumneutrinosRemoved_StringMode.c_str();
	      fRootHists.VisibleEnergy_NumneutrinosRemoved_HistMode[Type][mode] = new TH1D(VisibleEnergy_NumneutrinosRemoved_NameMode,VisibleEnergy_NumneutrinosRemoved_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_NumneutrinosPassed_StringMode                     = VisibleEnergy_NumneutrinosPassed_String + " " +  std::to_string(mode);
	      const char* VisibleEnergy_NumneutrinosPassed_NameMode                     = VisibleEnergy_NumneutrinosPassed_StringMode.c_str();
	      fRootHists.VisibleEnergy_NumneutrinosPassed_HistMode[Type][mode] = new TH1D(VisibleEnergy_NumneutrinosPassed_NameMode,VisibleEnergy_NumneutrinosPassed_NameMode,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	    }
	  }

	  if(fConfig.FillIntTypeHistograms){
	    for(int inttype=1000; inttype<1101; ++inttype){
	      std::string  TrueNumber_StringIntType                     = TrueNumber_String + " " +  std::to_string(inttype);
	      const char* TrueNumber_NameIntType                     = TrueNumber_StringIntType.c_str();
	      fRootHists.TrueNumber_HistIntType[Type][inttype] = new TH1D(TrueNumber_NameIntType,TrueNumber_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_BeforeSel_StringIntType                     = VisibleEnergy_BeforeSel_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_BeforeSel_NameIntType                     = VisibleEnergy_BeforeSel_StringIntType.c_str();
	      fRootHists.VisibleEnergy_BeforeSel_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_BeforeSel_NameIntType,VisibleEnergy_BeforeSel_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_AfterSel_StringIntType                     = VisibleEnergy_AfterSel_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_AfterSel_NameIntType                     = VisibleEnergy_AfterSel_StringIntType.c_str();
	      fRootHists.VisibleEnergy_AfterSel_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_AfterSel_NameIntType,VisibleEnergy_AfterSel_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_AfterSelOne_StringIntType                     = VisibleEnergy_AfterSelOne_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_AfterSelOne_NameIntType                     = VisibleEnergy_AfterSelOne_StringIntType.c_str();
	      fRootHists.VisibleEnergy_AfterSelOne_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_AfterSelOne_NameIntType,VisibleEnergy_AfterSelOne_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_AfterSelExtra_StringIntType                     = VisibleEnergy_AfterSelExtra_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_AfterSelExtra_NameIntType                     = VisibleEnergy_AfterSelExtra_StringIntType.c_str();
	      fRootHists.VisibleEnergy_AfterSelExtra_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_AfterSelExtra_NameIntType,VisibleEnergy_AfterSelExtra_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      std::string  VisibleEnergy_MCTruth_StringIntType                     = VisibleEnergy_MCTruth_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MCTruth_NameIntType                     = VisibleEnergy_MCTruth_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MCTruth_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MCTruth_NameIntType,VisibleEnergy_MCTruth_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);



	      	      std::string  VisibleEnergy_FVRemoved_StringIntType                     = VisibleEnergy_FVRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_FVRemoved_NameIntType                     = VisibleEnergy_FVRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_FVRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_FVRemoved_NameIntType,VisibleEnergy_FVRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_FVPassed_StringIntType                     = VisibleEnergy_FVPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_FVPassed_NameIntType                     = VisibleEnergy_FVPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_FVPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_FVPassed_NameIntType,VisibleEnergy_FVPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_GtrOneShowerRemoved_StringIntType                     = VisibleEnergy_GtrOneShowerRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_GtrOneShowerRemoved_NameIntType                     = VisibleEnergy_GtrOneShowerRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_GtrOneShowerRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_GtrOneShowerRemoved_NameIntType,VisibleEnergy_GtrOneShowerRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_GtrOneShowerPassed_StringIntType                     = VisibleEnergy_GtrOneShowerPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_GtrOneShowerPassed_NameIntType                     = VisibleEnergy_GtrOneShowerPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_GtrOneShowerPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_GtrOneShowerPassed_NameIntType,VisibleEnergy_GtrOneShowerPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      std::string  VisibleEnergy_OneShowerECutRemoved_StringIntType                     = VisibleEnergy_OneShowerECutRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_OneShowerECutRemoved_NameIntType                     = VisibleEnergy_OneShowerECutRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_OneShowerECutRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_OneShowerECutRemoved_NameIntType,VisibleEnergy_OneShowerECutRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_OneShowerECutPassed_StringIntType                     = VisibleEnergy_OneShowerECutPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_OneShowerECutPassed_NameIntType                     = VisibleEnergy_OneShowerECutPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_OneShowerECutPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_OneShowerECutPassed_NameIntType,VisibleEnergy_OneShowerECutPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_OneShowerResidualRemoved_StringIntType                     = VisibleEnergy_OneShowerResidualRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_OneShowerResidualRemoved_NameIntType                     = VisibleEnergy_OneShowerResidualRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_OneShowerResidualRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_OneShowerResidualRemoved_NameIntType,VisibleEnergy_OneShowerResidualRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_OneShowerResidualPassed_StringIntType                     = VisibleEnergy_OneShowerResidualPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_OneShowerResidualPassed_NameIntType                     = VisibleEnergy_OneShowerResidualPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_OneShowerResidualPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_OneShowerResidualPassed_NameIntType,VisibleEnergy_OneShowerResidualPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	      std::string  VisibleEnergy_ConversionGapRemoved_StringIntType                     = VisibleEnergy_ConversionGapRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ConversionGapRemoved_NameIntType                     = VisibleEnergy_ConversionGapRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ConversionGapRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ConversionGapRemoved_NameIntType,VisibleEnergy_ConversionGapRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ConversionGapPassed_StringIntType                     = VisibleEnergy_ConversionGapPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ConversionGapPassed_NameIntType                     = VisibleEnergy_ConversionGapPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ConversionGapPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ConversionGapPassed_NameIntType,VisibleEnergy_ConversionGapPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_dEdxRemoved_StringIntType                     = VisibleEnergy_dEdxRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_dEdxRemoved_NameIntType                     = VisibleEnergy_dEdxRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_dEdxRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_dEdxRemoved_NameIntType,VisibleEnergy_dEdxRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_dEdxPassed_StringIntType                     = VisibleEnergy_dEdxPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_dEdxPassed_NameIntType                     = VisibleEnergy_dEdxPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_dEdxPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_dEdxPassed_NameIntType,VisibleEnergy_dEdxPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_LengthRemoved_StringIntType                     = VisibleEnergy_LengthRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_LengthRemoved_NameIntType                     = VisibleEnergy_LengthRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_LengthRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_LengthRemoved_NameIntType,VisibleEnergy_LengthRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_LengthPassed_StringIntType                     = VisibleEnergy_LengthPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_LengthPassed_NameIntType                     = VisibleEnergy_LengthPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_LengthPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_LengthPassed_NameIntType,VisibleEnergy_LengthPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_OpeningAngleRemoved_StringIntType                     = VisibleEnergy_OpeningAngleRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_OpeningAngleRemoved_NameIntType                     = VisibleEnergy_OpeningAngleRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_OpeningAngleRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_OpeningAngleRemoved_NameIntType,VisibleEnergy_OpeningAngleRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_OpeningAnglePassed_StringIntType                     = VisibleEnergy_OpeningAnglePassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_OpeningAnglePassed_NameIntType                     = VisibleEnergy_OpeningAnglePassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_OpeningAnglePassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_OpeningAnglePassed_NameIntType,VisibleEnergy_OpeningAnglePassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerEnergyRemoved_StringIntType                     = VisibleEnergy_ShowerEnergyRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerEnergyRemoved_NameIntType                     = VisibleEnergy_ShowerEnergyRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerEnergyRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerEnergyRemoved_NameIntType,VisibleEnergy_ShowerEnergyRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerEnergyPassed_StringIntType                     = VisibleEnergy_ShowerEnergyPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerEnergyPassed_NameIntType                     = VisibleEnergy_ShowerEnergyPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerEnergyPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerEnergyPassed_NameIntType,VisibleEnergy_ShowerEnergyPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerDensityGradientRemoved_StringIntType                     = VisibleEnergy_ShowerDensityGradientRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerDensityGradientRemoved_NameIntType                     = VisibleEnergy_ShowerDensityGradientRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerDensityGradientRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerDensityGradientRemoved_NameIntType,VisibleEnergy_ShowerDensityGradientRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerDensityGradientPassed_StringIntType                     = VisibleEnergy_ShowerDensityGradientPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerDensityGradientPassed_NameIntType                     = VisibleEnergy_ShowerDensityGradientPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerDensityGradientPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerDensityGradientPassed_NameIntType,VisibleEnergy_ShowerDensityGradientPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerDensityPowerRemoved_StringIntType                     = VisibleEnergy_ShowerDensityPowerRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerDensityPowerRemoved_NameIntType                     = VisibleEnergy_ShowerDensityPowerRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerDensityPowerRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerDensityPowerRemoved_NameIntType,VisibleEnergy_ShowerDensityPowerRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerDensityPowerPassed_StringIntType                     = VisibleEnergy_ShowerDensityPowerPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerDensityPowerPassed_NameIntType                     = VisibleEnergy_ShowerDensityPowerPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerDensityPowerPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerDensityPowerPassed_NameIntType,VisibleEnergy_ShowerDensityPowerPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerTrackLengthRemoved_StringIntType                     = VisibleEnergy_ShowerTrackLengthRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerTrackLengthRemoved_NameIntType                     = VisibleEnergy_ShowerTrackLengthRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerTrackLengthRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerTrackLengthRemoved_NameIntType,VisibleEnergy_ShowerTrackLengthRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerTrackLengthPassed_StringIntType                     = VisibleEnergy_ShowerTrackLengthPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerTrackLengthPassed_NameIntType                     = VisibleEnergy_ShowerTrackLengthPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerTrackLengthPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerTrackLengthPassed_NameIntType,VisibleEnergy_ShowerTrackLengthPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_ShowerTrackWidthRemoved_StringIntType                     = VisibleEnergy_ShowerTrackWidthRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerTrackWidthRemoved_NameIntType                     = VisibleEnergy_ShowerTrackWidthRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerTrackWidthRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerTrackWidthRemoved_NameIntType,VisibleEnergy_ShowerTrackWidthRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_ShowerTrackWidthPassed_StringIntType                     = VisibleEnergy_ShowerTrackWidthPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_ShowerTrackWidthPassed_NameIntType                     = VisibleEnergy_ShowerTrackWidthPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_ShowerTrackWidthPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_ShowerTrackWidthPassed_NameIntType,VisibleEnergy_ShowerTrackWidthPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_MaxTrackLengthRemoved_StringIntType                     = VisibleEnergy_MaxTrackLengthRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MaxTrackLengthRemoved_NameIntType                     = VisibleEnergy_MaxTrackLengthRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MaxTrackLengthRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MaxTrackLengthRemoved_NameIntType,VisibleEnergy_MaxTrackLengthRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_MaxTrackLengthPassed_StringIntType                     = VisibleEnergy_MaxTrackLengthPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MaxTrackLengthPassed_NameIntType                     = VisibleEnergy_MaxTrackLengthPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MaxTrackLengthPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MaxTrackLengthPassed_NameIntType,VisibleEnergy_MaxTrackLengthPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_MaxTrackPIDARemoved_StringIntType                     = VisibleEnergy_MaxTrackPIDARemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MaxTrackPIDARemoved_NameIntType                     = VisibleEnergy_MaxTrackPIDARemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MaxTrackPIDARemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MaxTrackPIDARemoved_NameIntType,VisibleEnergy_MaxTrackPIDARemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_MaxTrackPIDAPassed_StringIntType                     = VisibleEnergy_MaxTrackPIDAPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MaxTrackPIDAPassed_NameIntType                     = VisibleEnergy_MaxTrackPIDAPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MaxTrackPIDAPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MaxTrackPIDAPassed_NameIntType,VisibleEnergy_MaxTrackPIDAPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_MaxTrackLengthPIDARemoved_StringIntType                     = VisibleEnergy_MaxTrackLengthPIDARemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MaxTrackLengthPIDARemoved_NameIntType                     = VisibleEnergy_MaxTrackLengthPIDARemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MaxTrackLengthPIDARemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MaxTrackLengthPIDARemoved_NameIntType,VisibleEnergy_MaxTrackLengthPIDARemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_MaxTrackLengthPIDAPassed_StringIntType                     = VisibleEnergy_MaxTrackLengthPIDAPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MaxTrackLengthPIDAPassed_NameIntType                     = VisibleEnergy_MaxTrackLengthPIDAPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MaxTrackLengthPIDAPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MaxTrackLengthPIDAPassed_NameIntType,VisibleEnergy_MaxTrackLengthPIDAPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_MVARemoved_StringIntType                     = VisibleEnergy_MVARemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MVARemoved_NameIntType                     = VisibleEnergy_MVARemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MVARemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MVARemoved_NameIntType,VisibleEnergy_MVARemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_MVAPassed_StringIntType                     = VisibleEnergy_MVAPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_MVAPassed_NameIntType                     = VisibleEnergy_MVAPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_MVAPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_MVAPassed_NameIntType,VisibleEnergy_MVAPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_NeutrinoPdgCodeRemoved_StringIntType                     = VisibleEnergy_NeutrinoPdgCodeRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_NeutrinoPdgCodeRemoved_NameIntType                     = VisibleEnergy_NeutrinoPdgCodeRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_NeutrinoPdgCodeRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NeutrinoPdgCodeRemoved_NameIntType,VisibleEnergy_NeutrinoPdgCodeRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_NeutrinoPdgCodePassed_StringIntType                     = VisibleEnergy_NeutrinoPdgCodePassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_NeutrinoPdgCodePassed_NameIntType                     = VisibleEnergy_NeutrinoPdgCodePassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_NeutrinoPdgCodePassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NeutrinoPdgCodePassed_NameIntType,VisibleEnergy_NeutrinoPdgCodePassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

	      std::string  VisibleEnergy_NumneutrinosRemoved_StringIntType                     = VisibleEnergy_NumneutrinosRemoved_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_NumneutrinosRemoved_NameIntType                     = VisibleEnergy_NumneutrinosRemoved_StringIntType.c_str();
	      fRootHists.VisibleEnergy_NumneutrinosRemoved_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NumneutrinosRemoved_NameIntType,VisibleEnergy_NumneutrinosRemoved_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	      std::string  VisibleEnergy_NumneutrinosPassed_StringIntType                     = VisibleEnergy_NumneutrinosPassed_String + " " +  std::to_string(inttype);
	      const char* VisibleEnergy_NumneutrinosPassed_NameIntType                     = VisibleEnergy_NumneutrinosPassed_StringIntType.c_str();
	      fRootHists.VisibleEnergy_NumneutrinosPassed_HistIntType[Type][inttype] = new TH1D(VisibleEnergy_NumneutrinosPassed_NameIntType,VisibleEnergy_NumneutrinosPassed_NameIntType,fRootHists.ebins,fRootHists.emin,fRootHists.emax);



	    }
	  }
	}
      }
    }
  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelectionReco)

