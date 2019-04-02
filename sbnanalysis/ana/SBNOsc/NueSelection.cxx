//C++ Includes 
#include <iostream>
#include <vector>

//Root Includes 
#include <TH2D.h>

//Framework Includes 
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"

//Larsoft Includes 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCFlux.h"

//SBN Includes 
#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

    NueSelection::NueSelection() : SelectionBase(), EventCounter(0), NuCount(0) {}

    
    void NueSelection::Initialize(fhicl::ParameterSet* config) {

      hello();

      fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("NueSelection");
      
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
      fConfig.minLengthExitingTrack           = pconfig.get<double>("minLengthExitingTrack",100); 
      fConfig.trackVisibleEnergyThreshold     = pconfig.get<double>("trackVisibleEnergyThreshold",0.1);
      fConfig.showerVisibleEnergyThreshold    = pconfig.get<double>("showerVisibleEnergyThreshold",200); 
      fConfig.showerEnergyDistortion          = pconfig.get<double>("showerEnergyDistortion",1); 
      fConfig.trackEnergyDistortion           = pconfig.get<double>("trackEnergyDistortion",1); 
      fConfig.photonVisibleEnergyThreshold    = pconfig.get<double>("photonVisibleEnergyThreshold",100);
      fConfig.leptonEnergyDistortionContained = pconfig.get<double>("leptonEnergyDistortionContained",1); 
      fConfig.nueEfficency                    = pconfig.get<double>("nueEfficencyWeight",0.8);
      fConfig.vtxEnergyCut                    = pconfig.get<double>("vtxEnergyCut",20); 
      fConfig.photonConvLenghCut              = pconfig.get<double>("photonConvLenghCut",2); 
      fConfig.dEdxPhotonCut                   = pconfig.get<double>("dEdxPhotonCut",0.06);
      fConfig.GlobalWeight                    = pconfig.get<double>("GlobalWeight",1);
      fConfig.POTWeight                       = pconfig.get<double>("POTWeight",1);
      fConfig.CosmicWeight                    = pconfig.get<double>("CosmicWeight",0.05);
      fConfig.CosmicVolumeRadius              = pconfig.get<double>("CosmicVolumeRadius",15);
      fConfig.UniformWeights = pconfig.get<std::vector<std::string>>("UniformWeights", {});

      fConfig.Detector               = pconfig.get<std::string>("Detector"," "); 
      
      fConfig.ApplyNueEfficiency     = pconfig.get<bool>("ApplyNueEfficiency",false);
      fConfig.ApplyElectronEnergyCut = pconfig.get<bool>("ApplyElectronEnergyCut",false);
      fConfig.ApplyFVCut             = pconfig.get<bool>("ApplyFVCut",false);
      fConfig.ApplyAVCut             = pconfig.get<bool>("ApplyAVCut",false);
      fConfig.ApplyPhotonEnergyCut   = pconfig.get<bool>("ApplyPhotonEnergyCut",false); 
      fConfig.ApplyConversionGapCut  = pconfig.get<bool>("ApplyConversionGapCut",false);
      fConfig.ApplydEdxCut           = pconfig.get<bool>("ApplydEdxCut",false);
      fConfig.ApplyMuonLenghtCut     = pconfig.get<bool>("ApplyMuonLenghtCut",false);
      fConfig.ApplyKMECCut           = pconfig.get<bool>("ApplyKMECCut",false);
      fConfig.ApplyCosmicCylinderCut = pconfig.get<bool>("ApplyCosmicCylinderCut",false);
      fConfig.ApplyCosmicFVCut       = pconfig.get<bool>("ApplyCosmicFVCut",false);

      fConfig.IncludeCosmics         = pconfig.get<bool>("IncludeCosmics",true);
      fConfig.IncludeDirt            = pconfig.get<bool>("IncludeDirt",true);  
      fConfig.CosmicsOnly            = pconfig.get<bool>("Cosmics Only",false);   
      fConfig.DirtOnly               = pconfig.get<bool>("DirtOnly",false); 
      fConfig.Verbose                = pconfig.get<bool>("Verbose",false);  
      
      //Set up the selection histograms 
      fOutputFile->cd();

      //Initilise the histograms 
      InitialiseHistograms();

      //Time the CPU.
      c_start = std::clock();

      MCTruthCounter =0;
    }
    
    
    void NueSelection::Finalize(){

      //Read the End time
      c_end = std::clock(); 
      long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
      std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

      std::cout << "Neutrinos selected: " <<  NuCount << " in " << EventCounter << " readout events" << std::endl;

      std::cout << "MCTruth Events count was: " <<  MCTruthCounter << std::endl;

      //Set up the selection histograms                                                  
      fOutputFile->cd();    
      gDirectory->mkdir("Histograms");
      fOutputFile->cd("Histograms");

      fEventHist->Write();

      fRootHists.VisibleEnergy_FNKPDecay_Hist->Write();
      fRootHists.VisibleEnergy_FNKMDecay_Hist->Write();
      fRootHists.VisibleEnergy_MuDecays_Hist->Write();

      for(int i=0; i<fRootHists.HistTypes.size(); ++i){
	fRootHists.TrueEnergy_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.TrueEnergyAll_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.CCQEEnergy_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_AVCut_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_FVCut_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_EnergyCut_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_PhotonEnergyCut_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_ConversionGapCut_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_MuLenghtCut_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_NCCut_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_LeptonPlusPhotonCut_Hist[fRootHists.HistTypes[i]]->Write();
	fRootHists.VisibleEnergy_Selection_Hist[fRootHists.HistTypes[i]]->Write();
      }

      fRootHists.VisibleEnergy_CosmicFVCut_Hist->Write();
      fRootHists.VisibleEnergy_CosmicClyinderCut_Hist->Write();
      fRootHists.VisibleEnergy_CosmicdEdxCut_Hist->Write();
      fRootHists.VisibleEnergy_CosmicWeightCut_Hist->Write();
      fRootHists.VisibleEnergy_CosmicEnergyCut_Hist->Write();
      fRootHists.VisibleEnergy_CosmicSelection_Hist->Write();
    }

    bool NueSelection::ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco){
      if (EventCounter % 10 == 0) {
        std::cout << "NueSelection: Processing event " << EventCounter << " "
                  << "(" << NuCount << " neutrinos selected)"
                  << std::endl;
      }

      std::cout << "New Event" << std::endl;

      bool selected = false;

      EventCounter++;

      //Lets keep track 
      fEventHist->Fill(1);

      //Grab a neutrino datat product from the event
      auto const& mctruths =  *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
      
      //Get tracks and showers 
      auto const& mctracks  = *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
      auto const& mcshowers = *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

      //Get the flux info 
      auto const& mcflux = ev.getValidHandle<std::vector<simb::MCFlux> >(fTruthTag);
      
      //Get a map of the all the particles (not just the final state ones).
      std::map<int,const simb::MCParticle*> mcparticles;
      auto const& mcparticle_list = *ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);
      for(auto const &mcparticle:  mcparticle_list){
	
	mcparticles[mcparticle.TrackId()] = &mcparticle;

		if(fConfig.Verbose){
		  std::cout << "MC Particle with track ID: " << mcparticle.TrackId() << " has pdgcode: " << mcparticle.PdgCode() << " with energy: " << mcparticle.E() << std::endl;
		}
      }
      
      float totalbnbweight = 1;

      //Remove for debugging perposes
      //      for (auto const &mct: mctracks) {
      // 	double mass = PDGMass(mct.PdgCode());
      // 	std::cout << "Track with id: " << mct.TrackID() << " has pdgcode: " << mct.PdgCode() << " and energy: " << (mct.Start().E() - mass) / 1000 << std::endl;
      // }
      // for (auto const &mcs: mcshowers) {
      // 	double mass = PDGMass(mcs.PdgCode());
      // 	std::cout << "Shower with id: " << mcs.TrackID() << " has pdgcode: " << mcs.PdgCode() << " and energy: " << (mcs.Start().E() - mass) / 1000 << std::endl;
      // }

    
      //Iterate through the neutrinos
      for (size_t i=0; i<mctruths.size(); i++) {

	//lazy pea ranch
	if(fConfig.CosmicsOnly){continue;}  

	auto const& mctruth = mctruths.at(i);  

	++MCTruthCounter;

	if(fConfig.Verbose){
	  std::cout << "##########################################################################################################" << std::endl;
	  std::cout << mctruth << std::endl;
	}
	

        const simb::MCNeutrino& nu = mctruth.GetNeutrino();
      
	//KMec didn't exist in the past might not be correct but remove for proposal
	bool pass_kMEC = !(fConfig.ApplyKMECCut && nu.Mode() == simb::kMEC); 
	if(!pass_kMEC){
	  if(fConfig.Verbose){std::cout << "Event Was kMEC. Event was not selected" << std::endl;}
	  return false;
	}
      
	//holder for the electron track ID
	int leptontrackID = -1;
	float lepton_energy = -99999;
	
	//Calculate the lepton track ID
	for(auto const &mcparticle_data: mcparticles){
	  const simb::MCParticle*  mcparticle = mcparticle_data.second;
	  if(TMath::Abs(mcparticle->PdgCode()) == 11){
	    //Check to see if its from the vertex.
	    if(isFromNuVertex(mctruth,mcparticle)){
	      //Say the lepton is the one with the highest energy.
	      if(mcparticle->E() > lepton_energy){
		lepton_energy = mcparticle->E();
		leptontrackID = mcparticle->TrackId();
	      }
	    }
	  }
	}

	if(fConfig.Verbose){
	  std::cout << "leptontrackID: " << leptontrackID << " lepton_energy: " << lepton_energy << std::endl;
	}
	    
	//Check if the neutrino is contained.
	bool containedinfv = containedInFV(nu.Nu().Position().Vect());

	//Set the values required to calculate the track and shower energys + distortion 
	VisibleEnergyCalculator calculator;
	calculator.lepton_pdgid = 11;
	calculator.track_threshold =  fConfig.trackVisibleEnergyThreshold;
	calculator.shower_energy_distortion = fConfig.showerEnergyDistortion;
	calculator.track_energy_distortion = fConfig.trackEnergyDistortion;
	calculator.lepton_energy_distortion_contained = fConfig.leptonEnergyDistortionContained;
	calculator.lepton_energy_distortion_leaving_A = 1; //Not required here.
	calculator.lepton_energy_distortion_leaving_B = 1; //Not required here. 
	calculator.lepton_contained = containedinfv; // We don't worrry at the moment if the electron is contained at this moment.
	calculator.lepton_index = leptontrackID;

	if(fConfig.Verbose){
	  if(mctruth.GetNeutrino().CCNC() == simb::kCC){
	    std::cout << "is charged current" << std::endl;
	  }
	  else if(mctruth.GetNeutrino().CCNC() == simb::kNC) {
	    std::cout << "is neutral current" << std::endl; 
	  }
	  else{
	    std::cout << "Somehow other. Maybe cosmic if you have some cosmics" << std::endl;
	  }
	}
	    

	// build the interaction
        Event::Interaction interaction = truth[i];

	//Calculate the Energy 
	const std::vector<double> visible_energy = FlavourEnergyDeposition(rand, mctruth, mcparticles,calculator);
	double Hadronic_energy = visible_energy[0];
	double Shower_energy   = visible_energy[1];
	double Leptonic_energy = visible_energy[2];
	//	double Nue_energy = Hadronic_energy + Shower_energy + Leptonic_energy;

	//Weighting for efficiency of selection.
	double weight = 1;

	///Weight with the POT
	weight *= fConfig.POTWeight;

	//Add the Global weightings if any.
	weight *= fConfig.GlobalWeight;

	//Apply Uniform weight for bnb correction 
	if(fConfig.Detector != "SBND"){
	  for (auto const &key: fConfig.UniformWeights) {
	    weight *= interaction.weights.at(key)[0];
	    std::cout << "weight: " << weight << std::endl;
	    totalbnbweight *=  interaction.weights.at(key)[0];     
	  }
	}
	
	//Print out the decay type. Find enums here: http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/
	int initnu =  mcflux->at(i).fntype;
	int fnd = mcflux->at(i).fndecay;

 	//Make a nue interaction 
	  NueSelection::NueInteraction intInfo({Hadronic_energy, Shower_energy, Leptonic_energy, weight,initnu,leptontrackID,11,fnd}); 

	//Selection Criteria 
	//Check for Tracks - Remove any track that come from the neutrino if 1m 
        //Find electron vertex candiates - Remove 80% due to fake dEdx cut 
        //If pi0 is made 2 photons come out as candiates in both are above 100 MeV remove.
        //If hadrons produce have more than 50 MeV this is visible if all photons pair produce more than 3cm away from the vertex remove.
        //Remove 94% pion of events as dEdx cut
        //Event could be CC0Pi Muons Interaction with photon background which mimics the CC1Pi Electron interatactions. these need to go throught the cuts above as well.  
        //Run the Selection

	//Selection. 
	bool selection = Select(ev, mctruth, mctracks, i, mcparticles, intInfo);

	//Make the reco interaction event. 
	Event::RecoInteraction reco_interaction(interaction, i);
        reco_interaction.reco_energy = intInfo.GetNueEnergy();
	reco_interaction.weight = intInfo.weight; 
	
	if(selection){

	  //Make the reco interaction event. 
	  Event::RecoInteraction reco_interaction(interaction, i);
	  reco_interaction.reco_energy = intInfo.GetNueEnergy();
	  reco_interaction.weight = intInfo.weight; 
          reco.push_back(reco_interaction);
	  selected =  selection;
	  
	  //Fill the histograms
	  FillHistograms(fRootHists.VisibleEnergy_Selection_Hist,nu,intInfo,dirtevent);
	  FillHistograms(fRootHists.TrueEnergy_Hist,nu,intInfo,nu.Nu().E(),dirtevent);
	  NuCount++;

	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    if(intInfo.fnd == 5 || intInfo.fnd == 6 || intInfo.fnd == 7){fRootHists.VisibleEnergy_FNKPDecay_Hist->Fill(intInfo.GetNueEnergy());}
	    if(intInfo.fnd == 8 || intInfo.fnd == 9 || intInfo.fnd == 10){fRootHists.VisibleEnergy_FNKMDecay_Hist->Fill(intInfo.GetNueEnergy());}  
	    if(intInfo.fnd == 11 || intInfo.fnd == 12){fRootHists.VisibleEnergy_MuDecays_Hist->Fill(intInfo.GetNueEnergy());} 
	  }
        }
	
	if(fConfig.Verbose){
	  
	  if(selection){
	    std::cout << "Event Was Selected with energy:" << intInfo.GetNueEnergy()<<"done" << std::endl;
	    if(mctruth.GetNeutrino().CCNC() == simb::kNC){std::cout << "Was NC event: " << nu.Nu().PdgCode() << std::endl;}
	  }

	  PrintInformation(mctruth, intInfo);
	}
	FillHistograms(fRootHists.TrueEnergyAll_Hist,nu,intInfo,nu.Nu().E(),dirtevent);
      }

      //########################
      //### Cosmic Selection ###
      //########################

      //I really don't want to split but lets split  
      if(fConfig.IncludeCosmics || !fConfig.DirtOnly){

	for(auto const &mcs: mcshowers){ 
	  
	  //Weighting for efficiency of selection.
	  double weightcos = 1;
	  
	  ///Weight with the POT
	  weightcos *= fConfig.POTWeight;
	  
	  //Add the Global weightings if any.
	  weightcos *= fConfig.GlobalWeight;
	  
	  //Account for bnb weights?
	  weightcos *= totalbnbweight;

	  //Make a nue interaction 
	  NueSelection::NueInteraction intInfo({0, 0, 0, weightcos,-9999,-9999,-9999}); 
	  SelectCosmic(mcs, mcparticles, intInfo);
	}
      }

      return selected;
    }
    
    bool NueSelection::Select(const gallery::Event& ev, const simb::MCTruth& mctruth, const std::vector<sim::MCTrack>& mctracks, unsigned i, std::map<int, const simb::MCParticle*>& mcparticles, NueSelection::NueInteraction& intInfo){
      
      //Neutrino Info 
      const simb::MCNeutrino& nu = mctruth.GetNeutrino();

      // bool isCC;
      // if(mctruth.GetNeutrino().CCNC() == simb::kCC){
      // 	isCC = true;
      // }
      // else{
      // 	isCC = false;
      // }

      //#########################
      //### Active Volume Cut ###
      //#########################

      //Resetting direvent bool
      dirtevent = false;
      // pass activevolume cut                                                                     
      bool pass_AV = passAV(nu.Nu().Position().Vect());
      if(!pass_AV){
	if(fConfig.Verbose){std::cout << "Failed the AV Cut" << std::endl;}

	if(fConfig.IncludeDirt){dirtevent = true;}
	else{return false;}
      }
      else{
	FillHistograms(fRootHists.VisibleEnergy_AVCut_Hist,nu,intInfo,dirtevent); 
	if(fConfig.DirtOnly){return false;}
      }
      if(fConfig.Verbose){std::cout << "Passed the AV Cut" << std::endl;}

      //######################
      //### Reco Efficency ###
      //######################

      //if(isCC){
      //Check to see if the event passes the 80% ID efficiency cut from the proposal
      //bool pass_NueIDEfficiencyCut = NueIDEfficiencyCut(nu);
      // if(!pass_NueIDEfficiencyCut){
      intInfo.weight *= fConfig.nueEfficency;
      

      //#######################
      //### Beam Backgrounds ##
      //#######################
      
      //Check the Neutral current cuts. Get the photons from pi0s
      std::vector<int> pi_zeros =  findNeutralPions(mcparticles, mctruth);

      //Check the number of photons 
      std::vector<int> photons  =  findPhotons(pi_zeros, mcparticles, mctruth);

      if(fConfig.Verbose){std::cout << "Number of pi_zeros: " << pi_zeros.size() << " Number of photons: " << photons.size() << std::endl;}

      //Check to see if there is actually a shower
      if(photons.size() == 0 && intInfo.leptontrackID < 0){
	if(fConfig.Verbose){std::cout << "No lepton or photon = no shower. Event Not Selected" << std::endl;}
	return false;
      }

      //Single photons can mimic the Nue CC events they come from NC events and misidentified Numu CC events.
      if(photons.size() > 0){
	
	int photon_trackID = -1;
	//Check to see if there are two photons are above 200 MeV and so are detectable.
	bool pass_photonEnergyCut = passPhotonEnergyCut(photons,mcparticles,photon_trackID);

	if(!pass_photonEnergyCut){
	  if(fConfig.Verbose){std::cout << "Failed the Photon Energy Cut. There are more than one photons in the event. Event not Selected" << std::endl;}
	  return false;
	}
      	if(fConfig.Verbose){std::cout << "Passed the Photon Energy Cut" << std::endl;}

	//Unfortunatly we must remove CC events where a photon exists (at this stage one photon exists)
	if(intInfo.leptontrackID > 0 && photon_trackID !=-99999 && intInfo.leptonic_energy*1000 < fConfig.showerVisibleEnergyThreshold){
	  FillHistograms(fRootHists.VisibleEnergy_LeptonPlusPhotonCut_Hist,nu,intInfo,dirtevent);
	  if(fConfig.Verbose){std::cout << "Failed becuase lepton + shower exists. Event Not Selected" << std::endl;}
	  return false;
	}
	if(fConfig.Verbose){std::cout << "Passed the lepton + photon cut" << std::endl;}

	//Continue if one photon has been matched if we havn't matched then we are a CC event with no visible photons in the FV.
	if(photon_trackID !=-99999){

	  //If we have one photon candandiate choose that to be the energy of the lepton 
	  VisibleEnergyCalculator calculator;
	  calculator.lepton_pdgid = 22;
	  calculator.track_threshold =  fConfig.trackVisibleEnergyThreshold;
	  calculator.shower_energy_distortion = fConfig.showerEnergyDistortion;
	  calculator.track_energy_distortion = fConfig.trackEnergyDistortion;
	  calculator.lepton_energy_distortion_contained = fConfig.leptonEnergyDistortionContained;
	  calculator.lepton_energy_distortion_leaving_A = 1; //Not required here.
	  calculator.lepton_energy_distortion_leaving_B = 1; //Not required here. 
	  calculator.lepton_contained = true; // We don't worry at the moment if the electron is contained at this moment.
	  calculator.lepton_index = photon_trackID;

	  //Add the energy 
	  intInfo.leptonic_energy = smearLeptonEnergy(rand,mcparticles[photon_trackID],calculator);
	  intInfo.leptonpdgID = 22; 

	  FillHistograms(fRootHists.VisibleEnergy_PhotonEnergyCut_Hist,nu,intInfo,dirtevent);

	  if(!dirtevent){

	    //If hadrons produce have more than 50 MeV vertex is visible. If all photons pair produce more than 3cm away from the vertex remove.
	    bool pass_conversionGapCut  = passConversionGapCut(mcparticles,photon_trackID,intInfo.hadronic_energy,nu);
	    if(!pass_conversionGapCut){
	      if(fConfig.Verbose){std::cout << "Failed the photon conversion gap cut. Event not selected." << std::endl;}
	      return false;
	    }
	    if(fConfig.Verbose){std::cout << "Passed the Photon Conversion gap Cut" << std::endl;}
	    
	    FillHistograms(fRootHists.VisibleEnergy_ConversionGapCut_Hist,nu,intInfo,dirtevent);
	    
	    //Check if it is a numu CC where the muon is less than 1m
	    bool pass_muLenghtCut = passMuLengthCut(mctracks, mctruth);
	    if(!pass_muLenghtCut){
	      if(fConfig.Verbose){std::cout << "Failed the muon length cut. Event not selected" << std::endl;}
	      return false;
	    }
	    if(fConfig.Verbose){std::cout << "Passed the muon min length Cut" << std::endl;}
	    FillHistograms(fRootHists.VisibleEnergy_MuLenghtCut_Hist,nu,intInfo,dirtevent); 
	  }
	
	
	  //Remove 94% on a dEdx cut                                       
	  bool pass_dEdxCut = passdEdxCut(22);
	  if(!pass_dEdxCut){intInfo.weight *= fConfig.dEdxPhotonCut;}
	  FillHistograms(fRootHists.VisibleEnergy_NCCut_Hist,nu,intInfo,dirtevent);
	}


	//No photon was matched in the FV the track id is not set remove event.
	if(photon_trackID ==-99999){
	  if(fConfig.Verbose){std::cout << "No Photon in FV. Event not selected" << std::endl;}
	  return false;
	}
      }

      //##################
      //### Energy Cut ###
      //##################

      //Check the electron is not removed due to the 200 MeV cut. Just check we havn't matched 
      bool pass_eEnergyCut = passeEnergyCut(intInfo.leptonic_energy);
      if(!pass_eEnergyCut){
	if(fConfig.Verbose){std::cout << "Shower was too low in energy. Event not Selected" << std::endl;}
	return false;
      }
      FillHistograms(fRootHists.VisibleEnergy_EnergyCut_Hist,nu,intInfo,dirtevent);
      if(fConfig.Verbose){std::cout << "Passed the Energy Cut." << std::endl;}
 

      //##########################
      //### Fiducal Volume cut ###
      //##########################
      
      // pass fiducial volume cut (already passed on photon events. )                                             
      if(intInfo.leptonpdgID == 11){
	bool pass_FV = passFV(nu.Nu().Position().Vect());
	if(!pass_FV){
	  if(fConfig.Verbose){std::cout << "Failed the FV cut. Event not selected" << std::endl;}
	  return false;
	}
	if(fConfig.Verbose){std::cout << "Passed the FV Cut" << std::endl;}
	//Should have no dirts here.
	FillHistograms(fRootHists.VisibleEnergy_FVCut_Hist,nu,intInfo,dirtevent);  
      }

      return true;
    }
   
    //Check if the point is the fiducal volume.
    bool NueSelection::containedInFV(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& FV: fConfig.fiducial_volumes) {
	if (FV.Contain(p)) return true;
      }
      return false;
    }

    //Check if the point is in the Active volume.
    bool NueSelection::containedInAV(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& AV: fConfig.active_volumes) {
	if (AV.Contain(p)) return true;
      }
      return false;
    }


    //Is the neutrino a electron neutrino if so apply ID efficiency. 
    bool NueSelection::NueIDEfficiencyCut(const simb::MCNeutrino& nu){
      
      if(nu.Nu().PdgCode() == 12){
	  return true;
      }
      
      return false;
    }

    //Remove Events where the lepton is below 200 MeV in energy.
    bool NueSelection::eEnergyCut(double leptonenergy){
      
      if(fConfig.ApplyElectronEnergyCut){
	if(leptonenergy*1000 < fConfig.showerVisibleEnergyThreshold){
	  return false;
	}
      }
      return true;
    }

    //Weight events by the dEdx cut.
    bool NueSelection::dEdxCut(int pdgcode){
      
      //Apply weight for photon evenets 
      if(pdgcode==22){return false;}
      
      return true;
    }

    //Function to find the neutral pions from a neutral current interaction.
    std::vector<int> NueSelection::findNeutralPions(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth){
      
      //Pions
      std::vector<int> pions;

      //Loop over the map and identify the pions associated to the vertex. 
      for(std::map<int, const simb::MCParticle*>::iterator mcparticle=mcparticles.begin(); mcparticle!=mcparticles.end(); ++mcparticle){

	//Get the data
	const simb::MCParticle* mcparticle_data = mcparticle->second; 

	//Check the particle is a pion  
	if(mcparticle_data->PdgCode() != 111){continue;}

	//Check to see if the pion came from the vertex
	if(isFromNuVertex(mctruth,mcparticle_data)){
	  pions.push_back(mcparticle_data->TrackId());
	}
      }
      return pions;
    }

    //Function to get photons from neutral pions in a neutral current interaction 
    std::vector<int> NueSelection::findPhotons(std::vector<int>& pi_zeros, std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth){

      //Photons
      std::vector<int> photons;

      //Loop over the pions and see if they have any daughters. 
      for(int pion=0; pion<pi_zeros.size(); ++pion){
	
	//Get the pion
	const simb::MCParticle* mcparticle = mcparticles[pi_zeros.at(pion)]; 

	if(mcparticle->NumberDaughters() == 0){continue;}

	//Loop of the daughters and id photons. 
	for(int d=0; d<mcparticle->NumberDaughters(); ++d){

	  const simb::MCParticle* daughter = mcparticles[mcparticle->Daughter(d)];

	  //Are we a photon? (very rare we are not I think) or energy has to be above threshold
	  if(daughter->PdgCode() == 22 && daughter->E()*1000 > fConfig.photonVisibleEnergyThreshold){
	    photons.push_back(daughter->TrackId());
	  }
	}
      }

      //Look for any misc photons.
      for(std::map<int, const simb::MCParticle*>::iterator mcparticle=mcparticles.begin(); mcparticle!=mcparticles.end(); ++mcparticle){

	const simb::MCParticle* mcparticle_data = mcparticle->second;

	//Check the particle is a photons and is above threshold.  
	if(mcparticle_data->PdgCode() != 22 ||  mcparticle_data->E()*1000 < fConfig.photonVisibleEnergyThreshold){continue;}


	//Check we have not added the photon already.
	if(std::find(photons.begin(),photons.end(),mcparticle_data->TrackId()) != photons.end()){continue;}

	//Might be needed.  
	// //Best check that we do not come from another photon.
	// bool motherused = false;
	// int motherID = mcparticle_data->Mother();
	// const simb::MCParticle* mcparticle_mother = mcparticle_data;

	// while(motherID != 0){

	//   //The mother might not be in the list so stop
	//   if(mcparticles.find(motherID) != mcparticles.end()){break;}

	//   mcparticle_mother = mcparticles[motherID];
	//   motherID = mcparticle_mother->Mother();
	  
	//   if(std::find(photons.begin(),photons.end(),mcparticle_data->TrackId()) != photons.end()){motherused = true; break;}
	// }
	// if(motherused){continue;} 

	//Add the photon.
	//	if(isFromNuVertex(mctruth,mcparticle_mother)){
	//  photons.push_back(mcparticle_data->TrackId());
	//}

	//Add the photon 
	if(isFromNuVertex(mctruth,mcparticle_data)){   
	  photons.push_back(mcparticle_data->TrackId());
	}
      }

      return photons;
    }
    
    //Look to see if there any photons in the fiducial volume in above the energy cut. Check there is not a complementary one in the active volume.
    bool NueSelection::PhotonEnergyCut(std::vector<int>& photons, std::map<int, const simb::MCParticle*>& mcparticles, int &photonTrackID){
      
      //Initialise the photons above the cut 
      std::vector<int> photons_abovecut;
      int              photons_inAV = 0;

      for(int p=0; p<photons.size(); ++p){
	const simb::MCParticle* photon = mcparticles[photons[p]];

	//Check the photons in the active volume/
	if(containedInAV(photon->EndPosition().Vect())){
	  ++photons_inAV;
	}
		
	//Check the photons energy.  
	if(photon->E()*1000 < fConfig.showerVisibleEnergyThreshold){continue;}

	//Check if the photon is within the fiducal volume.
	if(containedInFV(photon->EndPosition().Vect())){
	  photons_abovecut.push_back(photons[p]);
	}
      }

      //Count the number of photons above threshold from the neutrino in the FV. 
      if(photons_abovecut.size() > 1){
	photonTrackID = -99999;
	return false;
      }

      //If there is no photon in the FV it still might be a nue CC event so pass -- Personally think it should have been AV 
      if(photons_abovecut.size() == 0){
	photonTrackID = -99999;
	return true;
      }

      //There might be one in the active volume as well we can remove the event if this is the case 
      if(photons_inAV > 1){
	photonTrackID = -99999;
        return false;
      }

      //Return the boolean and the track ID 
      photonTrackID = photons_abovecut[0];
      return true;
    }

    //Function to check if the photon passes the conversion gap cut. We can only get to this stage with one photon        
    bool NueSelection::ConversionGapCut(std::map<int, const simb::MCParticle*>& mcparticles, int photonTrackID, const float& hadronic_energy, const simb::MCNeutrino& nu){
      
      //Calculate the conversion lengh of the photon.First Get the end position of the photon
      TVector3 EndPos   = mcparticles[photonTrackID]->EndPosition().Vect(); 

      //Get the Vertex position
      TVector3 VtxPos = nu.Nu().Trajectory().Position(0).Vect();

      //Get the conversion length.
      double ConvLength = (EndPos - VtxPos).Mag();

      if(fConfig.Verbose){
	std::cout << "Photon Conversation length is "<< ConvLength << " cm cut is:"  <<  fConfig.photonConvLenghCut << std::endl;
	std::cout << "hadronic_energy is: " << hadronic_energy*1000 << " MeV cut is: " << fConfig.vtxEnergyCut  << std::endl;
      }

      //Check the cuts.
      if(hadronic_energy*1000 > fConfig.vtxEnergyCut && ConvLength > fConfig.photonConvLenghCut){
	return false;
      }

      return true;
    }

    //Check to see if a muon exists 
    bool NueSelection::MuLengthCut(const std::vector<sim::MCTrack>& mctrack_list, const simb::MCTruth& mctruth){

      //Find the track index which is associated to the muon 
      int track_ind = -1;
      for (int i = 0; i < mctrack_list.size(); i++) {
	if (isFromNuVertex(mctruth, mctrack_list[i]) && abs(mctrack_list[i].PdgCode()) == 13 && mctrack_list[i].Process() == "primary") {
	  if (track_ind == -1 || mctrack_list[track_ind].Start().E() < mctrack_list[i].Start().E()) {
	    track_ind = i;
	  }
	}
      }

      //No Muon. No problem pass the event. 
      if(track_ind == -1){return true;}

      //Get the track in question.
      const sim::MCTrack& track = mctrack_list[track_ind];

      //Find the length of the muon.
      double track_length = (track.Start().Position().Vect() - track.End().Position().Vect()).Mag();

      if(fConfig.Verbose){
	std::cout << "Muon Track Lengh is: " << track_length  << " cm" << std::endl;
      }

      //check if it is contained
      if(track_length > fConfig.minLengthExitingTrack || !containedInAV(track.End().Position().Vect())){return false;}
      
      return true;
    }

  
    //Fill the histograms.... can't think of a better comment.
    void NueSelection::FillHistograms(std::map<std::string,TH1D*>& HistMap, const simb::MCNeutrino& nu, NueSelection::NueInteraction& intInfo, bool& booldirtevent){
      
      double Energy = intInfo.GetNueEnergy();
      
      FillHistograms(HistMap,nu,intInfo,Energy,booldirtevent);

    }

    void NueSelection::FillHistograms(std::map<std::string,TH1D*>& HistMap, const simb::MCNeutrino& nu, NueSelection::NueInteraction& intInfo,double Energy, bool& booldirtevent){
      
      if(dirtevent){
	if(nu.CCNC() == simb::kCC){
	  
	  //Check to see if we are Numu background
	  if(nu.Nu().PdgCode() == 14){
	    std::cout << "filling dirt" << std::endl;
	    HistMap["DirtNuMu"]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["DirtInNuE"]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["DirtOscNuE"]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
	else{ 
	  //Neutral current histogtam fill
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["DirtNCNuMu"]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["DirtNCInNuE"]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["DirtNCOscNuE"]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
      }
      else{
	//Check if we are charged current
	if(nu.CCNC() == simb::kCC){
	  
	  //Check to see if we are Numu background
	  if(nu.Nu().PdgCode() == 14){
	    std::cout << "filling" << std::endl;
	    HistMap["NuMu"]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["InNuE"]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["OscNuE"]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
	else{ 
	  //Neutral current histogtam fill
	  if(nu.Nu().PdgCode() == 14){
	    HistMap["NCNuMu"]->Fill(Energy,intInfo.weight);
	    return;
	  }
	  
	  //Check if we are an intrinsic electron or oscillated 
	  if(nu.Nu().PdgCode() == intInfo.initnu){
	    //Then we are intrinsic electrons nus 
	    HistMap["NCInNuE"]->Fill(Energy,intInfo.weight);
	  }
	  else{
	    if(nu.Nu().PdgCode() != intInfo.initnu){
	      //Then we are oscilallated electrons
	      HistMap["NCOscNuE"]->Fill(Energy,intInfo.weight);
	    }
	  }
	}
      }
      return;
    }

    
    void NueSelection::PrintInformation(const simb::MCTruth& mctruth, NueSelection::NueInteraction& intInfo){
      
      std::cout << "#############################################" << std::endl;
      std::cout << "######### Truth Neutrino Information ########" << std::endl; 

      std::cout << mctruth.GetNeutrino() << std::endl;
      std::cout << "Positon X: " << mctruth.GetNeutrino().Nu().Vx() << " Y: " << mctruth.GetNeutrino().Nu().Vy() << " Z: " << mctruth.GetNeutrino().Nu().Vz() << " T: " << mctruth.GetNeutrino().Nu().T() << std::endl;
      
      std::cout << "######### Reco Neutrino Information ########"        << std::endl;
      std::cout << "Hadronic Energy:    "     << intInfo.hadronic_energy << std::endl;
      std::cout << "Shower Like Energy: "     << intInfo.shower_energy   << std::endl;
      std::cout << "Lepton Energy:          " << intInfo.leptonic_energy << std::endl;
      std::cout << "Neutrino Visible Energy " << intInfo.GetNueEnergy()  << std::endl;
      std::cout << "#############################################"       << std::endl;
      
    } 

    void NueSelection::InitialiseHistograms(){

      fEventHist = new TH1I("EventHist", "EventHist",2,0,1);
      fRootHists.VisibleEnergy_FNKPDecay_Hist = new TH1D("VisibleEnergy_FNKPDecay_Hist","VisibleEnergy_FNKPDecay_Hist",fRootHists.ebins,fRootHists.emin,fRootHists.emax); 
      fRootHists.VisibleEnergy_FNKMDecay_Hist = new TH1D("VisibleEnergy_FNKMDecay_Hist","VisibleEnergy_FNKMDecay_Hist",fRootHists.ebins,fRootHists.emin,fRootHists.emax); 
      fRootHists.VisibleEnergy_MuDecays_Hist = new TH1D(".VisibleEnergy_MuDecays_Hist",".VisibleEnergy_MuDecays_Hist",fRootHists.ebins,fRootHists.emin,fRootHists.emax);


      //Loop over the types of interactions and set the bins
      for(auto& Type: fRootHists.HistTypes){

	std::string  TrueEnergy_String                     = Type + " TrueEnergy";
	std::string  TrueEnergyAll_String                     = Type + " TrueEnergyAll";
	std::string  CCQEEnergy_String                     = Type + " CCQEEnergy";
	std::string  VisibleEnergy_String                  = Type + " VisibleEnergy";
	std::string  VisibleEnergy_AVCut_String            = Type + " VisibleEnergy_AVCut";
	std::string  VisibleEnergy_FVCut_String            = Type + " VisibleEnergy_FVCut";
	std::string  VisibleEnergy_EnergyCut_String        = Type + " VisibleEnergy_EnergyCut";
	std::string  VisibleEnergy_PhotonEnergyCut_String  = Type + " VisibleEnergy_PhotonEnergyCut";
	std::string  VisibleEnergy_ConversionGapCut_String = Type + " VisibleEnergy_ConversionGapCut";
	std::string  VisibleEnergy_MuLenghtCut_String      = Type + " VisibleEnergy_MuLenghtCut";
	std::string  VisibleEnergy_NCCut_String            = Type + " VisibleEnergy_NCCut";
	std::string  VisibleEnergy_Selection_String        = Type + " VisibleEnergy_Selection";

	std::string  VisibleEnergy_LeptonPlusPhotonCut_String        = Type + " VisibleEnergy_LeptonPlusPhotonCut";

	
	const char* TrueEnergy_Name                     = TrueEnergy_String.c_str();
	const char* TrueEnergyAll_Name                  = TrueEnergyAll_String.c_str();
	const char* CCQEEnergy_Name                     = CCQEEnergy_String.c_str();
	const char* VisibleEnergy_Name                  = VisibleEnergy_String.c_str();
	const char* VisibleEnergy_AVCut_Name            = VisibleEnergy_AVCut_String.c_str();
	const char* VisibleEnergy_FVCut_Name            = VisibleEnergy_FVCut_String.c_str();
	const char* VisibleEnergy_EnergyCut_Name        = VisibleEnergy_EnergyCut_String.c_str();
	const char* VisibleEnergy_PhotonEnergyCut_Name  = VisibleEnergy_PhotonEnergyCut_String.c_str();
	const char* VisibleEnergy_ConversionGapCut_Name = VisibleEnergy_ConversionGapCut_String.c_str();
	const char* VisibleEnergy_MuLenghtCut_Name      = VisibleEnergy_MuLenghtCut_String.c_str();
	const char* VisibleEnergy_NCCut_Name            = VisibleEnergy_NCCut_String.c_str();
	const char* VisibleEnergy_Selection_Name        = VisibleEnergy_Selection_String.c_str();

	const char* VisibleEnergy_LeptonPlusPhotonCut_Name        = VisibleEnergy_LeptonPlusPhotonCut_String.c_str();

	
	fRootHists.TrueEnergy_Hist[Type] = new TH1D(TrueEnergy_Name,TrueEnergy_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.TrueEnergyAll_Hist[Type] = new TH1D(TrueEnergyAll_Name,TrueEnergyAll_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.CCQEEnergy_Hist[Type] = new TH1D(CCQEEnergy_Name,CCQEEnergy_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_Hist[Type] = new TH1D(VisibleEnergy_Name,VisibleEnergy_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_AVCut_Hist[Type] = new TH1D(VisibleEnergy_AVCut_Name,VisibleEnergy_AVCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_FVCut_Hist[Type] = new TH1D(VisibleEnergy_FVCut_Name,VisibleEnergy_FVCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_EnergyCut_Hist[Type] = new TH1D(VisibleEnergy_EnergyCut_Name,VisibleEnergy_EnergyCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_PhotonEnergyCut_Hist[Type] = new TH1D(VisibleEnergy_PhotonEnergyCut_Name,VisibleEnergy_PhotonEnergyCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_ConversionGapCut_Hist[Type] = new TH1D(VisibleEnergy_ConversionGapCut_Name,VisibleEnergy_ConversionGapCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_MuLenghtCut_Hist[Type] = new TH1D(VisibleEnergy_MuLenghtCut_Name,VisibleEnergy_MuLenghtCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_NCCut_Hist[Type] = new TH1D(VisibleEnergy_NCCut_Name,VisibleEnergy_NCCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);
	fRootHists.VisibleEnergy_Selection_Hist[Type] = new TH1D(VisibleEnergy_Selection_Name,VisibleEnergy_Selection_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);


	fRootHists.VisibleEnergy_LeptonPlusPhotonCut_Hist[Type] = new TH1D(VisibleEnergy_LeptonPlusPhotonCut_Name,VisibleEnergy_LeptonPlusPhotonCut_Name,fRootHists.ebins,fRootHists.emin,fRootHists.emax);

      }

      fRootHists.VisibleEnergy_CosmicFVCut_Hist = new TH1D("VisibleEnergy_CosmicFVCut","VisibleEnergy_CosmicFVCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicClyinderCut_Hist = new TH1D("VisibleEnergy_CosmicClyinderCut","VisibleEnergy_CosmicClyinderCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicdEdxCut_Hist = new TH1D("VisibleEnergy_CosmicdEdxCut","VisibleEnergy_CosmicdEdxCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicWeightCut_Hist = new TH1D("VisibleEnergy_CosmicWeightCut","VisibleEnergy_CosmicWeightCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicEnergyCut_Hist = new TH1D("VisibleEnergy_CosmicEnergyCut","VisibleEnergy_CosmicEnergyCut",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
      fRootHists.VisibleEnergy_CosmicSelection_Hist = new TH1D("Cosmic VisibleEnergy_Selection","Cosmic VisibleEnergy_Selection",fRootHists.ebins,fRootHists.emin,fRootHists.emax);
    }

    bool NueSelection::SelectCosmic(const sim::MCShower& mcs, std::map<int, const simb::MCParticle*>& mcparticles, NueSelection::NueInteraction& intInfo){

      dirtevent = false;

      //Check the shower is not a neutrino shower. 
      if(mcs.Origin() ==  simb::kBeamNeutrino){return false;}
      
      if(fConfig.Verbose){std::cout << "Is a Cosmic" << std::endl;}

      //Vertex for the cosmic shower start.
      TVector3 vertex;

      //If we have one photon candndiate choose that to be the energy of the lepton 
      VisibleEnergyCalculator calculator;
      calculator.lepton_pdgid = 22;
      calculator.track_threshold =  fConfig.trackVisibleEnergyThreshold;
      calculator.shower_energy_distortion = fConfig.showerEnergyDistortion;
      calculator.track_energy_distortion = fConfig.trackEnergyDistortion;
      calculator.lepton_energy_distortion_contained = fConfig.leptonEnergyDistortionContained;
      calculator.lepton_energy_distortion_leaving_A = 1; //Not required here.
      calculator.lepton_energy_distortion_leaving_B = 1; //Not required here. 
      calculator.lepton_contained = true; // We don't worry at the moment if the electron is contained at this moment.
      calculator.lepton_index = mcs.TrackID();
      
      //Add the energy 
      intInfo.leptonic_energy = smearLeptonEnergy(rand,mcparticles[mcs.TrackID()],calculator);
      intInfo.leptonpdgID = mcs.PdgCode(); 
      
      // bool pass_eEnergyCut = passeEnergyCut(intInfo.leptonic_energy);
      // if(!pass_eEnergyCut){
      // 	if(fConfig.Verbose){std::cout << "Shower was too low in energy. Event not Selected" << std::endl;}
      // 	return false;
      // }
      
      fRootHists.VisibleEnergy_CosmicEnergyCut_Hist->Fill(intInfo.leptonic_energy);
      if(fConfig.Verbose){std::cout << "Passed the Cosmic Energy Cut." << std::endl;}
           
      //Check if it is in the fiducal volume.
      bool Pass_CosmicInFV = PassCosmicInFV(mcs,vertex);
      if(!Pass_CosmicInFV){
	if(fConfig.Verbose){std::cout << "Cosmic shower not in the FV removed" << std::endl;}
	return false;
      }
      fRootHists.VisibleEnergy_CosmicFVCut_Hist->Fill(intInfo.leptonic_energy);

      //Check if the cosmic passes the cylinder cut.
      bool Pass_CosmicCylinderCut = PassCosmicCylinderCut(mcs,vertex,mcparticles);
      if(!Pass_CosmicCylinderCut){
	if(fConfig.Verbose){std::cout << "Cosmic shower in the cyclinder" << std::endl;}
	return false;
      }
      fRootHists.VisibleEnergy_CosmicClyinderCut_Hist->Fill(intInfo.leptonic_energy);

      //dEdx weight 
      bool pass_dEdxCut = passdEdxCut(mcs.PdgCode());
      if(!pass_dEdxCut){intInfo.weight *= fConfig.dEdxPhotonCut;}
      fRootHists.VisibleEnergy_CosmicdEdxCut_Hist->Fill(intInfo.leptonic_energy);

      //Cosmic global weight for time and CRTs 
      intInfo.weight *= fConfig.CosmicWeight;
      fRootHists.VisibleEnergy_CosmicWeightCut_Hist->Fill(intInfo.leptonic_energy);
                
      return true;
    }



    //See if a shower is in the FV.
    bool NueSelection::CosmicInFV(const sim::MCShower& mcs, TVector3& vertex){
      
      if(TMath::Abs(mcs.PdgCode()) == 11){
	vertex = mcs.Start().Position().Vect();
      }
      else if(mcs.PdgCode() == 22){
	vertex = mcs.End().Position().Vect();
      }
      else{
	std::cerr << "what kind of shower is this? :S" << std::endl;
	return false;
      }
      return containedInFV(vertex);
    }

    bool NueSelection::CosmicCylinderCut(const sim::MCShower& mcs, TVector3& vertex, std::map<int, const simb::MCParticle*>& mcparticles){
      
      //Get the mother track. If it doesn't have a mother we cannot do this. 
      if(mcs.MotherTrackID() == -1){return true;}
      
      const simb::MCParticle* mother = mcparticles[mcs.MotherTrackID()];
      if(fConfig.Verbose){std::cout << "cosmic mother is pdgcode: " << mother->PdgCode() << " track id: " << mother->TrackId() << " and energy: " << mother->E() << std::endl;} 
      
      //Cant do this if its not a muon 
      if(mother->PdgCode() != 13){return true;}
      
      //Muon will never be not in the FV with the shower being in this case
      //Check if the vertex is within a cylinder of the muon. Assuming the muon is just a straight line then work out mag from start to vector and star and end
      TVector3 StartPos   = mother->Position().Vect();
      TVector3 EndPos     = mother->EndPosition().Vect();
      
      double   StartVtxMag = (StartPos - vertex).Mag();
      double   StartEndMag = (StartPos - EndPos).Mag();
	   
      //What is the angle
      double AngleStart = (vertex - StartPos).Angle(EndPos - StartPos);
      if(AngleStart > TMath::Pi()/2){
	//angle is above the start of the track
	std::cout << "Dom Check the angle" << std::endl;
	if(fConfig.Verbose){std::cout << "Shower above start of the track" << std::endl;}
	return true;
      }
	  
      double AngleEnd = (vertex - EndPos).Angle(StartPos - EndPos);
      if(AngleEnd > TMath::Pi()/2){
	//angle is below the end of the track                                                                                                                              
	if(fConfig.Verbose){std::cout << "Shower below end of the track" << std::endl;}
	std::cout << "Dom Check the angle for end" << std::endl;
	return true;
      }
      
      //Simple trig to work out the distance from the vertex to the track now.
      double DistFromMother = StartVtxMag*TMath::Sin(AngleStart);
      if(fConfig.Verbose){std::cout << "Angle from the start is: " << AngleStart << " Hypot is: " << StartVtxMag << std::endl;} 
      
      if(DistFromMother > fConfig.CosmicVolumeRadius){return true;}
      
      return false;
    }
    
  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelection)

