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

//SBN Includes 
#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

    NueSelection::NueSelection() : SelectionBase(), EventCounter(0), NuCount(0) {}

    
    void NueSelection::Initialize(fhicl::ParameterSet* config) {


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

      fConfig.ApplyElectronEnergyCut = pconfig.get<bool>("ApplyElectronEnergyCut", true);
      fConfig.NueIDEfficiency        = pconfig.get<bool>("NueIDEfficiency",true);       
      fConfig.doFVCut                = pconfig.get<bool>("doFVCut",true);

      hello();
    }
    
    
    void NueSelection::Finalize() {}
    
    
    bool NueSelection::ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco){
      if (EventCounter % 10 == 0) {
        std::cout << "NueSelection: Processing event " << EventCounter << " "
                  << "(" << NuCount << " neutrinos selected)"
                  << std::endl;
      }
      EventCounter++;

      std::cout << "New Event" << std::endl;

      
      //Grab a neutrino datat product from the event
      auto const& mctruths =  *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
      
      //Get tracks and showers 
      auto const& mctracks =  *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
      auto const& mcshowers = *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);
      
      //Get a map of the all the particles (not just the final state ones).
      std::map<int,const simb::MCParticle*> mcparticles;
      auto const& mcparticle_list = *ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);
      for(auto const &mcparticle:  mcparticle_list){
	
	trueParticles[mcparticle.TrackId()] = &mcparticle;

	if(fConfig.Verbose){
	  std::cout << "MC Particle with track ID: " << mcparticle.TrackId() << " has pdgcode: " << mcparticle.PdgCode() << " with energy: " << mcparticle.E() << std::endl;
	}
      }


      //Iterate through the neutrinos
      for (size_t i=0; i<mctruths.size(); i++) {
        auto const& mctruth = mctruths.at(i);
        const simb::MCNeutrino& nu = mctruth.GetNeutrino();

	//Set the values required to calculate the track and shower energys + distortion 
	VisibleEnergyCalculator calculator;
	calculator.lepton_pdgid = 13;
	calculator.track_threshold =  fConfig.trackVisibleEnergyThreshold;
	calculator.shower_energy_distortion = fConfig.showerEnergyDistortion;
	calculator.track_energy_distortion = fConfig.trackEnergyDistortion;
	calculator.lepton_energy_distortion_contained = fConfig.leptonEnergyDistortionContained;
	calculator.lepton_energy_distortion_leaving_A = fConfig.leptonEnergyDistortionLeavingA;
	calculator.lepton_energy_distortion_leaving_B = fConfig.leptonEnergyDistortionLeavingB;	

	//Check if the neutrino is contained.
	bool containedinfv = containedInFV(nu.Nu().Position().Vect());
	calculator.lepton_contained = containedinfv;

	if(!containedinfv){return false;}

	std::cout << "\n\nINTERACTION:\n";
	std::cout << "MODE: " << mctruth.GetNeutrino().Mode() << std::endl;
	std::cout << "CC: " << mctruth.GetNeutrino().CCNC() << std::endl;

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


	for(int i=0;i< mctruth.NParticles();++i){
	  std::cout << "Particle with track ID: " <<  mctruth.GetParticle(i).TrackId() << " has pdgcode: " << mctruth.GetParticle(i).PdgCode() << " with energy: " << mctruth.GetParticle(i).E() << " daughter num: " << mctruth.GetParticle(i).NumberDaughters() << "FirstDaughter: " << mctruth.GetParticle(i).FirstDaughter() << " last: " <<  mctruth.GetParticle(i).LastDaughter() << std::endl;
	}

	// build the interaction
        Event::Interaction interaction = truth[i];

	//Calculate the Energy 
	const std::vector<double> visible_energy = FlavourEnergyDeposition(mctruth, mctracks, mcshowers,calculator);
	double Hadronic_energy = visible_energy[0];
	double Shower_energy   = visible_energy[1];
	double Leptonic_energy = visible_energy[2];
	double Nue_energy = Hadronic_energy + ShowerEnergy + LeptonicEnergy;

	//Weighting for efficiency of selection.
	double weight = 1;

	//Make a nue interaction 
	NueInteraction intInfo(Hadronic_energy,Shower_energy,Leptonic_energy,Nue_energy,weight); 

	//Make the reco interaction event. 
	Event::RecoInteraction reco_interaction(interaction, i);
        reco_interaction.reco_energy = visible_energy[0];

        //Run the Selection
	bool selection = Select(ev, mctruth, mctracks, fConfig, mcparticles, weight);

        //Check for Tracks - Remove any track that come from the neutrino if 1m 
        //Find electron vertex candiates - Remove 80% due to fake dEdx cut 
        //If pi0 is made 2 photons come out as candiates in both are above 100 MeV remove.
        //If hadrons produce have more than 50 MeV this is visible if all photons pair produce more than 3cm away from the vertex remove.
        //Remove 94% pion of events as dEdx cut
        //Event could be CC0Pi Muons Interaction with photon background which mimics the CC1Pi Electron interatactions. these need to go throught the cuts above as well.  
      
        if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
          Event::RecoInteraction interaction(truth[i], i);
          reco.push_back(interaction);
        }
      }
      
      if (true) {
        NuCount++;
      }
      return true;
    }
    
    bool NueSelection::Select(const gallery::Event& ev, const simb::MCTruth& mctruth, std::vector<const simb::MCTrack>& mctracks, unsigned i, std::map<int, const simb::MCParticle*>& mcparticles, NueInteraction& intInfo){
      
      //Neutrino Info 
      const simb::MCNeutrino& nu = mctruth.GetNeutrino();

      //#########################
      //### Fiducal Volume cut ##
      //#########################

      // pass fiducial volume cut                                                                                                                           
      bool pass_FV = passFV(nu.Nu().Position().Vect());
      if(!pass_FV){return false;}
      
      
      bool isCC;
      if(mctruth.GetNeutrino().CCNC() == simb::kCC){
	isCC = true;
      }
      else{
	isCC = false;
      }

      //##################
      //### Signal Nue ###
      //##################

      if(isCC){

	//Check to see if the event passes the 80% ID efficiency cut from the proposal
	bool pass_NueIDEfficiencyCut = passNueIDEfficiencyCut(intInfo.weight);
	if(!pass_NueIDEfficiencyCut){weight *= fConfig.nueEfficencyWeight}
	
	//Check the electron is not removed due to the 200 MeV cut. 
	bool pass_eEnergyCut = passeEnergyCut(intInfo.leptonic_energy);
	if(!pass_eEnergyCut){return false;}
      }
      
      //#######################
      //### Beam Backgrounds ##
      //#######################
      
      //Check the Neutral current cuts. Get the photons from pi0s
      std::vector<int> pi_zeros =  findNeutralPions(mcparticles, mctruth);
      //Check the number of photons 
      std::vector<int> photons  =  findPhotons(pi_zeros,mcparticles);

      //Single photons can mimic the Nue CC events they come from NC events and misidentified Numu CC events.
      if(photons.size() > 0){
	
	int photon_trackID;
	//Check to see if there are two photons are above 200 MeV and so are detectable.
	bool pass_photonEnergyCut = passPhotonEnergyCut(photons,mcparticles,photon_trackID);
	if(pass_secondphotonEnergyCut){return false;}

	//If hadrons produce have more than 50 MeV vertex is visible. If all photons pair produce more than 3cm away from the vertex remove.
	bool pass_conversionGapCut      = passConversionGapCut(photon_trackID,mcparticles,hadronic_energy);
	if(pass_conversionGapCut){return false;}

	//Remove 94% on a dEdx cut 
	bool pass_dEdxCut = passdEdxCut(nu);
	if(pass_dEdxCut){weight *= fConfig.dEdxPhotonCut;}

	//Check if it is a numu CC where the muon is less than 1m
	bool pass_muLenghtCut = passMuLengthCut(mctracks);
	if(pass_muLenghtCut){return false;}

	//Things to do. Set booleans functions to turn off cuts.
	//Add histograms.
	//Correct the lepton energy calculation. 
	//Set the reco interaction info. 

      }
      
      return true;
    }
  
    bool NueSelection::containedInFV(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& FV: fConfig.fiducial_volumes) {
	if (FV.Contain(p)) return true;
      }
      return false;
    }

    bool NueSelection::containedInAV(const TVector3 &v) {
      geoalgo::Point_t p(v);
      for (auto const& AV: fConfig.active_volumes) {
	if (AV.Contain(p)) return true;
      }
      return false;
    }


    //Is the neutrino a electron neutrino if so apply ID efficiency 
    bool NueSelection::passNueIDEfficiencyCut(const simb::MCNeutrino& nu){
      
      if(nu.Nu().PdgCode() == 12){
	  return false;
      }
      
      return true;
    }

    //Remove Events below 200 MeV.
    bool NueSelection::eEnergyCut(double leptonenergy){
      
      if(nu.Nu().PdgCode() == 12 && fConfig.ApplyElectronEnergyCut){
	if(leptonenergy < 200){
	  return false;
	}
      }
      return true;
    }

    //Weight events by the dEdx cut.
    bool NueSelection::passdEdxCut(const simb::MCNeutrino& nu){
      
      if(nu.CCNC() == simb::kCC){
	return false;
      }
      return true;
    }

    //Function to find the neutral pions from a neutral current interaction.
    std::vector<int> NueSelection::findNeutralPions(std::map<int, const simb::MCParticle*>& mcparticles, const simb::MCTruth& mctruth){
      
      //Pions
      std::vector<int> pions

      //Loop over the map and identify the pions associated to the vertex. 
      for(std::map<int, const simb::MCParticle*>::iterator mcparticle=mcparticles.begin(); mcparticle!=mcparticles.end(); ++mcparticle){

	//Get the data
	const simb::MCParticle* mcparticle_data = mcparticle->second; 

	//Check the particle is a pion  
	if(mcparticle_data->PdgCode() != 111){continue;}

	//Check to see if the pion came from the vertex
	if(isFromNuVertex(mctruth,mcparticle)){
	  pions.push_back(mcparticle_data->TrackId());
	}
      }
      return pions;
    }

    //Function to get photons from neutral pions in a neutral current interaction 
    std::vector<int> NueSelection::findPhotons(std::vector<int> pi_zeros, std::map<int, const simb::MCParticle*>& mcparticles){

      //Photons
      std::vector<int> photons;

      //Loop over the pions and see if they have any daughters. 
      for(int pion=0; pion<pions.size(); ++pion){
	
	//Get the pion
	const simb::MCParticle* mcparticle[pi_zeros.at(pion)]; 

	if(mcparticle->NumberDaughters() == 0){continue;}

	//Loop of the daughters and id photons. 
	for(int d=0; d<mcparticle->NumberDaughters(); ++d){

	  const simb::MCParticle* daughter = mcparticles[mcparticle->Daughter(d)];

	  //Are we a photon? (very rare we are not I think) or energy has to be above threshold
	  if(daughter->PdgCode() == 22 || daugther->E() > fConfig.ShowerEnergyThreshold){
	    photons.push_back(daughter->trackId());
	  }
	}
      }

      //Look for any misc photons.
      for(std::map<int, const simb::MCParticle*>::iterator mcparticle=mcparticles.begin(); mcparticle!=mcparticles.end(); ++mcparticle){

	const simb::MCParticle* mcparticle_data = mcparticle->second;

	//Check the particle is a photons and is above threshold.  
	if(mcparticle_data->PdgCode() != 22 || daugther->E() > fConfig.ShowerEnergyThreshold){continue;}

	//Check we have not added the photon already.
	if(std::find(photons.begin,photons.end(),mcparticle_data->TrackId()) != photons.end()){continue;}

	//Best check that we do not come from another photon.
	bool motherused = false;
	int motherID = mcparticle_data->Mother();
	const simb::MCParticle* mcparticle_mother = mcparticle_data;
	while(motherID != 0){
	  mcparticle_mother = mcparticles[motherID];
	  motherID = mcparticle_mother->Mother();
	  if(std::find(photons.begin,photons.end(),mcparticle_data->TrackId()) != photons.end()){motherused = true; break;}
	}
	if(motherused){continue;} 

	//Add the photon.
	if(isFromNuVertex(mctruth,mcparticle_mother)){
	  photons.push_back(mcparticle->TrackID());
	}
      }
      return photons;
    }

    bool NueSelection::passPhotonEnergyCut(std::vector<int> photons, std::map<int, const simb::MCParticle*>& mcparticles, int& photonTrackID){
      
      //Initialise the photons above the cut 
      std::vector<int> photons_abovecut;
      int              photons_inAV;

      for(int p=0; p<photons.size(); ++p){
	const simb::MCParticle* photon = mcparticles[p];
	
	//Check the photons energy.  
	if(photon->E() < fConfig.photonEnergyCut){continue;}

	//Check if the photon is within the fiducal volume.
	if(containedInFV(photon.Position().Vect())){
	  photons_abovecut.push_back(photons[p]);
	}
	
	//Check the photons in the active volume/
	if(containedInAV(photon.Position().Vect())){
	  ++photons_inAV;
	}
      }

      //Count the number of photons above threshold from the neutrino in the FV. 
      if(photons_abovecut_inFV.size() != 1){
	photontrackID = -99999;
	return true;
      }

      //There might be one in the active volume as well we can remove the event if this is the case 
      if(photons_inAV > 1){
	photontrackID = -99999;
        return true;
      }

      //Return the boolean and the track ID 
      photontrackID = photons_abovecut[0];
      return false;
    }

    //Function to check if the photon passes the conversion gap cut. We can only get to this stage with one photon        
    bool NueSelection::passConversionGapCut(std::map<int, const simb::MCParticle*>& mcparticles, int& photonTrackID, const int& hadronicE, const simb::MCNeutrino& nu){
      
      //Calculate the conversion lengh of the photon.First Get the end position of the photon
      TVector3 EndPos   = mcparticles[photonTrackID]->EndPosition().Vect(); 

      //Get the Vertex position
      TVector3 VtxPos = nu.Trajectory().Position(0).Vect();

      //Get the conversion length.
      double ConvLength = (EndPos - VtxPos).Mag();

      //Check the cuts.
      if(hardronicE > fConfig.VtxEnergyCut && ConvLength > fConfig.PhotonConvLenghCut){
	return true;
      }

      return false;
    }

    //Check to see if a muon exists 
    bool NueSelection::passMuLengthCut(std::vector<const sim::MCTrack>& mctrack_list){

      //Find the track index which is associated to the muon 
      int track_ind = -1;
      for (int i = 0; i < mctrack_list.size(); i++) {
	if (isFromNuVertex(mctruth, mctrack_list[i]) && abs(mctrack_list[i].PdgCode()) == 13 && mctrack_list[i].Process() == "primary") {
	  if (track_ind == -1 || mctrack_list[track_ind].Start().E() < mctrack_list[i].Start().E()) {
	    track_ind = i;
	  }
	}
      }

      if(track_ind = -1){return true;}

      //Get the track in question.
      const simb::MCTruth& track = mctrack_list[track_ind];

      //Find the length of the muon.
      double track_length = (track.Start().Position().Vect() - track.End().Position().Vect()).Mag();
      
      if(track_length > fConfig.minLengthExitingTrack){return true;}
      
      return false;
    }

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelection)

