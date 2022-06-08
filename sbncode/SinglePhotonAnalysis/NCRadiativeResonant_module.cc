////////////////////////////////////////////////////////////////////////
// Class:       NCRadiativeResonant
// Plugin Type: filter (art v2_05_00)
// File:        NCRadiativeResonant_module.cc
//
// Generated at Fri Jun 23 10:33:44 2017 by Robert Murrells using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////
// Disable the FillTree function  
//			by Keng Jun. 2022


#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "TTree.h"

#include <memory>

class NCRadiativeResonant : public art::EDFilter {

  TTree * ftree;

  int frun;
  int fsubrun;
  int fevent;
  int fnu_pdg;
  int fccnc;
  int fmode;
  int finteraction_type;
  int fis_nc_delta_radiative;
  int fparent_status_code;
  int fparent_pdg;

public:
  explicit NCRadiativeResonant(fhicl::ParameterSet const & p);

  NCRadiativeResonant(NCRadiativeResonant const &) = delete;
  NCRadiativeResonant(NCRadiativeResonant &&) = delete;
  NCRadiativeResonant & operator = (NCRadiativeResonant const &) = delete;
  NCRadiativeResonant & operator = (NCRadiativeResonant &&) = delete;

  void cout_stuff(art::Event & e, bool passed);
  void FillTree(art::Event & e, 
		size_t const mct_index,
		size_t const parent_index,
		bool const is_nc_delta_radiative);
  void Reset();
  bool filter(art::Event & e) override;

};


//CHECK ????? FIX THIS?
NCRadiativeResonant::NCRadiativeResonant(fhicl::ParameterSet const & p) :
  art::EDFilter(p),
    ftree(nullptr) {

  if(true) {

    art::ServiceHandle<art::TFileService> tfs;
    ftree = tfs->make<TTree>("NCRadiativeResonantFilter", "");

    ftree->Branch("run", &frun, "run/I");
    ftree->Branch("subrun", &fsubrun, "subrun/I");
    ftree->Branch("event", &fevent, "event/I");
    ftree->Branch("nu_pdg", &fnu_pdg, "nu_pdg/I");
    ftree->Branch("ccnc", &fccnc, "ccnc/I");
    ftree->Branch("mode", &fmode, "mode/I");
    ftree->Branch("is_nc_delta_radiative", &fis_nc_delta_radiative, "is_nc_delta_radiative/I");
    ftree->Branch("parent_status_code", &fparent_status_code, "parent_status_code/I");
    ftree->Branch("parent_pdg", &fparent_pdg, "parent_pdg/I");

  }

}


//select resonant events, it is not necessary a delta event.
void NCRadiativeResonant::cout_stuff(art::Event & e, bool passed = false) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  std::cout << passed << "\n"
	    << "===========Summary of Truth info===============\n";
  for(simb::MCTruth const & mct : *ev_mct) {
      std::cout<<std::setw(9)<<"track ID ";
	  std::cout<<std::setw(12)<<"Pdgcode ";
	  std::cout<<std::setw(11)<<"mother ID ";
	  std::cout<<std::setw(11)<<"Status code "<<std::endl;
	  std::cout << "----------------------------\n";
    for(int i = 0; i < mct.NParticles(); ++i) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      std::cout<<std::setw(9)<< mcp.TrackId(); 
	  std::cout<<std::setw(12)<< mcp.PdgCode();
	  std::cout<<std::setw(11)<< mcp.Mother();
	  std::cout<<std::setw(11)<< mcp.StatusCode() << "\n";
    }
  }

}


void NCRadiativeResonant::Reset() {

  frun = -1;
  fsubrun = -1;
  fevent = -1;
  fnu_pdg = 0;
  fccnc = -1;
  fmode = -2;
  finteraction_type = -2;
  fis_nc_delta_radiative = -1;
  fparent_status_code = -1; 
  fparent_pdg = 0;

}


void NCRadiativeResonant::FillTree(art::Event & e, 
				size_t mct_index,
				size_t const parent_index,
				bool const is_nc_delta_radiative) {

  Reset();

  frun = e.id().run();
  fsubrun = e.id().subRun();
  fevent = e.id().event();

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
	  e.getValidHandle<std::vector<simb::MCTruth>>("generator");
  simb::MCNeutrino const & mcn = ev_mct->at(mct_index).GetNeutrino(); 

  fnu_pdg = mcn.Nu().PdgCode();
  fccnc = mcn.CCNC();
  fmode = mcn.Mode();
  finteraction_type = mcn.InteractionType();

  fis_nc_delta_radiative = is_nc_delta_radiative;

  if(parent_index != SIZE_MAX) {
    fparent_status_code = ev_mct->at(mct_index).GetParticle(parent_index).StatusCode();
    fparent_pdg = ev_mct->at(mct_index).GetParticle(parent_index).PdgCode();
  }    

  ftree->Fill();

}


//true - pass; false - reject;
bool NCRadiativeResonant::filter(art::Event & e) {

	art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
		e.getValidHandle<std::vector<simb::MCTruth>>("generator");


	//looking for resonant photons...
	//	- Particles from NC interaction
	//	- Particles contian photons with StatusCode 1
	//		- the parent of a photon is not delta0 & delta+ ---> pass
	//		- if the parent of a photon is a photon
	//			- the grandparent of a photon is not delta0 & delta+ ---> pass
	for(size_t ith = 0; ith < ev_mct->size(); ++ith) {//loop through MCTruth 

		simb::MCTruth const & mct = ev_mct->at(ith);
//		std::cout<<"Looking at the "<<ith<<"th mctruth particle."<<std::endl;
		if(mct.GetNeutrino().CCNC() != 1) continue;//only look at nc particles;

		std::vector<size_t> exiting_photon_parents;
		for(int jth = 0; jth < mct.NParticles(); ++jth) {//loop through MCParticles
			simb::MCParticle const & mcp = mct.GetParticle(jth);
			
			//CHEKC hardcode, TPC filter:
			if(abs(mcp.Vx())>210 ||  abs(mcp.Vy())>210||mcp.Vz()>510 || mcp.Vz()<-1){
			std::cout<<"OUTSIDE TPC x y z ="<<mcp.Vx()<<" "<<mcp.Vy()<<" "<<mcp.Vz()<<std::endl;
			exit(0);
			}
			
			if(mcp.TrackId() != jth) {
				std::cout << "ERROR: " << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nTrackId does not match index\n";
				exit(1);
			}
			if(!(mcp.StatusCode() == 1 && mcp.PdgCode() == 22)) continue;//next, if this is not a photon with StatusCode 1
			exiting_photon_parents.push_back(mcp.Mother());//collect photon's parents
		}
		
		//looking for delta radiative decay, i.e. mother is a delta with statuscode 3 & daughter is a photon. We don't need exiting deltas.
		std::vector<size_t> in_nucleus_photons;
		for(size_t const s : exiting_photon_parents) {
			simb::MCParticle const & mcp = mct.GetParticle(s);
			//2114 delta0; 2214 delta+
			if(abs(mcp.PdgCode()) != 2114 && abs(mcp.PdgCode()) != 2214 ) {//Want mcp as any particles but not delta0+
//				if(ftree){ 
//					FillTree(e, ith, mcp.Mother(), true);
//				}
				cout_stuff(e,true);
				std::cout<<"\nYES a resonant evt: "<<mcp.PdgCode()<<std::endl;
				return true;
			}else if(mcp.PdgCode() == 22) {//collect photons;
				in_nucleus_photons.push_back(mcp.Mother());
			}
		}

		for(size_t const s : in_nucleus_photons) {
			simb::MCParticle const & mcp = mct.GetParticle(s);
			if(abs(mcp.PdgCode()) != 2114 && abs(mcp.PdgCode()) != 2214) {
//				if(ftree) FillTree(e, ith, s, true);
				cout_stuff(e,true);
				std::cout<<"\nYES a resonant evt with 2 photons: "<<mcp.PdgCode()<<std::endl;
				return true;
			}
		}

	}

//	if(ftree) FillTree(e, 0, SIZE_MAX, false);
	return false;

}


DEFINE_ART_MODULE(NCRadiativeResonant)
