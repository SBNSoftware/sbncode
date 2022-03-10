////////////////////////////////////////////////////////////////////////
// Class:       NCRadiativeResonant
// Plugin Type: filter (art v2_05_00)
// File:        NCRadiativeResonant_module.cc
//
// Generated at Fri Jun 23 10:33:44 2017 by Robert Murrells using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

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


NCRadiativeResonant::NCRadiativeResonant(fhicl::ParameterSet const & p) :
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


void NCRadiativeResonant::cout_stuff(art::Event & e, bool passed = false) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  std::cout << passed << "\n"
	    << "==========================\n";
  for(simb::MCTruth const & mct : *ev_mct) {
    std::cout << "----------------------------\n";
    for(int i = 0; i < mct.NParticles(); ++i) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      std::cout <<"FULL: "<< mcp.TrackId() << " " << mcp.PdgCode() << " " << mcp.Mother() << " " << mcp.StatusCode() << "\n";
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


bool NCRadiativeResonant::filter(art::Event & e) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  cout_stuff(e,true);

  for(size_t i = 0; i < ev_mct->size(); ++i) {

    simb::MCTruth const & mct = ev_mct->at(i);
    if(mct.GetNeutrino().CCNC() != 1) continue;

    std::vector<size_t> exiting_photon_parents;
    for(int i = 0; i < mct.NParticles(); ++i) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      if(mcp.TrackId() != i) {
	std::cout << "ERROR: " << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nTrackId does not match index\n";
	exit(1);
      }
      if(!(mcp.StatusCode() == 1 && mcp.PdgCode() == 22)) continue;
      exiting_photon_parents.push_back(mcp.Mother());
    }

    std::vector<size_t> in_nucleus_photons;
    for(size_t const s : exiting_photon_parents) {
      simb::MCParticle const & mcp = mct.GetParticle(s);
      if(abs(mcp.PdgCode()) != 2114 && abs(mcp.PdgCode()) != 2214 ) {
	if(ftree) FillTree(e, i, mcp.PdgCode(), true);
    std::cout<<"YES: "<<mcp.PdgCode()<<std::endl;
	return true;
      }
      else if(mcp.PdgCode() == 22) {
	in_nucleus_photons.push_back(mcp.Mother());
      }
    }

    for(size_t const s : in_nucleus_photons) {
      simb::MCParticle const & mcp = mct.GetParticle(s);
      if(abs(mcp.PdgCode()) != 2114 && abs(mcp.PdgCode()) != 2214) {
	if(ftree) FillTree(e, i, s, true);
    std::cout<<"YES: "<<mcp.PdgCode()<<std::endl;
	return true;
      }
    }

  }

  if(ftree) FillTree(e, 0, SIZE_MAX, false);
  return false;

}


DEFINE_ART_MODULE(NCRadiativeResonant)
