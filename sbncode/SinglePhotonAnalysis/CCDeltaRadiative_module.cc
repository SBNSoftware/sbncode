////////////////////////////////////////////////////////////////////////
// Class:       CCDeltaRadiative
// Plugin Type: filter (art v2_05_00)
// File:        CCDeltaRadiative_module.cc
//
// Selects CC events with a photon from a Delta radiative decay.
// Mirrors NCDeltaRadiative_module.cc but requires CCNC == 0 (charged
// current) and checks all four Delta charge states (1114, 2114, 2214, 2224).
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
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

#include <memory>

#include "nusimdata/SimulationBase/MCTruth.h"

#include "TTree.h"

class CCDeltaRadiative : public art::EDFilter {

  TTree * ftree;

  int frun;
  int fsubrun;
  int fevent;
  int fnu_pdg;
  int fccnc;
  int fmode;
  int finteraction_type;
  int fis_cc_delta_radiative;
  int fparent_status_code;
  int fparent_pdg;

public:
  explicit CCDeltaRadiative(fhicl::ParameterSet const & p);

  CCDeltaRadiative(CCDeltaRadiative const &) = delete;
  CCDeltaRadiative(CCDeltaRadiative &&) = delete;
  CCDeltaRadiative & operator = (CCDeltaRadiative const &) = delete;
  CCDeltaRadiative & operator = (CCDeltaRadiative &&) = delete;

  void cout_stuff(art::Event & e, bool passed);
  void FillTree(art::Event & e,
    size_t const mct_index,
    size_t const parent_index,
    bool const is_cc_delta_radiative);
  void Reset();
  bool filter(art::Event & e) override;

};


CCDeltaRadiative::CCDeltaRadiative(fhicl::ParameterSet const & p) :
  art::EDFilter(p),
  ftree(nullptr) {

    art::ServiceHandle<art::TFileService> tfs;
    ftree = tfs->make<TTree>("CCDeltaRadiativeFilter", "");

    ftree->Branch("run", &frun, "run/I");
    ftree->Branch("subrun", &fsubrun, "subrun/I");
    ftree->Branch("event", &fevent, "event/I");
    ftree->Branch("nu_pdg", &fnu_pdg, "nu_pdg/I");
    ftree->Branch("ccnc", &fccnc, "ccnc/I");
    ftree->Branch("mode", &fmode, "mode/I");
    ftree->Branch("is_cc_delta_radiative", &fis_cc_delta_radiative, "is_cc_delta_radiative/I");
    ftree->Branch("parent_status_code", &fparent_status_code, "parent_status_code/I");
    ftree->Branch("parent_pdg", &fparent_pdg, "parent_pdg/I");

}


void CCDeltaRadiative::cout_stuff(art::Event & e, bool passed = false) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  std::cout << passed << "\n"
      << "==========================\n";
  for(simb::MCTruth const & mct : *ev_mct) {
    std::cout << "----------------------------\n";
    for(int i = 0; i < mct.NParticles(); ++i) {
      simb::MCParticle const & mcp = mct.GetParticle(i);
      std::cout << mcp.TrackId() << " " << mcp.PdgCode() << " " << mcp.Mother() << " " << mcp.StatusCode() << "\n";
    }
  }

}


void CCDeltaRadiative::Reset() {

  frun = -1;
  fsubrun = -1;
  fevent = -1;
  fnu_pdg = 0;
  fccnc = -1;
  fmode = -2;
  finteraction_type = -2;
  fis_cc_delta_radiative = -1;
  fparent_status_code = -1;
  fparent_pdg = 0;

}


void CCDeltaRadiative::FillTree(art::Event & e,
        size_t const mct_index,
        size_t const parent_index,
        bool const is_cc_delta_radiative) {

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

  fis_cc_delta_radiative = is_cc_delta_radiative;
  if(parent_index != SIZE_MAX) {
    fparent_status_code = ev_mct->at(mct_index).GetParticle(parent_index).StatusCode();
    fparent_pdg = ev_mct->at(mct_index).GetParticle(parent_index).PdgCode();
  }

  ftree->Fill();

}


// Returns true if |pdg| matches any Delta charge state (Δ-, Δ0, Δ+, Δ++).
static bool isDelta(int pdg) {
  int apdg = std::abs(pdg);
  return apdg == 1114 || apdg == 2114 || apdg == 2214 || apdg == 2224;
}


bool CCDeltaRadiative::filter(art::Event & e) {

  art::ValidHandle<std::vector<simb::MCTruth>> const & ev_mct =
    e.getValidHandle<std::vector<simb::MCTruth>>("generator");

  for(size_t i = 0; i < ev_mct->size(); ++i) {

    simb::MCTruth const & mct = ev_mct->at(i);
    // Require charged-current interaction (CCNC == 0).
    if(mct.GetNeutrino().CCNC() != 0) continue;

    std::vector<size_t> exiting_photon_parents;
    for(int j = 0; j < mct.NParticles(); ++j) {
      simb::MCParticle const & mcp = mct.GetParticle(j);
      if(mcp.TrackId() != j) {
        std::cout << "ERROR: " << __LINE__ << " " << __PRETTY_FUNCTION__ << "\nTrackId does not match index\n";
        return false;
      }
      if(!(mcp.StatusCode() == 1 && mcp.PdgCode() == 22)) continue;
      exiting_photon_parents.push_back(mcp.Mother());
    }

    std::vector<size_t> in_nucleus_photons;
    for(size_t const s : exiting_photon_parents) {
      simb::MCParticle const & mcp = mct.GetParticle(s);
      // Require parent vertex inside TPC.
      if(std::abs(mcp.Vx()) > 210 || std::abs(mcp.Vy()) > 210 || mcp.Vz() > 510 || mcp.Vz() < -1) {
        std::cout << "OUTSIDE TPC x y z =" << mcp.Vx() << " " << mcp.Vy() << " " << mcp.Vz() << std::endl;
        return false;
      }
      if(isDelta(mcp.PdgCode())) {
        if(ftree) FillTree(e, i, s, true);
        return true;
      }
      else if(mcp.PdgCode() == 22) {
        in_nucleus_photons.push_back(mcp.Mother());
      }
    }

    for(size_t const s : in_nucleus_photons) {
      simb::MCParticle const & mcp = mct.GetParticle(s);
      if(isDelta(mcp.PdgCode())) {
        if(ftree) FillTree(e, i, s, true);
        return true;
      }
    }

  }

  if(ftree) FillTree(e, 0, SIZE_MAX, false);
  return false;

}


DEFINE_ART_MODULE(CCDeltaRadiative)
