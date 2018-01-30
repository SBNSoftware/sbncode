#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// includes for random draws from gaussian
#include <random>

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "uboone/EventWeight/MCEventWeight.h"
#include "gallery/Event.h"

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TTree.h>

#include "TSUtil.h"
#include "TSSelection.h"

namespace ana {
namespace lee_truth_selection {
void TSSelection::setShowerEnergyResolution(float res, bool by_percent) {
  _shower_energy_resolution = res;
  _shower_energy_by_percent = by_percent;
  _shower_energy_distribution = std::normal_distribution<float>(0., res);
}

void TSSelection::setTrackEnergyResolution(float res, bool by_percent) {
  _track_energy_resolution = res;
  _track_energy_by_percent = by_percent;
  _track_energy_distribution = std::normal_distribution<float>(0., res);
}

void TSSelection::setShowerAngleResolution(float res, bool by_percent) {
  _shower_angle_resolution = res;
  _shower_angle_by_percent = by_percent;
  _shower_angle_distribution = std::normal_distribution<float>(0., res);
}

void TSSelection::setTrackAngleResolution(float res, bool by_percent) {
  _track_angle_resolution = res;
  _track_angle_by_percent = by_percent;
  _track_angle_distribution = std::normal_distribution<float>(0., res);
}

// the number of protons you pass in is the signal region that will be changed
// options:
// zero protons, one proton, any number of protons
void TSSelection::setAcceptP(bool b, int n_protons) {
  if (n_protons == 1)
    _accept_1p = b;
  else if (n_protons > 1)
    _accept_np = b;
  else if (n_protons == 0)
    _accept_0p = b;
}

float TSSelection::nextTrackEnergyDistortion(float this_energy=0.) {
  if (_track_energy_resolution < 1e-4)
    return 0.;
  if (_track_energy_by_percent) {
    return _track_energy_distribution( _gen ) * this_energy; 
  } 
  else {
    return _track_energy_distribution( _gen );
  }
}

float TSSelection::nextShowerEnergyDistortion(float this_energy=0.) {
  if (_shower_energy_resolution < 1e-4)
    return 0.;
  if (_shower_energy_by_percent) {
    return _shower_energy_distribution( _gen) * this_energy;
  }
  else {
    return _shower_energy_distribution( _gen );
  }
}

float TSSelection::nextTrackAngleDistortion(float this_angle=0.) {
  if (_track_angle_resolution < 1e-4) 
    return 0.;
  if (_track_angle_by_percent) {
    return _track_angle_distribution(_gen) * this_angle;
  }
  else {
    return _track_angle_distribution(_gen);
  }
}

int TSSelection::nextParticleID(float energy, int true_pdgid) {
  if (!_particle_misid.is_set()) {
    return true_pdgid;
  }
  else {
    return _particle_misid.get(energy)->particle_id(true_pdgid, _random(_gen));
  }

}

float TSSelection::nextShowerAngleDistortion(float this_angle=0.) {
  if (_shower_angle_resolution < 1e-4) 
    return 0.;
  if (_shower_angle_by_percent) {
    return _shower_angle_distribution(_gen) * this_angle;
  }
  else {
    return _shower_angle_distribution(_gen);
  }
}

int TSSelection::get_np(std::vector<PIDParticle>& p) {
  int np = 0;
  for (size_t i=0; i<p.size(); i++) {
    if (p[i].pdg == 2212) {
      np++;
    }
  }
  return np;
}

int TSSelection::get_ntrk(std::vector<PIDParticle>& p) {
  int n_trk = 0;
  for (size_t i=0; i<p.size(); i++) {
    // test for charged pions (the only other possible track?)
    if (p[i].pdg == 2212 || p[i].pdg == 211) {
      n_trk++;
    }
  }
  return n_trk;
}

int TSSelection::get_nl(std::vector<PIDParticle>& p, int lpdg) { 
  int nl = 0;
  for (size_t i=0; i<p.size(); i++) {
    if (p[i].pdg == lpdg) {
      nl++;
    }
  }
  return nl;
}

bool TSSelection::pass_selection(std::vector<PIDParticle>& p, int lpdg, EventType t) {
  // Count protons and the chosen lepton type
  size_t np = (size_t)get_np(p);
  size_t nl = (size_t)get_nl(p, lpdg);
  size_t n_trk = (size_t)get_ntrk(p);

  bool pass_1l1p = nl == 1 && np == 1 && nl + np == p.size();
  bool pass_1lnp = nl == 1 && np >= 1 && nl + np == p.size();
  bool pass_1lntrk = nl == 1 && n_trk >= 1 && nl + n_trk == p.size();
  bool pass_1l0p = nl == 1 && np == 0 && nl + np == p.size();

  bool ret;
  switch(t) {
  case ANY: 
    ret = (pass_1l1p && _accept_1p) 
       || (pass_1lnp && _accept_np) 
       || (pass_1lntrk && _accept_ntrk)
       || (pass_1l0p && _accept_0p);
    break;

  case P0: ret = (pass_1l0p && _accept_0p); break;
  case PN: ret = (pass_1lnp && _accept_np); break;
  case P1: ret = (pass_1l1p && _accept_1p); break;
  case TRKN: ret = (pass_1lntrk && _accept_ntrk); break; 
  }
  return ret;
}


bool TSSelection::initialize(std::vector<std::string> input_files) {
  // initialize input files
  _input_files = input_files;

  // Initialize dataset identifier
  _dataset_id = -1;

  // set producer to null
  _ew_producer = "";

  // initialize resolutions to 0 (perfect resolution)
  _shower_energy_resolution = 0.;
  _track_energy_resolution = 0.;
  _shower_energy_by_percent = false;
  _track_energy_by_percent = false;
  _shower_energy_distribution = std::normal_distribution<float>(0.0, 0.0);
  _track_energy_distribution = std::normal_distribution<float>(0.0, 0.0);

  _shower_angle_resolution = 0.;
  _track_angle_resolution = 0.;
  _shower_angle_by_percent = false;
  _track_angle_by_percent = false;
  _shower_angle_distribution = std::normal_distribution<float>(0.0, 0.0);
  _track_angle_distribution = std::normal_distribution<float>(0.0, 0.0);

  _random = std::uniform_real_distribution<float>(0., 1.);
  _particle_misid = EnergyMap<PDGConfusionMatrix>();

  _accept_1p = true;
  _accept_ntrk = true;
  _accept_np = true;
  _accept_0p = true;

  // setting up random # stuff
  std::random_device rd;
  _gen = std::mt19937( rd() );;

  // Set up the output trees
  assert(_fout);
  assert(_fout->IsOpen());
  _fout->cd();

  _data = new OutputData;
  _tree = new TTree("data", "");
  _tree->Branch("np", &_data->np);
  _tree->Branch("n_trk", &_data->n_trk);
  _tree->Branch("nupdg", &_data->nupdg);
  _tree->Branch("enu", &_data->enu);
  _tree->Branch("q2", &_data->q2);
  _tree->Branch("w", &_data->w);
  _tree->Branch("q0", &_data->q0);
  _tree->Branch("q3", &_data->q3);
  _tree->Branch("int", &_data->int_type);
  _tree->Branch("mode", &_data->int_mode);
  _tree->Branch("ccnc", &_data->ccnc);
  _tree->Branch("eccqe", &_data->eccqe);
  _tree->Branch("eps", &_data->eps);
  _tree->Branch("ppdgs", &_data->ppdgs);
  _tree->Branch("elep", &_data->elep);
  _tree->Branch("thetalep", &_data->thetalep);
  _tree->Branch("philep", &_data->philep);
  _tree->Branch("lpdg", &_data->lpdg);
  _tree->Branch("lpid", &_data->lpid);
  _tree->Branch("llen", &_data->llen);
  _tree->Branch("lexit", &_data->lexit);
  _tree->Branch("bnbweight", &_data->bnbweight);
  _tree->Branch("dataset", &_data->dataset);
  _tree->Branch("weights", &_data->weights);
  _tree->Branch("event_number", &_data->event_number);

  _truth_data = new OutputData;
  _truth_data_tree = new TTree("truth_data", "");
  _truth_data_tree->Branch("np", &_truth_data->np);
  _truth_data_tree->Branch("n_trk", &_truth_data->n_trk);
  _truth_data_tree->Branch("nupdg", &_truth_data->nupdg);
  _truth_data_tree->Branch("enu", &_truth_data->enu);
  _truth_data_tree->Branch("q2", &_truth_data->q2);
  _truth_data_tree->Branch("w", &_truth_data->w);
  _truth_data_tree->Branch("q0", &_truth_data->q0);
  _truth_data_tree->Branch("q3", &_truth_data->q3);
  _truth_data_tree->Branch("int", &_truth_data->int_type);
  _truth_data_tree->Branch("mode", &_truth_data->int_mode);
  _truth_data_tree->Branch("ccnc", &_truth_data->ccnc);
  _truth_data_tree->Branch("eccqe", &_truth_data->eccqe);
  _truth_data_tree->Branch("eps", &_truth_data->eps);
  _truth_data_tree->Branch("ppdgs", &_truth_data->ppdgs);
  _truth_data_tree->Branch("elep", &_truth_data->elep);
  _truth_data_tree->Branch("thetalep", &_truth_data->thetalep);
  _truth_data_tree->Branch("philep", &_truth_data->philep);
  _truth_data_tree->Branch("lpdg", &_truth_data->lpdg);
  _truth_data_tree->Branch("lpid", &_truth_data->lpid);
  _truth_data_tree->Branch("llen", &_truth_data->llen);
  _truth_data_tree->Branch("lexit", &_truth_data->lexit);
  _truth_data_tree->Branch("bnbweight", &_truth_data->bnbweight);
  _truth_data_tree->Branch("dataset", &_truth_data->dataset);
  _truth_data_tree->Branch("weights", &_truth_data->weights);
  _truth_data_tree->Branch("event_number", &_truth_data->event_number);

  _truthtree = new TNtuple("truth", "", "nupdg:enu:ccnc:int:mode:w:q2:lpdg:elep:tlep:npip:npim:npi0:np:nn:fw:ttrk:rtrk:texit:tshr:rshr:sexit");
  _mectree = new TNtuple("mec", "", "nupdg:enu:ccnc:mode:w:q2:lpdg:tlep:ep0:ep1:ep2:ep3:ep4");

  _counts = std::map<EventType, EventCounts>();
  _counts[P0] = EventCounts();
  _counts[P1] = EventCounts();
  _counts[PN] = EventCounts();
  _counts[TRKN] = EventCounts();
  _counts[ANY] = EventCounts();

  _record_truth = true;
  _record_mec = true;

  _event_number = 0;

  return true;
}


bool TSSelection::run() {
  bool ret = true;
  for (gallery::Event ev(_input_files) ; !ev.atEnd(); ev.next()) {
    ret = ret && analyze(&ev);
  }
  return ret;
}

bool TSSelection::analyze(gallery::Event* ev) {
  // Get handles for event data
  art::InputTag gtruth_tag(_mct_producer);
  auto const& gtruth_list = \
    (*ev->getValidHandle<std::vector<simb::GTruth> >(gtruth_tag));

  std::vector<evwgh::MCEventWeight> eventweights_list;
  if (_ew_producer.size() > 0) { 
    art::InputTag eventweight_tag(_ew_producer);
    eventweights_list = \
    (*ev->getValidHandle<std::vector<evwgh::MCEventWeight> >(eventweight_tag));
  }

  art::InputTag mctruth_tag(_mct_producer);
  auto const& mctruth_list = \
    (*ev->getValidHandle<std::vector<simb::MCTruth> >(mctruth_tag));

  art::InputTag mcshower_tag(_mcshw_producer);
  auto const& mcshower_list = \
    (*ev->getValidHandle<std::vector<sim::MCShower> >(mcshower_tag));

  art::InputTag mctrack_tag(_mctrk_producer);
  auto const& mctrack_list = \
    (*ev->getValidHandle<std::vector<sim::MCTrack> >(mctrack_tag));

  // Sanity checks
  assert(_fout);
  assert(mctruth_list.size() == gtruth_list.size());

  // BNB flux weight
  double wbnb = 1.0;
  if (!eventweights_list.empty() &&
      eventweights_list[0].fWeight.find("bnbcorrection_FluxHist") != eventweights_list[0].fWeight.end()) {
    wbnb = eventweights_list[0].fWeight.at("bnbcorrection_FluxHist")[0];
  }

  // Loop through MC truth interactions
  for(size_t i=0; i<mctruth_list.size(); i++) {
    const simb::MCTruth& mctruth = mctruth_list.at(i);
    const simb::GTruth& gtruth = gtruth_list.at(i);

    // first collect the required information, get rid of events that fail 
    // the minimum energy cut
    for (int trial_ind = 0; trial_ind < _n_trials; trial_ind ++) {
      size_t ntracks = 0, nshowers = 0;

      // Keep track of event particle content (currently a little redundant)
      std::vector<PIDParticle> particles_found;
      std::vector<PIDParticle> particles_true;

      // Get vertex-associated contained tracks
      for (size_t j=0; j<mctrack_list.size(); j++) {
        const sim::MCTrack& mct = mctrack_list.at(j);
      
        // First do truth stuff
        
        // Track length
        double s = 0;
        TLorentzVector pos = mct.End().Position();
        for (long k=mct.size()-2; k>=0; k--) {
          s += (pos.Vect() - mct[k].Position().Vect()).Mag();
          pos = mct[k].Position();
        }

        // Apply track cuts for truth information
        if (goodTrack(mct, mctruth)) {
	  // don't apply energy distortion to the "true" particle data
	  particles_true.push_back({
	    mct.PdgCode(),
	    mct.PdgCode(),
	    mct.Start().Momentum(),
	    mct.Start().E() - tsutil::get_pdg_mass(mct.PdgCode()),
	    tsutil::eccqe(mct.Start().Momentum()),
	    s,
	    !tsutil::inFV(mct),
	    mct.TrackID()
	  });
          // truth info on # of tracks
          ntracks++;
        }

        // collect variables for cuts
        int pdg_true = mct.PdgCode();
        float this_angle = mct.Start().Momentum().Theta();
        float this_energy = mct.Start().E() - tsutil::get_pdg_mass(mct.PdgCode());
        float energy_distortion = nextTrackEnergyDistortion( this_energy );
        float angle_distortion = nextTrackAngleDistortion(this_angle);

        bool isEmpty = mct.empty();
        bool isFromNuVertex = tsutil::isFromNuVertex(mctruth, mct);
        bool isPrimaryProcess = mct.Process() == "primary";

        // do pdgid confusion
        int pdg_best = nextParticleID(this_energy + energy_distortion, pdg_true);

        // Un-PID "protons" that are too long or short
        if (pdg_best == 2212 && (s > 80 || s < 12)) {
          pdg_best = -888;
        }

        // TODO: What to do with this???
        // Call all remaining unmatched tracks protons
        //if (pdg_best == -999) {
        //  pdg_best = 2212;
        //}

        // 
        bool pass_cut = tsutil::is_shower_pdgid(pdg_best) ? 
        (goodShower(isFromNuVertex, isPrimaryProcess, this_energy + energy_distortion, pdg_true)) :
        (goodTrack(isEmpty, isFromNuVertex, isPrimaryProcess, this_energy + energy_distortion, pdg_true)); 

        if (pass_cut) {
          auto new_momentum = TLorentzVector(mct.Start().Momentum());
          new_momentum.SetTheta(this_angle + angle_distortion);
          particles_found.push_back({
            pdg_best,
            mct.PdgCode(),
            new_momentum,
            mct.Start().E() - tsutil::get_pdg_mass(mct.PdgCode()) + energy_distortion,
            tsutil::eccqe(mct.Start().Momentum(), energy_distortion, angle_distortion),
            s,
            !tsutil::inFV(mct),
            mct.TrackID()
          });
        }
      } 


      // Get vertex-associated contained showers
      for (size_t j=0; j<mcshower_list.size(); j++) {
        const sim::MCShower& mcs = mcshower_list.at(j);

        // Apply shower cuts
        if (goodShower(mcs, mctruth)) {
	  // don't apply energy distortion to the "true" particle data
	  particles_true.push_back({
	    mcs.PdgCode(),
	    mcs.PdgCode(),
	    mcs.Start().Momentum(),
	    mcs.Start().E() - tsutil::get_pdg_mass(mcs.PdgCode()),
	    tsutil::eccqe(mcs.Start().Momentum()),
	    -1,
	    !tsutil::inFV(mcs),
	    mcs.TrackID()
	  });
	  
	  // keep track of truth # of nshowers
	  nshowers++;
        }

        int pdg_true = mcs.PdgCode();
        float this_energy = mcs.Start().E() - tsutil::get_pdg_mass(mcs.PdgCode());
        float energy_distortion = nextShowerEnergyDistortion( this_energy );
        float this_angle = mcs.Start().Momentum().Theta();
        float angle_distortion = nextShowerAngleDistortion(this_angle);

        bool isEmpty = false; // Is this the best way to set this?
        bool isFromNuVertex = tsutil::isFromNuVertex(mctruth, mcs);
        bool isPrimaryProcess = mcs.Process() == "primary";

        // get pdg from matrix
        int pdg_best = nextParticleID(this_energy + energy_distortion, pdg_true); 

        bool pass_cut = tsutil::is_shower_pdgid(pdg_best) ? 
        (goodShower(isFromNuVertex, isPrimaryProcess, this_energy + energy_distortion, pdg_true)) :
        (goodTrack(isEmpty, isFromNuVertex, isPrimaryProcess, this_energy + energy_distortion, pdg_true)); 

        if (pass_cut) {
	  TLorentzVector new_momentum(mcs.Start().Momentum());
	  new_momentum.SetTheta(this_angle + angle_distortion);
	  particles_found.push_back({
	    pdg_best,
	    mcs.PdgCode(),
	    new_momentum,
	    mcs.Start().E() + energy_distortion - tsutil::get_pdg_mass(mcs.PdgCode()),
	    tsutil::eccqe(mcs.Start().Momentum(), energy_distortion, angle_distortion),
	    -1,
	    !tsutil::inFV(mcs),
	    mcs.TrackID()
	  });
        }
      }

      // Now cut based on the signal regions ntrack/0p/1p/np
      // Each event can only have one lepton and is classified whether the lepton is 
      // a muon or electron

      // get Event counts
      EventType event_types[5] = {PN, P0, P1, TRKN, ANY};
      for (auto t: event_types) {
        bool f_1e = pass_selection(particles_found, 11, t); 
        bool t_1e = pass_selection(particles_true, 11, t); 
        bool f_1m = pass_selection(particles_found, 13, t); 
        bool t_1m = pass_selection(particles_true, 13, t); 
        _counts[t].fill(f_1e, t_1e, f_1m, t_1m);
      }

      // Classify the event (found/true 1l/1m)
      // "True" good_event here means there are one true l and one true p that pass the
      // track/shower cuts (i.e. are in within this specific signal definition).
      bool f_1e = pass_selection(particles_found, 11);
      bool t_1e = pass_selection(particles_true, 11);
      bool f_1m = pass_selection(particles_found, 13);
      bool t_1m = pass_selection(particles_true, 13);

      // Where have all the muons gone?
      //if (t_1m && !f_1m) {
      //  std::cout << "1m1p missed" << std::endl;
      //  std::cout << "true: n=" << particles_true.size() << ": ";
      //  for (size_t i=0; i<particles_true.size(); i++) {
      //    std::cout << particles_true[i].pdg << "/" << particles_true[i].pdgtrue << "(" << particles_true[i].evis << ") ";
      //  }
      //  std::cout << std::endl;
      
      //  std::cout << "found: n=" << particles_found.size() << ": ";
      //  for (size_t i=0; i<particles_found.size(); i++) {
      //    std::cout << particles_found[i].pdg << "/" << particles_found[i].pdgtrue << "(" << particles_found[i].evis << ") ";
      //  }
      //  std::cout << std::endl;
      //}
     
      // Print out PID information mis-IDs
      if ((f_1e && !t_1e) || (f_1m && !t_1m)) {
        std::cout << "true: " << mctruth.GetNeutrino().Nu().E() * 1000
                  << "[" << mctruth.GetNeutrino().InteractionType() << "] ";
        for (size_t k=0; k<particles_true.size(); k++) {
          std::cout << particles_true[k] << " ";
        }
        std::cout << std::endl;
        std::cout << "est: " << ntracks << " tracks, "
                             << nshowers << " showers; ";
        for (size_t k=0; k<particles_found.size(); k++) {
          std::cout << particles_found[k] << " ";
        }
        std::cout << std::endl;
      }
      unsigned int lep_id = 0;
      // Write event to output tree for found 1l1p events
      if (f_1e || f_1m) {
        double eccqe=-1, elep=-1, thetalep=-1, philep=-1, lpdg=-1, lpid=-1, llen=-1, lexit=-1;

        std::vector<double> eps;
        std::vector<int> ppdgs;
        for (size_t k=0; k<particles_found.size(); k++) {
          if (particles_found[k].pdg == 2212) {
            eps.push_back( particles_found[k].evis );
            ppdgs.push_back( particles_found[k].pdgtrue );
          }
          else {
            lep_id = particles_found[k].id;
            eccqe = particles_found[k].eccqe;
            elep = particles_found[k].evis;
            thetalep = particles_found[k].p.Theta();
            philep = particles_found[k].p.Phi();
            lpid = particles_found[k].pdg;
            lpdg = particles_found[k].pdgtrue;
            llen = particles_found[k].len;
            lexit = particles_found[k].exiting;
          }
        }

        const simb::MCNeutrino& nu = mctruth.GetNeutrino();
        const simb::MCParticle& pnu = nu.Nu();
        const simb::MCParticle& plep = nu.Lepton();
        TLorentzVector xp = (pnu.Momentum() - plep.Momentum());

        std::map<std::string, std::vector<double> > wgh;
        if (!eventweights_list.empty()) {
          wgh = eventweights_list[0].fWeight;
        }

        _data->np = get_np(particles_found);
        _data->n_trk = get_ntrk(particles_found);
	_data->nupdg = nu.Nu().PdgCode();
	_data->enu = nu.Nu().E();
	_data->q2 = nu.QSqr();
	_data->w = nu.W();
	_data->q0 = xp.E();
	_data->q3 = xp.Vect().Mag();
	_data->int_type = nu.InteractionType();
	_data->int_mode = nu.Mode();
	_data->ccnc = nu.CCNC();
	_data->eccqe = eccqe;
	_data->eps = eps;
	_data->ppdgs = ppdgs;
	_data->elep = elep;
	_data->thetalep = thetalep;
	_data->philep = philep;
	_data->lpdg = lpdg;
	_data->lpid = lpid;
	_data->llen = llen;
	_data->lexit = lexit;
	_data->bnbweight = wbnb;
	_data->dataset = _dataset_id;
	_data->weights = &wgh;
        _data->event_number = _event_number;
	_tree->Fill();
      }
      // block to contain filling truth data tree -- may include if statement later
      {
        double eccqe=-1, elep=-1, thetalep=-1, philep=-1, lpdg=-1, lpid=-1, llen=-1, lexit=-1;
        std::vector<double> eps;
        std::vector<int> ppdgs;
        for (size_t k=0; k<particles_true.size(); k++) {
          if (particles_true[k].pdg == 2212) {
            eps.push_back( particles_true[k].evis );
            ppdgs.push_back( particles_true[k].pdgtrue );
          }
          else {
            if (lep_id == particles_true[k].id) {
	      eccqe = particles_true[k].eccqe;
	      elep = particles_true[k].evis;
	      thetalep = particles_true[k].p.Theta();
	      philep = particles_true[k].p.Phi();
	      lpid = particles_true[k].pdg;
	      lpdg = particles_true[k].pdgtrue;
	      llen = particles_true[k].len;
	      lexit = particles_true[k].exiting;
            }
          }
        }

        const simb::MCNeutrino& nu = mctruth.GetNeutrino();
        const simb::MCParticle& pnu = nu.Nu();
        const simb::MCParticle& plep = nu.Lepton();
        TLorentzVector xp = (pnu.Momentum() - plep.Momentum());

        std::map<std::string, std::vector<double> > wgh;
        if (!eventweights_list.empty()) {
          wgh = eventweights_list[0].fWeight;
        }

        _truth_data->np = get_np(particles_true);
        _truth_data->n_trk = get_ntrk(particles_true);
	_truth_data->nupdg = nu.Nu().PdgCode();
	_truth_data->enu = nu.Nu().E();
	_truth_data->q2 = nu.QSqr();
	_truth_data->w = nu.W();
	_truth_data->q0 = xp.E();
	_truth_data->q3 = xp.Vect().Mag();
	_truth_data->int_type = nu.InteractionType();
	_truth_data->int_mode = nu.Mode();
	_truth_data->ccnc = nu.CCNC();
	_truth_data->eccqe = eccqe;
	_truth_data->eps = eps;
	_truth_data->ppdgs = ppdgs;
	_truth_data->elep = elep;
	_truth_data->thetalep = thetalep;
	_truth_data->philep = philep;
	_truth_data->lpdg = lpdg;
	_truth_data->lpid = lpid;
	_truth_data->llen = llen;
	_truth_data->lexit = lexit;
	_truth_data->bnbweight = wbnb;
	_truth_data->dataset = _dataset_id;
	_truth_data->weights = &wgh;
        _truth_data->event_number = _event_number;
	_truth_data_tree->Fill();
      }


      // Fill the event truth tree
      if (_record_truth) {
        float vtt[22] = {
	  (float) mctruth.GetNeutrino().Nu().PdgCode(),
	  (float) mctruth.GetNeutrino().Nu().E(),
	  (float) mctruth.GetNeutrino().CCNC(),
	  (float) mctruth.GetNeutrino().InteractionType(),
	  (float) mctruth.GetNeutrino().Mode(),
	  (float) mctruth.GetNeutrino().W(),
	  (float) mctruth.GetNeutrino().QSqr(),
	  (float) mctruth.GetNeutrino().Lepton().PdgCode(),
	  (float) mctruth.GetNeutrino().Lepton().E(),
	  (float) gtruth.fGint,
	  (float) gtruth.fNumPiPlus,
	  (float) gtruth.fNumPiMinus,
	  (float) gtruth.fNumPi0,
	  (float) gtruth.fNumProton,
	  (float) gtruth.fNumNeutron,
	  (float) wbnb,
	  (float) mctrack_list.size(),
	  (float) ntracks,
	  (float) 0,
	  (float) mcshower_list.size(),
	  (float) nshowers,
	  (float) 0
        };
        _truthtree->Fill(vtt);
      }


      if (_record_mec) {
        // Fill tree for MEC events with sorted proton energies
        if (mctruth.GetNeutrino().Mode() == simb::kMEC) {
          std::vector<float> epmec(5, -1);
          for (size_t k=0; k<mctrack_list.size(); k++) {
            const sim::MCTrack& t = mctrack_list[k];
            if (t.PdgCode() == 2212 && t.Process() == "primary") {
              double ke = t.Start().E() - tsutil::get_pdg_mass(t.PdgCode());
              epmec.push_back(ke);
            }
          }

          std::sort(epmec.begin(), epmec.end(), std::greater<>());

          float v2[13] = {
	    (float) mctruth.GetNeutrino().Nu().PdgCode(),
	    (float) mctruth.GetNeutrino().Nu().E(),
	    (float) mctruth.GetNeutrino().CCNC(),
	    (float) mctruth.GetNeutrino().Mode(),
	    (float) mctruth.GetNeutrino().W(),
	    (float) mctruth.GetNeutrino().QSqr(),
	    (float) mctruth.GetNeutrino().Lepton().PdgCode(),
	    (float) mctruth.GetNeutrino().Lepton().E(),
	    (float) epmec[0],
	    (float) epmec[1],
	    (float) epmec[2],
	    (float) epmec[3],
	    (float) epmec[4]
          };
          _mectree->Fill(v2);
        }
      }
    }
  }
  return true;
}


bool TSSelection::finalize() {
  // Print out statistics
  std::cout << "1e true: " << _counts[ANY].true_1e
            << ", good: " << _counts[ANY].good_1e
            << ", miss: " << _counts[ANY].miss_1e
            << std::endl;

  std::cout << "1e eff: "
            << 1.0 * (_counts[ANY].good_1e + _counts[ANY].miss_1e) / _counts[ANY].true_1e
            << std::endl;

  std::cout << "1e pur: "
            << 1.0 * _counts[ANY].good_1e / (_counts[ANY].good_1e + _counts[ANY].miss_1e)
            << std::endl;

  std::cout << "1m true: " << _counts[ANY].true_1m
            << ", good: " << _counts[ANY].good_1m
            << ", miss: " << _counts[ANY].miss_1m << std::endl;

  std::cout << "1m eff: "
            << 1.0 * (_counts[ANY].good_1m + _counts[ANY].miss_1m) / _counts[ANY].true_1m
            << std::endl;

  std::cout << "1m pur: "
            << 1.0 * _counts[ANY].good_1m / (_counts[ANY].good_1m + _counts[ANY].miss_1m)
            << std::endl;

  std::cout << "SHOWER, TRACK ENERGY RESOLUTION: " << _shower_energy_resolution << " "<< _track_energy_resolution << std::endl;
 std::cout << "ACCEPT NP " << _accept_np << std::endl;
 std::cout << "ACCEPT NTRK " << _accept_ntrk << std::endl;

  // record header data
  _fout->cd();
  HeaderData header;
  HeaderData *to_header = &header;
  TTree *header_tree = new TTree("header", "");
  header_tree->Branch("track_producer", &to_header->track_producer);
  header_tree->Branch("fw_producer", &to_header->fw_producer);
  header_tree->Branch("ew_producer", &to_header->ew_producer);
  header_tree->Branch("mct_producer", &to_header->mct_producer);
  header_tree->Branch("mcf_producer", &to_header->mcf_producer);
  header_tree->Branch("mctrk_producer", &to_header->mctrk_producer);
  header_tree->Branch("mcshw_producer", &to_header->mcshw_producer);

  header_tree->Branch("shower_energy_resolution", &to_header->shower_energy_resolution);
  header_tree->Branch("shower_energy_by_percent", &to_header->shower_energy_by_percent);
  header_tree->Branch("track_energy_resolution", &to_header->track_energy_resolution);
  header_tree->Branch("track_energy_by_percent", &to_header->track_energy_by_percent);

  header_tree->Branch("shower_angle_resolution", &to_header->shower_angle_resolution);
  header_tree->Branch("shower_angle_by_percent", &to_header->shower_angle_by_percent);
  header_tree->Branch("track_angle_resolution", &to_header->track_angle_resolution);
  header_tree->Branch("track_angle_by_percent", &to_header->track_angle_by_percent);
 
  header_tree->Branch("n_trials", &to_header->n_trials);
  
  header_tree->Branch("accept_1p", &to_header->accept_1p);
  header_tree->Branch("accept_0p", &to_header->accept_0p);
  header_tree->Branch("accept_np", &to_header->accept_np);
  header_tree->Branch("accept_ntrk", &to_header->accept_ntrk);
  header_tree->Branch("input_files", &to_header->input_files);

  header.track_producer = _track_producer;
  header.fw_producer = _fw_producer;
  header.ew_producer = _ew_producer;
  header.mct_producer = _mct_producer;
  header.mcf_producer = _mcf_producer;
  header.mctrk_producer = _mctrk_producer;
  header.mcshw_producer = _mcshw_producer;

  header.shower_energy_resolution = _shower_energy_resolution;
  header.shower_energy_by_percent = _shower_energy_by_percent;
  header.track_energy_resolution = _track_energy_resolution;
  header.track_energy_by_percent = _track_energy_by_percent;

  header.shower_angle_resolution = _shower_angle_resolution;
  header.shower_angle_by_percent = _shower_angle_by_percent;
  header.track_angle_resolution = _track_angle_resolution;
  header.track_angle_by_percent = _track_angle_by_percent;

  header.n_trials = _n_trials;

  header.accept_1p = _accept_1p;
  header.accept_0p = _accept_0p;
  header.accept_np = _accept_np;
  header.accept_ntrk = _accept_ntrk;

  header.input_files = _input_files;

  header.counts = std::map<EventType, EventCounts>(_counts);

  header_tree->Fill();

  // Write output to ROOT file
  if (_fout) {
    _fout->cd();
    _tree->Write();
    _truthtree->Write();
    _mectree->Write();
    header_tree->Write();
  }

  return true;
}


std::ostream& operator<<(std::ostream& os, const TSSelection::PIDParticle& dt) {
  os << dt.pdg << "(" << dt.evis << ")";
  return os;
}

}
}
