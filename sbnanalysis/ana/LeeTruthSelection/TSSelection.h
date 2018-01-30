/**
 * \file TSSelection.h
 * \brief A truth-based 1l1p event selection
 * \author A. Mastbaum <mastbaum@uchicago.edu>
 */

#ifndef SBNANA_TSSELECTION_H
#define SBNANA_TSSELECTION_H

#include <map>
#include <string>
#include <random>

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "gallery/Event.h"
#include "TSUtil.h"

class TDatabasePDG;
class TFile;
class TH2F;
class TNtuple;
class TTree;

namespace ana {
namespace lee_truth_selection {
/**
 * \class TSSelection
 * \brief Truth-based selection approximating 1l1p
 */
class TSSelection {
public:
  enum EventType {
    P0, P1, PN, TRKN, ANY
  };

  struct EventCounts {
    size_t true_1e;
    size_t good_1e;
    size_t miss_1e;
    size_t true_1m;
    size_t good_1m;
    size_t miss_1m;


    inline void fill(bool f_1e, bool t_1e, bool f_1m, bool t_1m) {
      if (t_1e) true_1e ++;
      if (f_1e && t_1e) good_1e++;
      if (f_1e && !t_1e) miss_1e++;
      if (t_1m) true_1m++;
      if (f_1m && t_1m) good_1m++;
      if (f_1m && !t_1m) miss_1m++;
    }
  };

  // A structure to hold temporary track/shower data during processing
  struct PIDParticle {
    int pdg;
    int pdgtrue;
    TLorentzVector p;
    double evis;
    double eccqe;
    double len;
    bool exiting;
    unsigned int id;

    // Output stream operator to print a PIDParticle
    friend std::ostream& operator<<(std::ostream& os, const PIDParticle& dt);
  };

  // Class which implements particle-id confusion
  //
  // add(): configure with the rates to id particle "true_pdgid" as particle "test_pdgid"
  // check(): check that every particle in the matrix has a chance of being id's as a particle of 1
  // get(): get the rate that "true_pdgid" is id'd as "test_pdgid"
  // set(): setup a confusion matrix with a previously set up map of particles -> id rates and a set of considered particles
  struct PDGConfusionMatrix {
    // map[true_pdg, fake_pdg] = the % chance to id a particle
    // w/ 'true_pdg' as a the particle 'fake_pdg'
    std::map<std::tuple<int, int>, float> _map;
    std::set<int> _particle_pdgids;
    PDGConfusionMatrix() { 
      _particle_pdgids = std::set<int>();
      _map = std::map<std::tuple<int, int>, float>();
    }

    int particle_id(int true_pdgid, float chance) {
      float accumulator = 0.;
      for (int pdg: _particle_pdgids) {
        accumulator += get(true_pdgid, pdg);
        if (accumulator > chance) {
          return pdg;
        }
      } 
      assert(false);
      return -1;
    }

    void add(int true_pdgid, int test_pdgid, float id_rate) {
      _particle_pdgids.insert(true_pdgid);
      _particle_pdgids.insert(test_pdgid);
      auto key = std::tuple<int, int>(true_pdgid, test_pdgid);
      _map.insert(std::map<std::tuple<int, int>, float>::value_type(key, id_rate));
    }

    float get(int true_pdgid, int test_pdgid) {
      // assert that this particle has been considered
      assert(_particle_pdgids.find(true_pdgid) != _particle_pdgids.end());
      assert(_particle_pdgids.find(test_pdgid) != _particle_pdgids.end());
      
      auto key = std::tuple<int, int>(true_pdgid, test_pdgid);
      if (_map.count(key) > 0) {
        return _map[key];
      }
      else {
        return 0.;
      }
    }

    void set(std::map<std::tuple<int, int>, float> map, std::set<int> vec) {
      _map = map;
      _particle_pdgids = vec;
    }

    void check() {
      for (int true_pdg: _particle_pdgids) {
        float accumulator = 0.;
        for (int test_pdg: _particle_pdgids) {
          accumulator += get(true_pdg, test_pdg);
        }
        assert(abs(accumulator - 1.) < 1e-4);
      } 

    }
  };

  // list of typename T containing an instnace of T for different energy ranges
  //
  // set_energies(): sets the energy ranges being considered -- should be called before calling get()
  // get(): returns a pointer to the object T contained in the given energy bucket
  template<typename T> class EnergyMap {
    public:
    std::vector<float> *_energies;
    std::vector<T> *_objs;
    EnergyMap() {}

    inline void set_energies(std::vector<float> *energies) {
      _energies = energies;
      _objs = new std::vector<T>();
      for (float energy: *_energies) {
       (void) energy;
        T type;
        _objs->push_back(type);
      }
    }

    inline bool is_set() {
      return _energies != NULL;
    }

    T *get(float energy) {
      assert(_energies!= NULL);
      assert(_objs != NULL);
      for (unsigned i = 0; i < _energies->size(); i++) {
        float check_energy = (*_energies)[i];
        if (energy < check_energy) {
          return &(*_objs)[i];
        } 
      }
      assert(false);
      return NULL;
    }
  };

  TSSelection(): _verbose(false) {}

  bool initialize(std::vector<std::string> input_files);

  bool run();

  bool analyze(gallery::Event* ev);

  bool finalize();

  void setVerbose(bool b) { _verbose = b; }

  void setOutputFile(TFile *f) { _fout = f; }

  // Set the producers for data products
  void setFluxWeightProducer(std::string s) { _fw_producer = s; }
  void setEventWeightProducer(std::string s) { _ew_producer = s; }
  void setMCTruthProducer(std::string s) { _mct_producer = s; }
  void setMCFluxProducer(std::string s) { _mcf_producer = s; }
  void setMCShowerProducer(std::string s) { _mcshw_producer = s; }
  void setMCTrackProducer(std::string s) { _mctrk_producer = s; }

  // set the energy resolutions
  // bool is whether the resolution is by_percent or is absolute
  void setShowerEnergyResolution(float, bool);
  void setTrackEnergyResolution(float, bool);
  // get the energy distortion at that given energy generated from internal source of randomness 
  float nextShowerEnergyDistortion(float);
  float nextTrackEnergyDistortion(float);

  // same thing as energy, but for angle
  void setShowerAngleResolution(float, bool);
  void setTrackAngleResolution(float, bool);
  float nextShowerAngleDistortion(float);
  float nextTrackAngleDistortion(float);

  // get the next particle id based on confusion rates set by the PDGConfusionMatrix at that energy
  int nextParticleID(float energy, int true_pdgid);

  // add in particle id rates to the different energy ranges
  void addParticleIDRate(int true_pdgid, int test_pdgid, float rate, float energy)
      { _particle_misid.get(energy)->add(true_pdgid, test_pdgid, rate); }
  void setParticleIDEnergyRange(std::vector<float> range) 
      { _particle_misid.set_energies(new std::vector<float>(range)); }
  void checkParticleIDRates() {
    for (unsigned i = 0; i < _particle_misid._objs->size(); i++) {
      (*_particle_misid._objs)[i].check();
    } 
  }
  
  // setters 
  void setAcceptP(bool, int);
  void setAcceptNTrk(bool b) { _accept_ntrk = b; }
  // number of throws for different random stuff like energy smearing, particle mis-id's
  void setNTrials(int n) { _n_trials = n; }

  // Set a numeric dataset ID, which is written into the tree as a tag
  void setDatasetID(int id) { _dataset_id = id; }

  // Utility function to test if a list of particles is 1lip
  // If EventType is not set, will accept any event enabled by selection.
  // Otherwise, will only accept events inside specified EventType.
  bool pass_selection(std::vector<PIDParticle>& p, int lpdg, EventType t=ANY); 

  // Apply track cuts
  static inline bool goodTrack(const sim::MCTrack& t, const simb::MCTruth& truth, float energy_distortion=0., float angle_distortion=0.) {
    return (!t.empty() &&
            tsutil::isFromNuVertex(truth, t) &&
            t.Process() == "primary" &&
            pass_min_energy_cut(t.PdgCode(), t.Start().E() - tsutil::get_pdg_mass(t.PdgCode()) + energy_distortion));
  }

  static inline bool goodTrack(bool isEmpty, bool isFromNuVertex, bool isPrimaryProcess, int pdgid, float energy) {
    return !isEmpty && isFromNuVertex && isPrimaryProcess && pass_min_energy_cut(pdgid, energy);
  }

  // energy cuts which are stand in for checking if given particle hits 3 wires
  // energy cut definitions are from talk by Wes on Nov. 13 2017
  static inline bool pass_min_energy_cut(int pdgid, float energy) {
    if (abs(pdgid) == 13 || abs(pdgid) == 11 || abs(pdgid) == 211 || pdgid == 22) return energy > 20.;
    if (pdgid == 2212) return energy > 40.;
    if (pdgid == 111) return true;
    // we shouldn't be looking at any other particles
    assert(false);
    return false;
  }

  // Apply shower cuts
  static inline bool goodShower(const sim::MCShower& s, const simb::MCTruth& truth, float energy_distortion=0., float angle_distortion=0.) {
    return (tsutil::isFromNuVertex(truth, s) &&
            s.Process() == "primary" &&
            pass_min_energy_cut(s.PdgCode(), s.Start().E() - tsutil::get_pdg_mass(s.PdgCode()) + energy_distortion));
  }

  static inline bool goodShower(bool isFromNuVertex, bool isPrimaryProcess, int pdgid, float energy) {
    return isFromNuVertex && isPrimaryProcess && pass_min_energy_cut(pdgid, energy);
  }

  // get the event type given Particle info defined above
  int get_nl(std::vector<PIDParticle>& p, int lpdg);
  int get_ntrk(std::vector<PIDParticle>& p);
  int get_np(std::vector<PIDParticle>& p);

  // A structure used to hold TTree output
  struct OutputData {
    OutputData() {
      weights = NULL;
    }
    int np;
    int n_trk;
    int nupdg;
    double enu;
    double q2;
    double w;
    double q0;
    double q3;
    int int_type;
    int int_mode;
    bool ccnc;
    double eccqe;
    std::vector<double> eps;
    std::vector<int> ppdgs;
    double elep;
    double thetalep;
    double philep;
    int lpdg;
    int lpid;
    double llen;
    bool lexit;
    double bnbweight;
    int dataset;
    std::map<std::string, std::vector<double> >* weights;
    unsigned event_number;
  };


  // structure to hold bokkeeping data
  struct HeaderData {
    // Data product producers
    std::string track_producer;
    std::string fw_producer;
    std::string ew_producer;
    std::string mct_producer;
    std::string mcf_producer;
    std::string mctrk_producer;
    std::string mcshw_producer;
    
    // Optionally set some energy resolution
    float shower_energy_resolution;
    bool shower_energy_by_percent;
    float track_energy_resolution;
    bool track_energy_by_percent;
    // and angle distortion
    float shower_angle_resolution;
    bool shower_angle_by_percent;
    float track_angle_resolution;
    bool track_angle_by_percent;
    
    // numbers of things
    int n_trials;
 
    // turn on/off different types of selections
    bool accept_0p;
    bool accept_1p;
    bool accept_np;
    bool accept_ntrk;
    
    // input files
    std::vector<std::string> input_files;  
    // event counts
    std::map<EventType, EventCounts> counts;
  };


protected:
  // Data product producers
  std::string _track_producer;
  std::string _fw_producer;
  std::string _ew_producer;
  std::string _mct_producer;
  std::string _mcf_producer;
  std::string _mctrk_producer;
  std::string _mcshw_producer;

  // Counters for efficiency and purity calculations
  std::map<EventType, EventCounts> _counts;

  // Optionally set some energy resolution
  float _shower_energy_resolution;
  bool _shower_energy_by_percent;
  std::normal_distribution<float> _shower_energy_distribution;
  float _track_energy_resolution;
  bool _track_energy_by_percent;
  std::normal_distribution<float> _track_energy_distribution;
  // and angular resolution
  float _shower_angle_resolution;
  bool _shower_angle_by_percent;
  std::normal_distribution<float> _shower_angle_distribution;
  float _track_angle_resolution;
  bool _track_angle_by_percent;
  std::normal_distribution<float> _track_angle_distribution;

  // source of randomness for confusion stuff
  std::uniform_real_distribution<float> _random;

  // track-shower confusion
  //EnergyMap<float> _track_shower_confusion;

  // particle mis-id's
  EnergyMap<PDGConfusionMatrix> _particle_misid;

  // random stuff
  std::mt19937 _gen;

  // number of times random stuff happens per event
  int _n_trials;

  // turn on/off different types of selections
  bool _accept_1p;
  bool _accept_0p;
  bool _accept_np;
  bool _accept_ntrk;

  // keep track of event index/number
  unsigned _event_number;

  bool _verbose;  //!< Print verbose output
  int _dataset_id;  //!< An arbitrary numeric ID
  OutputData* _data;  //!< Output data
  TTree* _tree;  //!< Output tree

  OutputData* _truth_data;  //!< Output data of truth
  TTree* _truth_data_tree;  //!< Output tree of truth

  TNtuple* _truthtree;
  TNtuple* _mectree;

  // input files
  std::vector<std::string> _input_files;

  // output file
  TFile* _fout;
  // whether to record truth level data
  bool _record_truth;
  bool _record_mec;
};

}
}
#endif  // SBNANA_TSSELECTION_H

