#ifndef __sbnanalysis_ana_LEETruthSelection_MisID__
#define __sbnanalysis_ana_LEETruthSelection_MisID__

/**
 * \file MisID.h
 * \brief Particle mis-ID utilities for the LEE truth selection
 *
 * Author: A. Mastbaum, G. Putnam
 */

#include <algorithm>
#include "Util.h"

namespace ana {
  namespace LEETruthSelection {


/**
 * \class PDGConfusionMatrix
 * \brief Implements particle ID confusion
 */
class PDGConfusionMatrix {
public:
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

  /** Configure with the rates to id particle "true_pdgid" as particle "test_pdgid" */
  void add(int true_pdgid, int test_pdgid, float id_rate) {
    _particle_pdgids.insert(true_pdgid);
    _particle_pdgids.insert(test_pdgid);
    auto key = std::tuple<int, int>(true_pdgid, test_pdgid);
    _map.insert(std::map<std::tuple<int, int>, float>::value_type(key, id_rate));
  }

  /** Get the rate that "true_pdgid" is id'd as "test_pdgid" */
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

  /**
   * Setup a confusion matrix with a previously set up map of particles -> id
   * rates and a set of considered particles.
   */
  void set(std::map<std::tuple<int, int>, float> map, std::set<int> vec) {
    _map = map;
    _particle_pdgids = vec;
  }

  /** Check that every particle in the matrix has a chance of being id's as a particle of 1 */
  void check() {
    for (int true_pdg: _particle_pdgids) {
      float accumulator = 0.;
      for (int test_pdg: _particle_pdgids) {
        accumulator += get(true_pdg, test_pdg);
      }
      assert(abs(accumulator - 1.) < 1e-4);
    }
  }

protected:
  std::map<std::tuple<int, int>, float> _map;  //!< Percent chance to id a particle [true, fake]
  std::set<int> _particle_pdgids;
};


/**
 * \class EnergyMap
 *
 * List of typename T containing an instance of T for different energy ranges
 */
template<typename T>
class EnergyMap {
public:
  EnergyMap() {}

  /** Set the energy ranges being considered -- should be called before calling get() */
  inline void set_energies(std::vector<float> *energies) {
    _energies = energies;
    _objs = new std::vector<T>(energies.size());
  }

  inline bool is_set() {
    return _energies != NULL;
  }

  /** Returns a pointer to the object T contained in the given energy bucket */
  T* get(float energy) {
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

  std::vector<float> *_energies;
  std::vector<T> *_objs;
};

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
  
  }  // namespace LEETruthSelection
}  // namespace ana

#endif  // __sbnanalysis_ana_LEETruthSelection_MisID__

