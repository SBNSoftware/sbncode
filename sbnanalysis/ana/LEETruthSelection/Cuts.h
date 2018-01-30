#ifndef __sbnanalysis_ana_LEETruthSelection_Cuts__
#define __sbnanalysis_ana_LEETruthSelection_Cuts__

/**
 * \file Cuts.h
 * \brief Cut functions for the LEE truth selection
 *
 * Author: A. Mastbaum, G. Putnam
 */

#include <algorithm>
#include "Util.h"

namespace ana {
  namespace lee_truth_selection {

/** A functor to test a PIDParticle's PDG code. */
class PDGTest {
public:
  PDGTest(std::vector<int> _pdg) : pdg(_pdg) {}
  inline bool operator()(PIDParticle& p) {
    return std::find(pdg.begin(), pdg.end(), p.pdg) != pdg.end();
  }

private:
  std::vector<int> pdg;
};


/**
 * Count the number of protons in a list.
 *
 * \param p A list of identified particles
 * \returns The number of protons
 */
int GetNp(std::vector<PIDParticle>& p);


/**
 * Count the number of protons and pi+ in a list.
 *
 * \param p A list of identified particles
 * \returns The number of protons
 */
int GetNtracks(std::vector<PIDParticle>& p);


/**
 * Count the number of leptons in a list.
 *
 * \param p A list of identified particles
 * \param lpdg The lepton PDG code
 * \returns The number of protons
 */
int GetNlep(std::vector<PIDParticle>& p, int lpdg); 


/**
 * Utility function to test if a list of particles meets criteria.
 *
 */
bool TestSelection(std::vector<PIDParticle>& p, int lpdg, EventType t);


/**
 * Apply object cuts.
 *
 * \param t isFromNuVertex Is the track originating at the neutrino vertex?
 * \param isPrimaryProcess Is the track process tagged as "primary"
 * \param pdg Track (true) PDG code
 * \param ke Track start kinetic energy
 * \returns True if track passes all cuts
 */
bool GoodObject(bool isFromNuVertex, bool isPrimaryProcess, int pdg, float ke);


/**
 * Kinetic energy threshold.
 *
 * Depends on the (true) particle type.
 *
 * Energy cut definitions are from talk by W. Ketchum on Nov. 13 2017
 *
 * \param pdg The particle PDG code
 * \param ke Kinetic energy
 * \returns True if KE is over threshold
 */
bool KineticEnergyThreshold(int pdg, float ke);

  }  // namespace lee_truth_selection
}  // namespace ana

#endif  // __sbnanalysis_ana_LEETruthSelection_Cuts__

