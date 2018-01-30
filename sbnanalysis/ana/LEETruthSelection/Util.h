#ifndef __sbnanalysis_ana_LEETruthSelection_Utils__
#define __sbnanalysis_ana_LEETruthSelection_Utils__

/**
 * \file Util.h
 * \brief Utilities for the LEE truth selection
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

#include <string>
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "TDatabasePDG.h"

class TDatabasePDG;
class TH2F;

namespace ana {
  namespace lee_truth_selection {

/**
 * \struct PIDParticle
 * \brief A structure to hold temporary track/shower data during processing
 */
struct PIDParticle {
  int pdg;
  int pdgtrue;
  TLorentzVector p;
  double evis;
  double eccqe;
  double len;
  bool exiting;
  unsigned int id;

  /** Output stream operator to print a PIDParticle. */
  friend std::ostream& operator<<(std::ostream& os, const PIDParticle& dt);
};

/** Collection of utilities for the LEE truth selection */
namespace util {

/** A global ROOT PDG table. */
static TDatabasePDG* gPDGTable = NULL;



/**
 * Check if shower endpoint is within the fiducial volume.
 *
 * \param show The MCShower to test
 * \returns True if endpoint is within FV
 */
bool InFV(const sim::MCShower& show);


/**
 * Check if track is contained (N.B. just checking endpoint for now).
 *
 * \param show The MCTrack to test
 * \returns True if endpoint is within FV
 */
bool InFV(const sim::MCTrack& track);


/**
 * Check if shower is from vertex (within distance tolerance).
 *
 * \param mc The MCTruth containing the vertex
 * \param show The shower to test
 * \param distance The tolerance in cm
 * \returns True is shower start is within tolerance of vertex
 */
bool IsFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                    float distance=5.0);


/**
 * Check if track is from vertex (within distance tolerance).
 *
 * \param mc The MCTruth containing the vertex
 * \param track The track to test
 * \param distance The tolerance in cm
 * \returns True is shower start is within tolerance of vertex
 */
bool IsFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                    float distance=5.0);


/**
 * Calculate the CCQE energy from lepton four-momentum.
 *
 * \param v The lepton four-momentum
 * \param energy_distortion Amount to randomly smear the energy
 * \returns Energy using the CCQE calculation
 */
double ECCQE(const TLorentzVector& v, float energy_distortion=0.0);


/**
 * Get the mass for a particle or ion.
 *
 * \param pdg The integer PDG code
 * \returns The particle mass
 */
double GetPDGMass(const int pdg);


}  // namespace util
  }  // namespace lee_truth_selection
}  // namespace ana

#endif  // __sbnanalysis_ana_LEETruthSelection_Util__

