/**
 * \file TSUtil.h
 * \brief Utilities
 * \author A. Mastbaum <mastbaum@uchicago.edu>
 */

#ifndef SBNANA_TSUTIL_H
#define SBNANA_TSUTIL_H

#include <string>
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "TDatabasePDG.h"

class TDatabasePDG;
class TH2F;

namespace ana {
namespace lee_truth_selection {
namespace tsutil {

// A global ROOT PDG table
static TDatabasePDG* gPDGTable = NULL;


// Check if shower endpoint is within the fiducial volume
bool inFV(const sim::MCShower& show);


// Check if track is contained (N.B. just checking endpoint for now)
bool inFV(const sim::MCTrack& track);


// Check if shower is from vertex (within distance tolerance)
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                           float distance=5.0);


// Check if track is from vertex (within distance tolerance)
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                           float distance=5.0);


// Calculate the CCQE energy from lepton four-momentum
double eccqe(const TLorentzVector& v, float energy_distortion=0., float angle_distortion=0.);


// Get the mass for a particle or ion
double get_pdg_mass(const int pdg);

bool is_shower_pdgid(int pdg);

// Build a distribution of track dE/dx vs. residual range
TH2F* HistDEdx(const sim::MCTrack& t, const std::string name="htemp",
               int lowbin=2);

}  // namespace tsutil
}
}

#endif  // SBNANA_TSUTIL_H

