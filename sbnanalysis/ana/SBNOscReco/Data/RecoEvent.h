#ifndef _sbnumurecodata_RecoEvent_hh
#define _sbnumurecodata_RecoEvent_hh

#include <vector>
#include <map>

#include "TVector3.h"

#include "RecoParticle.h"
#include "RecoTrack.h"
#include "Match.h"

namespace numu {

struct RecoSlice {
  int primary_index;
  int primary_track_index;
  std::map<size_t, RecoParticle> particles;
  std::map<size_t, RecoTrack> tracks;
};


/**
*  Reconstruction information for each neutrino vertex.
*  Produced from both reconstruction and truth information
*/
struct RecoInteraction {
  RecoSlice slice; //!< Particle content of the interaction
  TVector3 position; //!< location of the vertex
  double nu_energy; //!< true/reconstructed neutrino energy
  TruthMatch match; //!< Info for mathing to truth
  int multiplicity;
  RecoTrack primary_track;
};


/** Reconstruction Information about Event */
struct RecoEvent {
  std::vector<RecoTrack> tracks;
  std::map<size_t, RecoTrack> true_tracks;
  std::vector<RecoInteraction> reco; //!< List of reconstructed vertices
  std::vector<RecoInteraction> truth; //!< List of truth vertices
};
}
#endif
