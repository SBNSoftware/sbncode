#ifndef _sbnumurecodata_RecoEvent_hh
#define _sbnumurecodata_RecoEvent_hh

#include <vector>
#include <map>

#include "TVector3.h"

#include "ana/SBNOscReco/Data/TrueParticle.h"
#include "ana/SBNOscReco/Data/RecoParticle.h"
#include "ana/SBNOscReco/Data/RecoTrack.h"
#include "ana/SBNOscReco/Data/TruthMatch.h"
#include "ana/SBNOscReco/Data/DetInfo.h"
#include "ana/SBNOscReco/Data/MCType.h"
#include "ana/SBNOscReco/Data/FlashMatch.h"

namespace numu {

/**
 *  Slice of TPC charge containing a list of "particles" (Reco objects)
 *  and TPC tracks (showers to be added). Also contains the primary 
 *  track of the slice
 */
struct RecoSlice {
  int primary_index; //!< Index of the primary particle of this slice. 
                     // If this slice is a neutrino slice, then the priamry particle is the neutrino
  int primary_track_index; //!< Index of the primary track
  std::map<size_t, RecoParticle> particles; //!< Map of particle index to particle information
  std::vector<size_t> tracks; //!< List of track indices contained in this slice
  FlashMatch flash_match; //!< Result of flash matching algorithm on this slice
};


/**
*  Reconstruction information for each neutrino vertex.
*  Produced from both reconstruction and truth information
*/
struct RecoInteraction {
  RecoSlice slice; //!< Particle content of the interaction
  TVector3 position; //!< location of the vertex
  float nu_energy; //!< true/reconstructed neutrino energy
  TruthMatch match; //!< Info for mathing to truth
  int multiplicity; //!< Number of tracks in this interaction
  RecoTrack primary_track; //!< Stores the primary track in the slice
};

/** Reconstruction Information about Event */
struct RecoEvent {
  std::map<size_t, RecoTrack> tracks; //!< Map of track indices to Track information.
  std::map<size_t, TrueParticle> particles; //!< Map of indices to True particle information
  std::vector<RecoInteraction> reco; //!< List of reconstructed vertices
  std::vector<CRTHit> in_time_crt_hits; //!< List of crt hits in time with the beam spill
  std::vector<FlashTriggerPrimitive> flash_trigger_primitives; //!< List of trigger primitives from optical detectors
  MCType type;
};
}
#endif
