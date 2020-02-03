#ifndef _sbnumurecodata_TruthMatch_hh
#define _sbnumurecodata_TruthMatch_hh

#include "TVector3.h"

#include "ana/SBNOscReco/Data/Mode.h"

namespace numu {
/**
 * Truth information on a reconstructed Track.
*/
struct TrackTruth {
  float total_deposited_energy;
  struct ParticleMatch {
    int G4ID;
    float energy;
  };
  std::vector<ParticleMatch> matches;

  int GetPrimaryMatchID() const {
    if (matches.size() == 0) return -1;
    // matches are sorted such that first match
    // has highest matched energy
    return matches[0].G4ID;
  }

  float Purity() const {
    if (matches.size() == 0) return 0.;
    return matches[0].energy / total_deposited_energy;
  }

};

/**
 * Truth information on a reconstructed slice.
 */
struct SliceTruth {
  InteractionMode mode;
  int interaction_id;
  float total_deposited_energy;
  float neutrino_match_energy;
  float cosmic_energy;
  float unmatched_energy;
  bool IsPrimary() const {
    return mode != numu::mCCNonPrimary && mode != numu::mNCNonPrimary;
  }
};

}
#endif
