////////////////////////////////////////////////////////////////////////
// \file    SRTrackTruth.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRACKTRUTH_H
#define SRTRACKTRUTH_H

#include <vector>

namespace caf
{
  /// Match from a reconstructed track to a true particle  */
  class ParticleMatch {
  public:
    ParticleMatch(); //!< Constructor
    int G4ID; //!< ID of the particle match, taken from G4 */
    float energy; //!< Total energy matching between reco track and true particle across all three planes [GeV]. NOTE: this energy is a sum of the depoisted energy as seen individually by each plane */
  };

  /// Truth matching information between a SRTrack and a list of SRTrueParticle objects
  /// matching is done using the sum of energies on each plane, so energy variables use the sum of energy
  /// seen by each plane
  //
  // TODO: should this be done differently? Use the best plane?
  class SRTrackTruth {
  public:
    SRTrackTruth(); //!< Constructor
    float total_deposited_energy; //!< True total deposited energy associated with this Track across all 3 planes [GeV]. NOTE: this energy is a sum of the depoisted energy as seen individually by each plane
  
    int                       nmatches; //!< Number of matches
    std::vector<ParticleMatch> matches; //!< List of particle matches, sorted by most energy matched */
    ParticleMatch bestmatch; //!< Best match for track (most energy). Same as index-0 of matches. Useful for columnar analyses.
  };
} // end namespace

#endif // SRTRACKTRUTH_H
//////////////////////////////////////////////////////////////////////////////
