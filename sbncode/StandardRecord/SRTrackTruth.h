////////////////////////////////////////////////////////////////////////
// \file    SRTrackTruth.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRACKTRUTH_H
#define SRTRACKTRUTH_H

#include <vector>

namespace caf
{
  /// Truth matching information between a SRTrack and a list of SRTrueParticle objects
  /// matching is done using the sum of energies on each plane, so energy variables use the sum of energy
  /// seen by each plane
  //
  // TODO: should this be done differently? Use the best plane?
  class SRTrackTruth {
  public:
    float total_deposited_energy; //!< True total deposited energy associated with this Track across all 3 planes [GeV]. NOTE: this energy is a sum of the depoisted energy as seen individually by each plane
  
    /// Match from a reconstructed track to a true particle  */
    class ParticleMatch {
    public:
      int G4ID; //!< ID of the particle match, taken from G4 */
      float energy; //!< Total energy matching between reco track and true particle across all three planes [GeV]. NOTE: this energy is a sum of the depoisted energy as seen individually by each plane */
    };

    std::vector<ParticleMatch> matches; //!< List of particle matches, sorted by most energy matched */

    /// Get the G4ID of the primary (most energy) particle match to this track (returns -1 if no match)
    int GetPrimaryMatchID() const;
    
    /// Get the purity of the energy associated with this track relative to the best-matched particle (0 if no match) 
    float Purity() const;
  };
} // end namespace

#endif // SRTRACKTRUTH_H
//////////////////////////////////////////////////////////////////////////////
