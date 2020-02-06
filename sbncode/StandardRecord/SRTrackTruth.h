////////////////////////////////////////////////////////////////////////
// \file    SRTrackTruth.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRACKTRUTH_H
#define SRTRACKTRUTH_H

#include <vector>

namespace caf
{
  class SRTrackTruth {
  public:
    float total_deposited_energy; //!< True total deposited energy associated with this Track [GeV]
  
    /// Match from a reconstructed track to a true particle  */
    class ParticleMatch {
    public:
      int G4ID; //!< ID of the particle match, taken from G4 */
      float energy; //!< Total energy matching between reco track and true particle [GeV] */
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
