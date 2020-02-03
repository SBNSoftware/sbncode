#ifndef _sbncode_PriamryTrack_hh
#define _sbncode_PriamryTrack_hh

#include <vector>
#include "../Data/RecoEvent.h"

// Do truth match with the NumuReco objects defined in Data/

namespace numu {
  /**
 * Algorithm to select the primary track -- get the longest one
 *
 * \param tracks All reconstructed tracks in the event
 * \param slice The RecoSlice to select the primary track for
 * \return The ID of the primary track. Returns -1 if no such track exists.
 */
 int SelectLongestTrack(const std::map<size_t, RecoTrack> &tracks, const RecoSlice &slice);

 /**
 * Algorithm to select the primary track -- get the longest one that has a better 
 * particle ID for muon than proton
 *
 * \param tracks All reconstructed tracks in the event
 * \param slice The RecoSlice to select the primary track for
 * \return The ID of the primary track. Returns -1 if no such track exists.
 */
 int SelectLongestIDdMuon(const std::map<size_t, RecoTrack> &tracks, const RecoSlice &slice);

 int SelectMuon(const std::map<size_t, RecoTrack> &tracks, const RecoSlice &slice);

}



#endif
