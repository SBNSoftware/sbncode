#ifndef _sbncode_TruthMatch_hh
#define _sbncode_TruthMatch_hh

#include <vector>
#include "../Data/RecoEvent.h"

// Do truth match with the NumuReco objects defined in Data/

namespace numu {
  /**
 *  Returns the truth match interaction in a RecoInteraction
 *
 *  \param truth The list of true interactions in the event
 *  \param reco_tracks The set of reconstructed tracks in the event
 *  \param reco The reconstructed interaction in the event. The existing content
 *              of the reco.match instance does not affect the algorithm
 *  \return The TruthMatch object for the RecoInteraction reco  
 */
  TruthMatch InteractionTruthMatch(const std::vector<RecoInteraction> &truth, const std::map<size_t, RecoTrack> &reco_tracks, const RecoInteraction &reco);

  /**
 * Corrects a list of reco interaction objects when some are matched to the
 * same truth neutrino interaction by making some of them non-primary
 *
 * \param event The reconstructed event
 * \param recos The list of reconstructed vertices. The "match" field of each
 *              RecoInteraction may be edited by this function
 *
 */
  void CorrectMultiMatches(const RecoEvent &event, std::vector<RecoInteraction> &recos);
}



#endif
