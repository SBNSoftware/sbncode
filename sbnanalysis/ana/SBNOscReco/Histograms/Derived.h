#ifndef _sbncode_Derived_hh
#define _sbncode_Derived_hh

#include <vector>
#include "../Data/RecoEvent.h"
#include "core/Event.hh"

// helper functions to calculate derived quantities from
// output classes of the reconstruction processing

namespace numu {
  /**
 * Distance between one interation vertex and a list of candidate matching vertices
 * \param truth The true neutrino interaction to match
 * \param candidates The list of candidate vertices matching to this one
 *
 * \return The minimum distance between the list of candidates and the vertex
 */
  float dist2Match(const event::Interaction &truth, const std::vector<numu::RecoInteraction> &candidates);

  /**
 * Get the completon of a reconstructed track matching to a true track
 *
 * \param truth_index The index of the true particle to match
 * \param event The reconstructed event
 *
 * \return The maximum completion of all reconstructed tracks to the true track. -1 
 *         if there is no matching reconstructed track
 */
  float trackMatchCompletion(unsigned truth_index, const numu::RecoEvent &event);
}



#endif
