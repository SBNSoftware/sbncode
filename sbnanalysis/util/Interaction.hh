#ifndef __UTILITY_INTERACTION_HH
#define __UTILITY_INTERACTION_HH

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "gallery/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "core/Event.hh"

namespace util {
/** Calculate CCQE energy from associated lepton information. Energy in GeV. 
 *
 * \param l_momentum Lepton momentum (in any units -- used only to get angle info)
 * \param l_energy Lepton energy in GeV
 *
 * \return CCQE energy in GeV.
 * */
double ECCQE(const TVector3& l_momentum, double l_energy);
}  // namespace util
#endif
