#ifndef __sbnanalysis_util_Interaction__
#define __sbnanalysis_util_Interaction__

#include <TVector3.h>

namespace util {

/**
 * Calculate CCQE energy from associated lepton information.
 *
 * \param l_momentum Lepton momentum (any units, used only to get angle info)
 * \param l_energy Lepton energy in GeV
 * \return CCQE energy in GeV.
 */
double ECCQE(const TVector3& l_momentum, double l_energy);

}  // namespace util

#endif  // __sbnanalysis_util_Interaction__

