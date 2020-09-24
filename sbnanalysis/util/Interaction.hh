#ifndef __sbnanalysis_util_Interaction__
#define __sbnanalysis_util_Interaction__

#include <TVector3.h>
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

namespace util {

/**
 * Calculate CCQE energy from associated lepton information.
 *
 * \param l_momentum Lepton momentum (any units, used only to get angle info)
 * \param l_energy Lepton energy in GeV
 * \return CCQE energy in GeV.
 */
double ECCQE(const TVector3& l_momentum, double l_energy);

/**
 *
 * Calculate length between two points contained in a list of
 * (non-overlapping) volumes. Units are abitrary, but must be consistent
 * between the three inputs.
 *
 * \param v0 The first point
 * \param v1 The second point
 * \param boxes A list of non-overlapping boxes through which the
 * contained length will be calculated.
 * \return The contained length in units of the inputs
 */
double ContainedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geo::BoxBoundedGeo> &boxes);

/**
 * Calculated length of an MCParticle trajectory through a list of
 * (non-overlapping) volumes.
 *
 * \param particle The particle to calculate the length for
 * \param active_volumes A list of non-overlapping volumes through which
 * the contained length will be calculated
 * \return The length of the MCParticle trajectory contained in the
 * active_volumes. Returns 0 if the list of volumes is empty or the
 * particle does not contain trajectory points.
 */
double MCParticleContainedLength(const simb::MCParticle &particle, const std::vector<geo::BoxBoundedGeo> &active_volumes);

/**
 * Calculate length of an MCParticle trajectory
 *
 * \param particle the MCParticle
 * \return The length travelled by the particle. Returns 0 if the
 * particle does not contain trajectory points. 
 */
double MCParticleLength(const simb::MCParticle &particle);

}  // namespace util

#endif  // __sbnanalysis_util_Interaction__

