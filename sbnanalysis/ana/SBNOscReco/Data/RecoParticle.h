#ifndef _sbnumurecodata_RecoParticle_hh
#define _sbnumurecodata_RecoParticle_hh

#include <vector>

#include "larcorealg/Geometry/BoxBoundedGeo.h"


namespace numu {
/** Reconstructed information about each particle. Internal struct used
* to combine information on each reconstructed particle. 
* */
struct RecoParticle {
  bool p_is_clear_cosmic; //!< Taken from Pandora metadata "is_clear_cosmic"
  bool p_is_neutrino; //!< Taken from Pandora metadata "is_neutrino"
  double p_nu_score; //!< Take from Pandora metadata "nu_score"
  std::vector<geo::Point_t> vertices; //!< List of vertices associated with the particle
  std::vector<size_t> daughters; //!< Daughters of the particle in the "particle flow". Value represents index into pandora information.
  size_t self; //!< Index into pandora information of this particle
};
}
#endif
