////////////////////////////////////////////////////////////////////////
// \author  Bruno Zamorano
// \date    February 2019
////////////////////////////////////////////////////////////////////////
#ifndef SRPARTICLE_H
#define SRPARTICLE_H

#include <TVector3.h>

namespace caf
{
  /// The SRParticle is a representation of final state pcle. info
  class SRParticle
    {
    public:
      SRParticle();
      ~SRParticle() {  };

      int          pdg;           ///< PDG code
      double       energy;        ///< True energy [GeV]
      TVector3     momentum;      ///< True momentum [GeV]
    };

} // end namespace

#endif // SRPARTICLE_H
//////////////////////////////////////////////////////////////////////////////
