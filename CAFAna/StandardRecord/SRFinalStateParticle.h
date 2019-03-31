#ifndef SRFINALSTATEPARTICLE_H
#define SRFINALSTATEPARTICLE_H

#include "TVector3.h"

namespace caf
{
  class SRFinalStateParticle
  {
  public:
    int pdg;  //!< PDG Code
    double energy;  //!< Energy
    TVector3 momentum;  //!< Three-momentum
  };
}

#endif
