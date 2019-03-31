#ifndef SRINTERACTION_H
#define SRINTERACTION_H

#include "StandardRecord/SRNeutrino.h"
#include "StandardRecord/SRFinalStateParticle.h"

namespace caf
{
  class SRInteraction
  {
  public:
    SRNeutrino neutrino;  //!< The neutrino
    SRFinalStateParticle lepton;  //!< The primary final state lepton
    std::vector<SRFinalStateParticle> finalstate; //!< Final state particles

    // std::map<std::string, std::vector<double> > weights; // TODO TODO TODO
  };
}

#endif
