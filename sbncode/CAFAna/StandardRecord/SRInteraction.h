#ifndef SRINTERACTION_H
#define SRINTERACTION_H

#include "StandardRecord/SRNeutrino.h"
#include "StandardRecord/SRFinalStateParticle.h"

namespace caf
{
  class SRWeight_t
  {
  public:
    size_t param_idx;  //!< Parameter name
    size_t universe;  //!< Universe index
    float value;  //!< Parameter value
    float weight;  //!< Weight
  };

  class SRInteraction
  {
  public:
    SRNeutrino neutrino;  //!< The neutrino
    SRFinalStateParticle lepton;  //!< The primary final state lepton
    std::vector<SRFinalStateParticle> finalstate; //!< Final state particles

    std::vector<SRWeight_t> weights;
  };
}

#endif
