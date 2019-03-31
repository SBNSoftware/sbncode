#ifndef SRRECOINTERACTION_H
#define SRRECOINTERACTION_H

#include "StandardRecord/SRInteraction.h"

namespace caf
{
  class SRRecoInteraction
  {
  public:
    SRInteraction truth; //!< Contains truth level information about interaction
    int truth_index;
    double reco_energy;
    double weight;  
  };
}


#endif
