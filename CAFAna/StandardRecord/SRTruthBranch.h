////////////////////////////////////////////////////////////////////////
// \author  Bruno Zamorano
// \date    February 2019
////////////////////////////////////////////////////////////////////////
#ifndef SRTRUTHBRANCH_H
#define SRTRUTHBRANCH_H

#include "StandardRecord/SRNeutrino.h"
#include "StandardRecord/SRParticle.h"

#include <vector>

namespace caf
{
  /// \brief Contains truth information for the slice for the parent
  /// neutrino/cosmic
  class SRTruthBranch
    {
    public:
      SRTruthBranch();
      ~SRTruthBranch();

      SRNeutrino neutrino;
      SRParticle lepton;
      std::vector<SRParticle> finalstate;
    };
  
} // end namespace

#endif // SRTRUTHBRANCH_H
//////////////////////////////////////////////////////////////////////////////
