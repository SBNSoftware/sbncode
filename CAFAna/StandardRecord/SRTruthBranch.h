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

      std::vector<SRNeutrino>   neutrino;   ///< implemented as a vector to maintain mc.nu structure, i.e. not a pointer, but with 0 or 1 entries. 
      std::vector<SRParticle>     lepton;
      std::vector<SRParticle> finalstate;
      void setDefault();

    };
  
} // end namespace

#endif // SRTRUTHBRANCH_H
//////////////////////////////////////////////////////////////////////////////
