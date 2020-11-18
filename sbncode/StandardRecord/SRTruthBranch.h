////////////////////////////////////////////////////////////////////////
// \file    SRTruthBranch.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRUTHBRANCH_H
#define SRTRUTHBRANCH_H

#include "sbncode/StandardRecord/SRTrueInteraction.h"

#include <vector>

namespace caf
{
  /// Vectors of reconstructed vertices found by various algorithms
  class SRTruthBranch
  {
  public:
    SRTruthBranch();
    ~SRTruthBranch();
            
    std::vector<SRTrueInteraction> nu;   ///< Vector of true nu or cosmic
    size_t                        nnu;   ///< Number of true nu or cosmic

    void fillSizes();
      
  };
  
} // end namespace

#endif // SRTRUTHBRANCH_H
////////////////////////////////////////////////////////////////////////////
