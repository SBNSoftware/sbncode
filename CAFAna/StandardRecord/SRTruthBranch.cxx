////////////////////////////////////////////////////////////////////////
// \brief   An SRTruthBranch contains vectors of SRTruth.  
//          It is intended for use in the Common Analysis File ROOT trees.
//
// \author  Dominick Rocco
// \date    November 2012
////////////////////////////////////////////////////////////////////////

#include "StandardRecord/SRTruthBranch.h"


namespace caf
{
  
  SRTruthBranch::SRTruthBranch():
  neutrino(),
  lepton(),
  finalstate()
  {  }
  
  SRTruthBranch::~SRTruthBranch() {}
  
  
  void SRTruthBranch::setDefault()
  {
  }
  
  
} // end namespace caf
////////////////////////////////////////////////////////////////////////
