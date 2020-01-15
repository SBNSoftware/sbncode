////////////////////////////////////////////////////////////////////////
// \file    SRSliceRecoBranch.cxx
////////////////////////////////////////////////////////////////////////

#include "SRSliceRecoBranch.h"


namespace caf 
{
  
  SRSliceRecoBranch::SRSliceRecoBranch() :
    ntrk(-1),
    //    nvtx(-1),
    nshw(-1)
  {  
  }
  
  SRSliceRecoBranch::~SRSliceRecoBranch()
  { 
  }

  void SRSliceRecoBranch::fillSizes()
  {
    
    //    nvtx   = vtx.size();
    ntrk   = trk.size();
    nshw   = shw.size();

  }
 
} // end namespace caf
////////////////////////////////////////////////////////////////////////
