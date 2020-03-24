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
    // nshw_em(-1),
    // nshw_pandora(-1)
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
    // nshw_em = shw_em.size();
    // nshw_pandora = shw_pandora.size();

  }
 
} // end namespace caf
////////////////////////////////////////////////////////////////////////
