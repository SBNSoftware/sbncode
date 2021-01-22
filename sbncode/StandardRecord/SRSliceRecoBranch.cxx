////////////////////////////////////////////////////////////////////////
// \file    SRSliceRecoBranch.cxx
////////////////////////////////////////////////////////////////////////

#include "sbncode/StandardRecord/SRSliceRecoBranch.h"


namespace caf 
{
  
  SRSliceRecoBranch::SRSliceRecoBranch() :
    ntrk(0),
    ntrk_split(0),
    //    nvtx(0),
    nshw(0)
    // nshw_em(0),
    // nshw_pandora(0)
  {  
  }
  
  SRSliceRecoBranch::~SRSliceRecoBranch()
  { 
  }

  void SRSliceRecoBranch::fillSizes()
  {
    
    //    nvtx   = vtx.size();
    ntrk   = trk.size();
    ntrk_split = trk_split.size();
    nshw   = shw.size();
    // nshw_em = shw_em.size();
    // nshw_pandora = shw_pandora.size();

  }
 
} // end namespace caf
////////////////////////////////////////////////////////////////////////
