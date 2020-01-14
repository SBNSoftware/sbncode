//////////////////////////////////////////////////////////////////////
// \file    FillReco.cxx
// \brief   Fill reco SR branches 
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "FillReco.h"

namespace caf
{

  //......................................................................
  void FillSliceVars(const recob::Slice& slice,
                     caf::SRSlice& srslice,
                     bool allowEmpty)
  {

    srslice.id           = slice.ID();
    srslice.charge       = slice.Charge();

  }
  //......................................................................

}
