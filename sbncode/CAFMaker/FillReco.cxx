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

  void FillTrackVars(const recob::Track& track,
                     caf::SRTrack& srtrack,
                     bool allowEmpty)
  {

    srtrack.npts         = track.CountValidPoints();
    srtrack.len          = track.Length();

  }
  //......................................................................

}
