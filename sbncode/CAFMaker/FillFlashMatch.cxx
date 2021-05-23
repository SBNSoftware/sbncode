//////////////////////////////////////////////////////////////////////
// \file    FillFlashMatch.cxx
// \brief   Fill flash match info from multiple algs
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "FillFlashMatch.h"

namespace caf
{


  //......................................................................
  void FillSliceFlashMatch(const sbn::SimpleFlashMatch* fmatch /* can be nullptr */,
                           caf::SRSlice& srslice,
                           bool allowEmpty)
  {
    if (fmatch == nullptr) {
      srslice.fmatch.present = false;
      return;
    }
    srslice.fmatch.present     = fmatch->present;
    srslice.fmatch.time        = fmatch->time;
    srslice.fmatch.charge.q    = fmatch->charge.q;
    srslice.fmatch.light.pe    = fmatch->light.pe;
    srslice.fmatch.score.total = fmatch->score.total;
    srslice.fmatch.score.y     = fmatch->score.y;
    srslice.fmatch.score.z     = fmatch->score.z;
    srslice.fmatch.score.rr    = fmatch->score.rr;
    srslice.fmatch.score.ratio = fmatch->score.ratio;
    srslice.fmatch.charge.center = fmatch->charge.center;
    srslice.fmatch.light.center  = fmatch->light.center;
  }

  //......................................................................
  void FillSliceFlashMatchA(const sbn::SimpleFlashMatch* fmatch /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatch == nullptr) {
      srslice.fmatch_a.present = false;
      return;
    }
    srslice.fmatch_a.present     = fmatch->present;
    srslice.fmatch_a.time        = fmatch->time;
    srslice.fmatch_a.charge.q    = fmatch->charge.q;
    srslice.fmatch_a.light.pe    = fmatch->light.pe;
    srslice.fmatch_a.score.total = fmatch->score.total;
    srslice.fmatch_a.score.y     = fmatch->score.y;
    srslice.fmatch_a.score.z     = fmatch->score.z;
    srslice.fmatch_a.score.rr    = fmatch->score.rr;
    srslice.fmatch_a.score.ratio = fmatch->score.ratio;
    srslice.fmatch_a.charge.center = fmatch->charge.center;
    srslice.fmatch_a.light.center  = fmatch->light.center;
  }

  //......................................................................
  void FillSliceFlashMatchB(const sbn::SimpleFlashMatch* fmatch /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatch == nullptr) {
      srslice.fmatch_b.present = false;
      return;
    }
    srslice.fmatch_b.present     = fmatch->present;
    srslice.fmatch_b.time        = fmatch->time;
    srslice.fmatch_b.charge.q    = fmatch->charge.q;
    srslice.fmatch_b.light.pe    = fmatch->light.pe;
    srslice.fmatch_b.score.total = fmatch->score.total;
    srslice.fmatch_b.score.y     = fmatch->score.y;
    srslice.fmatch_b.score.z     = fmatch->score.z;
    srslice.fmatch_b.score.rr    = fmatch->score.rr;
    srslice.fmatch_b.score.ratio = fmatch->score.ratio;
    srslice.fmatch_b.charge.center = fmatch->charge.center;
    srslice.fmatch_b.light.center  = fmatch->light.center;
  }

  //......................................................................

} // end namespace
