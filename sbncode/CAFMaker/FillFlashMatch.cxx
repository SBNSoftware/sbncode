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
    srslice.fmatch.present    = fmatch->present;
    srslice.fmatch.time       = fmatch->time;
    srslice.fmatch.chargeQ    = fmatch->charge.q;
    srslice.fmatch.lightPE    = fmatch->light.pe;
    srslice.fmatch.score      = fmatch->score.total;
    srslice.fmatch.scoreY     = fmatch->score.y;
    srslice.fmatch.scoreZ     = fmatch->score.z;
    srslice.fmatch.scoreRR    = fmatch->score.rr;
    srslice.fmatch.scoreRatio = fmatch->score.ratio;
    srslice.fmatch.scoreSlope = fmatch->score.slope;
    srslice.fmatch.scorePEToQ = fmatch->score.petoq;
    srslice.fmatch.chargeCenter = fmatch->charge.center;
    srslice.fmatch.chargeWidth = fmatch->charge.width;
    srslice.fmatch.lightCenter  = fmatch->light.center;
    srslice.fmatch.lightWidth  = fmatch->light.width;
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
    srslice.fmatch_a.present    = fmatch->present;
    srslice.fmatch_a.time       = fmatch->time;
    srslice.fmatch_a.chargeQ    = fmatch->charge.q;
    srslice.fmatch_a.lightPE    = fmatch->light.pe;
    srslice.fmatch_a.score      = fmatch->score.total;
    srslice.fmatch_a.scoreY     = fmatch->score.y;
    srslice.fmatch_a.scoreZ     = fmatch->score.z;
    srslice.fmatch_a.scoreRR    = fmatch->score.rr;
    srslice.fmatch_a.scoreRatio = fmatch->score.ratio;
    srslice.fmatch_a.scoreSlope = fmatch->score.slope;
    srslice.fmatch_a.scorePEToQ = fmatch->score.petoq;
    srslice.fmatch_a.chargeCenter = fmatch->charge.center;
    srslice.fmatch_a.chargeWidth  = fmatch->charge.width;
    srslice.fmatch_a.lightCenter  = fmatch->light.center;
    srslice.fmatch_a.lightWidth   = fmatch->light.width;
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
    srslice.fmatch_b.present    = fmatch->present;
    srslice.fmatch_b.time       = fmatch->time;
    srslice.fmatch_b.chargeQ    = fmatch->charge.q;
    srslice.fmatch_b.lightPE    = fmatch->light.pe;
    srslice.fmatch_b.score      = fmatch->score.total;
    srslice.fmatch_b.scoreY     = fmatch->score.y;
    srslice.fmatch_b.scoreZ     = fmatch->score.z;
    srslice.fmatch_b.scoreRR    = fmatch->score.rr;
    srslice.fmatch_b.scoreRatio = fmatch->score.ratio;
    srslice.fmatch_b.scoreSlope = fmatch->score.slope;
    srslice.fmatch_b.scorePEToQ = fmatch->score.petoq;
    srslice.fmatch_b.chargeCenter = fmatch->charge.center;
    srslice.fmatch_b.chargeWidth  = fmatch->charge.width;
    srslice.fmatch_b.lightCenter  = fmatch->light.center;
    srslice.fmatch_b.lightWidth   = fmatch->light.width;
  }

  //......................................................................

} // end namespace
