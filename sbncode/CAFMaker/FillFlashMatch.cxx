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
                           caf::SRFlashMatch& srflash,
                           bool allowEmpty)
  {
    if (fmatch == nullptr) {
      srflash.present = false;
      return;
    }
    srflash.present    = fmatch->present;
    srflash.time       = fmatch->time;
    srflash.chargeQ    = fmatch->charge.q;
    srflash.lightPE    = fmatch->light.pe;
    srflash.score      = fmatch->score.total;
    srflash.scoreY     = fmatch->score.y;
    srflash.scoreZ     = fmatch->score.z;
    srflash.scoreRR    = fmatch->score.rr;
    srflash.scoreRatio = fmatch->score.ratio;
    srflash.scoreSlope = fmatch->score.slope;
    srflash.scorePEToQ = fmatch->score.petoq;
    srflash.chargeCenter = fmatch->charge.center;
    srflash.chargeWidth = fmatch->charge.width;
    srflash.lightCenter  = fmatch->light.center;
    srflash.lightWidth  = fmatch->light.width;
  }

  //......................................................................

} // end namespace
