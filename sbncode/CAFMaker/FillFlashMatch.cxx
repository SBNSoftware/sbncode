//////////////////////////////////////////////////////////////////////
// \file    FillFlashMatch.cxx
// \brief   Fill flash match info from multiple algs
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "FillFlashMatch.h"

namespace caf
{


  //......................................................................
  void FillSliceFlashMatchA(const sbn::SimpleFlashMatch* fmatch /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatch == nullptr) {
      srslice.fmatch.present = false;
      return;
    }
    srslice.fmatch.present   = fmatch->mPresent;
    srslice.fmatch.time      = fmatch->mTime;
    srslice.fmatch.charge_q  = fmatch->mChargeQ;
    srslice.fmatch.light_pe  = fmatch->mLightPE;
    srslice.fmatch.score     = fmatch->mScore;
    srslice.fmatch.scr_y     = fmatch->mScr_y;
    srslice.fmatch.scr_z     = fmatch->mScr_z;
    srslice.fmatch.scr_rr    = fmatch->mScr_rr;
    srslice.fmatch.scr_ratio = fmatch->mScr_ratio;
    srslice.fmatch.chargeXYZ = fmatch->mChargeXYZ;
    srslice.fmatch.lightXYZ  = fmatch->mLightXYZ;
  }

  //......................................................................
  void FillSliceFlashMatch(const sbn::SimpleFlashMatch* fmatch /* can be nullptr */,
                           caf::SRSlice& srslice,
                           bool allowEmpty)
  {
    caf::SRFlashMatch flashmatch;
    if (fmatch == nullptr) {
      flashmatch.present = false;
      return;
    }
    flashmatch.present   = fmatch->mPresent;
    flashmatch.time      = fmatch->mTime;
    flashmatch.charge_q  = fmatch->mChargeQ;
    flashmatch.light_pe  = fmatch->mLightPE;
    flashmatch.score     = fmatch->mScore;
    flashmatch.scr_y     = fmatch->mScr_y;
    flashmatch.scr_z     = fmatch->mScr_z;
    flashmatch.scr_rr    = fmatch->mScr_rr;
    flashmatch.scr_ratio = fmatch->mScr_ratio;
    flashmatch.chargeXYZ = fmatch->mChargeXYZ;
    flashmatch.lightXYZ  = fmatch->mLightXYZ;
    srslice.fmatch_a = flashmatch;
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
    srslice.fmatch_b.present   = fmatch->mPresent;
    srslice.fmatch_b.time      = fmatch->mTime;
    srslice.fmatch_b.charge_q  = fmatch->mChargeQ;
    srslice.fmatch_b.light_pe  = fmatch->mLightPE;
    srslice.fmatch_b.score     = fmatch->mScore;
    srslice.fmatch_b.scr_y     = fmatch->mScr_y;
    srslice.fmatch_b.scr_z     = fmatch->mScr_z;
    srslice.fmatch_b.scr_rr    = fmatch->mScr_rr;
    srslice.fmatch_b.scr_ratio = fmatch->mScr_ratio;
    srslice.fmatch_b.chargeXYZ = fmatch->mChargeXYZ;
    srslice.fmatch_b.lightXYZ  = fmatch->mLightXYZ;
  }

  //......................................................................

} // end namespace
