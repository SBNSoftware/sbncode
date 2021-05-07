//////////////////////////////////////////////////////////////////////
// \file    FillFlashMatch.cxx
// \brief   Fill flash match info from multiple algs
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "FillFlashMatch.h"

namespace caf
{


  //......................................................................
  void FillSliceFlashMatchA(const sbn::SimpleFlashMatch &fmatch /* can be NULL */,
                           caf::SRSlice &srslice,
                           bool allowEmpty)
  {
    //if (fmatch.mPresent == true) {
      //std::cout << fmatch.mPresent << " " << fmatch.mTime << " " << fmatch.mScore << " " << fmatch.mScr_y << " " << fmatch.mScr_z << " " << fmatch.mScr_rr << " " << fmatch.mScr_ratio << " " << fmatch.mPE << std::endl;  
    srslice.fmatch.present = fmatch.mPresent;
    srslice.fmatch.time    = fmatch.mTime;
    srslice.fmatch.score   = fmatch.mScore;
    srslice.fmatch.scr_y   = fmatch.mScr_y;
    srslice.fmatch.scr_z   = fmatch.mScr_z;
    srslice.fmatch.scr_rr  = fmatch.mScr_rr;
    srslice.fmatch.scr_ratio = fmatch.mScr_ratio;
    srslice.fmatch.pe      = fmatch.mPE;
    srslice.fmatch.chargeXYZ = fmatch.mChargeXYZ;
    srslice.fmatch.lightXYZ = fmatch.mLightXYZ;
    //}
    //else {
    //srslice.fmatch.present = false;
    //}
  }

  //......................................................................
  void FillSliceFlashMatch(const sbn::SimpleFlashMatch &fmatch /* can be NULL */,
			   caf::SRSlice &srslice,
			   bool allowEmpty)
  {
    caf::SRFlashMatch flashmatch;
    
    //if (fmatch.mPresent == true) {
    flashmatch.present = fmatch.mPresent;
    flashmatch.time    = fmatch.mTime;
    flashmatch.score   = fmatch.mScore;
    flashmatch.scr_y   = fmatch.mScr_y;
    flashmatch.scr_z   = fmatch.mScr_z;
    flashmatch.scr_rr  = fmatch.mScr_rr;
    flashmatch.scr_ratio = fmatch.mScr_ratio;
    flashmatch.pe      = fmatch.mPE;
    flashmatch.chargeXYZ = fmatch.mChargeXYZ;
    flashmatch.lightXYZ = fmatch.mLightXYZ;
    //}
    //else {
    //flashmatch.present = false;
    //}
    srslice.fmatch_a = flashmatch;

  }

  //......................................................................
  void FillSliceFlashMatchB(const sbn::SimpleFlashMatch &fmatch /* can be NULL */,
			    caf::SRSlice &srslice,
			    bool allowEmpty)
  {
    //if (fmatch.mPresent == true) {
    srslice.fmatch_b.present = fmatch.mPresent;
    srslice.fmatch_b.time    = fmatch.mTime;
    srslice.fmatch_b.score   = fmatch.mScore;
    srslice.fmatch_b.scr_y   = fmatch.mScr_y;
    srslice.fmatch_b.scr_z   = fmatch.mScr_z;
    srslice.fmatch_b.scr_rr  = fmatch.mScr_rr;
    srslice.fmatch_b.scr_ratio = fmatch.mScr_ratio;
    srslice.fmatch_b.pe      = fmatch.mPE;
    srslice.fmatch_b.chargeXYZ = fmatch.mChargeXYZ;
    srslice.fmatch_b.lightXYZ = fmatch.mLightXYZ;
    //}
    //else {
    //srslice.fmatch_b.present = false;
    //}
  }

  //......................................................................

} // end namespace
