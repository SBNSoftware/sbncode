//////////////////////////////////////////////////////////////////////
// \file    FillFlashMatch.cxx
// \brief   Fill flash match info from multiple algs
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "FillFlashMatch.h"

namespace caf
{


  //......................................................................
  void FillSliceFlashMatchA(const anab::T0 *fmatch /* can be NULL */,
                           caf::SRSlice &srslice,
                           bool allowEmpty)
  {
    if (fmatch != NULL) {
      srslice.fmatch.present = true;
      srslice.fmatch.time    = fmatch->Time();
      srslice.fmatch.score   = fmatch->TriggerConfidence();
      srslice.fmatch.pe      = fmatch->TriggerType();
    }
    else {
      srslice.fmatch.present = false;
    }
  }

  //......................................................................
  void FillSliceFlashMatch(const anab::T0 *fmatch /* can be NULL */,
			   caf::SRSlice &srslice,
			   bool allowEmpty)
  {
    caf::SRFlashMatch flashmatch;

    if (fmatch != NULL) {
      flashmatch.present = true;
      flashmatch.time    = fmatch->Time();
      flashmatch.score   = fmatch->TriggerConfidence();
      flashmatch.pe      = fmatch->TriggerType();
    }
    else {
      flashmatch.present = false;
    }
    srslice.fmatch_a = flashmatch;

  }

  //......................................................................
  void FillSliceFlashMatchB(const anab::T0 *fmatch /* can be NULL */,
			    caf::SRSlice &srslice,
			    bool allowEmpty)
  {
    if (fmatch != NULL) {
      srslice.fmatch_b.present = true;
      srslice.fmatch_b.time    = fmatch->Time();
      srslice.fmatch_b.score   = fmatch->TriggerConfidence();
      srslice.fmatch_b.pe      = fmatch->TriggerType();
    }
    else {
      srslice.fmatch_b.present = false;
    }
  }

  //......................................................................

} // end namespace
