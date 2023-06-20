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
// OpFlashes

  void FillSliceFlashMatchOp(const sbn::SimpleFlashMatch* fmatchop /* can be nullptr */,
                           caf::SRSlice& srslice,
                           bool allowEmpty)
  {
    if (fmatchop == nullptr) {
      srslice.fmatchop.present = false;
      return;
    }
    srslice.fmatchop.present    = fmatchop->present;
    srslice.fmatchop.time       = fmatchop->time;
    srslice.fmatchop.chargeQ    = fmatchop->charge.q;
    srslice.fmatchop.lightPE    = fmatchop->light.pe;
    srslice.fmatchop.score      = fmatchop->score.total;
    srslice.fmatchop.scoreY     = fmatchop->score.y;
    srslice.fmatchop.scoreZ     = fmatchop->score.z;
    srslice.fmatchop.scoreRR    = fmatchop->score.rr;
    srslice.fmatchop.scoreRatio = fmatchop->score.ratio;
    srslice.fmatchop.scoreSlope = fmatchop->score.slope;
    srslice.fmatchop.scorePEToQ = fmatchop->score.petoq;
    srslice.fmatchop.chargeCenter = fmatchop->charge.center;
    srslice.fmatchop.chargeWidth = fmatchop->charge.width;
    srslice.fmatchop.lightCenter  = fmatchop->light.center;
    srslice.fmatchop.lightWidth  = fmatchop->light.width;
  }

  //......................................................................
  void FillSliceFlashMatchOpA(const sbn::SimpleFlashMatch* fmatchop /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatchop == nullptr) {
      srslice.fmatchop_a.present = false;
      return;
    }
    srslice.fmatchop_a.present    = fmatchop->present;
    srslice.fmatchop_a.time       = fmatchop->time;
    srslice.fmatchop_a.chargeQ    = fmatchop->charge.q;
    srslice.fmatchop_a.lightPE    = fmatchop->light.pe;
    srslice.fmatchop_a.score      = fmatchop->score.total;
    srslice.fmatchop_a.scoreY     = fmatchop->score.y;
    srslice.fmatchop_a.scoreZ     = fmatchop->score.z;
    srslice.fmatchop_a.scoreRR    = fmatchop->score.rr;
    srslice.fmatchop_a.scoreRatio = fmatchop->score.ratio;
    srslice.fmatchop_a.scoreSlope = fmatchop->score.slope;
    srslice.fmatchop_a.scorePEToQ = fmatchop->score.petoq;
    srslice.fmatchop_a.chargeCenter = fmatchop->charge.center;
    srslice.fmatchop_a.chargeWidth  = fmatchop->charge.width;
    srslice.fmatchop_a.lightCenter  = fmatchop->light.center;
    srslice.fmatchop_a.lightWidth   = fmatchop->light.width;
  }

  //......................................................................
  void FillSliceFlashMatchOpB(const sbn::SimpleFlashMatch* fmatchop /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatchop == nullptr) {
      srslice.fmatchop_b.present = false;
      return;
    }
    srslice.fmatchop_b.present    = fmatchop->present;
    srslice.fmatchop_b.time       = fmatchop->time;
    srslice.fmatchop_b.chargeQ    = fmatchop->charge.q;
    srslice.fmatchop_b.lightPE    = fmatchop->light.pe;
    srslice.fmatchop_b.score      = fmatchop->score.total;
    srslice.fmatchop_b.scoreY     = fmatchop->score.y;
    srslice.fmatchop_b.scoreZ     = fmatchop->score.z;
    srslice.fmatchop_b.scoreRR    = fmatchop->score.rr;
    srslice.fmatchop_b.scoreRatio = fmatchop->score.ratio;
    srslice.fmatchop_b.scoreSlope = fmatchop->score.slope;
    srslice.fmatchop_b.scorePEToQ = fmatchop->score.petoq;
    srslice.fmatchop_b.chargeCenter = fmatchop->charge.center;
    srslice.fmatchop_b.chargeWidth  = fmatchop->charge.width;
    srslice.fmatchop_b.lightCenter  = fmatchop->light.center;
    srslice.fmatchop_b.lightWidth   = fmatchop->light.width;
  }

  //......................................................................
// XARAPUCA
  void FillSliceFlashMatchARA(const sbn::SimpleFlashMatch* fmatchara /* can be nullptr */,
                           caf::SRSlice& srslice,
                           bool allowEmpty)
  {
    if (fmatchara == nullptr) {
      srslice.fmatchara.present = false;
      return;
    }
    srslice.fmatchara.present    = fmatchara->present;
    srslice.fmatchara.time       = fmatchara->time;
    srslice.fmatchara.chargeQ    = fmatchara->charge.q;
    srslice.fmatchara.lightPE    = fmatchara->light.pe;
    srslice.fmatchara.score      = fmatchara->score.total;
    srslice.fmatchara.scoreY     = fmatchara->score.y;
    srslice.fmatchara.scoreZ     = fmatchara->score.z;
    srslice.fmatchara.scoreRR    = fmatchara->score.rr;
    srslice.fmatchara.scoreRatio = fmatchara->score.ratio;
    srslice.fmatchara.scoreSlope = fmatchara->score.slope;
    srslice.fmatchara.scorePEToQ = fmatchara->score.petoq;
    srslice.fmatchara.chargeCenter = fmatchara->charge.center;
    srslice.fmatchara.chargeWidth = fmatchara->charge.width;
    srslice.fmatchara.lightCenter  = fmatchara->light.center;
    srslice.fmatchara.lightWidth  = fmatchara->light.width;
  }

  //......................................................................
  void FillSliceFlashMatchARAA(const sbn::SimpleFlashMatch* fmatchara /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatchara == nullptr) {
      srslice.fmatchara_a.present = false;
      return;
    }
    srslice.fmatchara_a.present    = fmatchara->present;
    srslice.fmatchara_a.time       = fmatchara->time;
    srslice.fmatchara_a.chargeQ    = fmatchara->charge.q;
    srslice.fmatchara_a.lightPE    = fmatchara->light.pe;
    srslice.fmatchara_a.score      = fmatchara->score.total;
    srslice.fmatchara_a.scoreY     = fmatchara->score.y;
    srslice.fmatchara_a.scoreZ     = fmatchara->score.z;
    srslice.fmatchara_a.scoreRR    = fmatchara->score.rr;
    srslice.fmatchara_a.scoreRatio = fmatchara->score.ratio;
    srslice.fmatchara_a.scoreSlope = fmatchara->score.slope;
    srslice.fmatchara_a.scorePEToQ = fmatchara->score.petoq;
    srslice.fmatchara_a.chargeCenter = fmatchara->charge.center;
    srslice.fmatchara_a.chargeWidth  = fmatchara->charge.width;
    srslice.fmatchara_a.lightCenter  = fmatchara->light.center;
    srslice.fmatchara_a.lightWidth   = fmatchara->light.width;
  }

  //......................................................................
  void FillSliceFlashMatchARAB(const sbn::SimpleFlashMatch* fmatchara /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatchara == nullptr) {
      srslice.fmatchara_b.present = false;
      return;
    }
    srslice.fmatchara_b.present    = fmatchara->present;
    srslice.fmatchara_b.time       = fmatchara->time;
    srslice.fmatchara_b.chargeQ    = fmatchara->charge.q;
    srslice.fmatchara_b.lightPE    = fmatchara->light.pe;
    srslice.fmatchara_b.score      = fmatchara->score.total;
    srslice.fmatchara_b.scoreY     = fmatchara->score.y;
    srslice.fmatchara_b.scoreZ     = fmatchara->score.z;
    srslice.fmatchara_b.scoreRR    = fmatchara->score.rr;
    srslice.fmatchara_b.scoreRatio = fmatchara->score.ratio;
    srslice.fmatchara_b.scoreSlope = fmatchara->score.slope;
    srslice.fmatchara_b.scorePEToQ = fmatchara->score.petoq;
    srslice.fmatchara_b.chargeCenter = fmatchara->charge.center;
    srslice.fmatchara_b.chargeWidth  = fmatchara->charge.width;
    srslice.fmatchara_b.lightCenter  = fmatchara->light.center;
    srslice.fmatchara_b.lightWidth   = fmatchara->light.width;
  }

  //......................................................................
// XARAPUCA OpFlash
  void FillSliceFlashMatchOpARA(const sbn::SimpleFlashMatch* fmatchopara /* can be nullptr */,
                           caf::SRSlice& srslice,
                           bool allowEmpty)
  {
    if (fmatchopara == nullptr) {
      srslice.fmatchopara.present = false;
      return;
    }
    srslice.fmatchopara.present    = fmatchopara->present;
    srslice.fmatchopara.time       = fmatchopara->time;
    srslice.fmatchopara.chargeQ    = fmatchopara->charge.q;
    srslice.fmatchopara.lightPE    = fmatchopara->light.pe;
    srslice.fmatchopara.score      = fmatchopara->score.total;
    srslice.fmatchopara.scoreY     = fmatchopara->score.y;
    srslice.fmatchopara.scoreZ     = fmatchopara->score.z;
    srslice.fmatchopara.scoreRR    = fmatchopara->score.rr;
    srslice.fmatchopara.scoreRatio = fmatchopara->score.ratio;
    srslice.fmatchopara.scoreSlope = fmatchopara->score.slope;
    srslice.fmatchopara.scorePEToQ = fmatchopara->score.petoq;
    srslice.fmatchopara.chargeCenter = fmatchopara->charge.center;
    srslice.fmatchopara.chargeWidth = fmatchopara->charge.width;
    srslice.fmatchopara.lightCenter  = fmatchopara->light.center;
    srslice.fmatchopara.lightWidth  = fmatchopara->light.width;
  }

  //......................................................................
  void FillSliceFlashMatchOpARAA(const sbn::SimpleFlashMatch* fmatchopara /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatchopara == nullptr) {
      srslice.fmatchopara_a.present = false;
      return;
    }
    srslice.fmatchopara_a.present    = fmatchopara->present;
    srslice.fmatchopara_a.time       = fmatchopara->time;
    srslice.fmatchopara_a.chargeQ    = fmatchopara->charge.q;
    srslice.fmatchopara_a.lightPE    = fmatchopara->light.pe;
    srslice.fmatchopara_a.score      = fmatchopara->score.total;
    srslice.fmatchopara_a.scoreY     = fmatchopara->score.y;
    srslice.fmatchopara_a.scoreZ     = fmatchopara->score.z;
    srslice.fmatchopara_a.scoreRR    = fmatchopara->score.rr;
    srslice.fmatchopara_a.scoreRatio = fmatchopara->score.ratio;
    srslice.fmatchopara_a.scoreSlope = fmatchopara->score.slope;
    srslice.fmatchopara_a.scorePEToQ = fmatchopara->score.petoq;
    srslice.fmatchopara_a.chargeCenter = fmatchopara->charge.center;
    srslice.fmatchopara_a.chargeWidth  = fmatchopara->charge.width;
    srslice.fmatchopara_a.lightCenter  = fmatchopara->light.center;
    srslice.fmatchopara_a.lightWidth   = fmatchopara->light.width;
  }

  //......................................................................
  void FillSliceFlashMatchOpARAB(const sbn::SimpleFlashMatch* fmatchopara /* can be nullptr */,
                            caf::SRSlice& srslice,
                            bool allowEmpty)
  {
    if (fmatchopara == nullptr) {
      srslice.fmatchopara_b.present = false;
      return;
    }
    srslice.fmatchopara_b.present    = fmatchopara->present;
    srslice.fmatchopara_b.time       = fmatchopara->time;
    srslice.fmatchopara_b.chargeQ    = fmatchopara->charge.q;
    srslice.fmatchopara_b.lightPE    = fmatchopara->light.pe;
    srslice.fmatchopara_b.score      = fmatchopara->score.total;
    srslice.fmatchopara_b.scoreY     = fmatchopara->score.y;
    srslice.fmatchopara_b.scoreZ     = fmatchopara->score.z;
    srslice.fmatchopara_b.scoreRR    = fmatchopara->score.rr;
    srslice.fmatchopara_b.scoreRatio = fmatchopara->score.ratio;
    srslice.fmatchopara_b.scoreSlope = fmatchopara->score.slope;
    srslice.fmatchopara_b.scorePEToQ = fmatchopara->score.petoq;
    srslice.fmatchopara_b.chargeCenter = fmatchopara->charge.center;
    srslice.fmatchopara_b.chargeWidth  = fmatchopara->charge.width;
    srslice.fmatchopara_b.lightCenter  = fmatchopara->light.center;
    srslice.fmatchopara_b.lightWidth   = fmatchopara->light.width;
  }

  //......................................................................



} // end namespace
