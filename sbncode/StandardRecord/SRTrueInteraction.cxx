////////////////////////////////////////////////////////////////////////
// \file    SRTrueInteraction.cxx
// \brief   SRTrueInteraction holds true interaction info.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "SRTrueInteraction.h"


namespace caf
{

  SRTrueInteraction::SRTrueInteraction():
    pdg(-1),
    interaction_id(-1),
    inttype(0),
    iscc(true),
    isvtxcont(std::numeric_limits<bool>::signaling_NaN()),
    E(std::numeric_limits<float>::signaling_NaN()),
    visE(std::numeric_limits<float>::signaling_NaN()),
    visEinslc(std::numeric_limits<float>::signaling_NaN()),
    visEcosmic(std::numeric_limits<float>::signaling_NaN()),
    eff(std::numeric_limits<float>::signaling_NaN()),
    pur(std::numeric_limits<float>::signaling_NaN()),
    time(std::numeric_limits<float>::signaling_NaN()),
    genweight(std::numeric_limits<float>::signaling_NaN()),
    xsec(std::numeric_limits<float>::signaling_NaN()),
    q2(std::numeric_limits<float>::signaling_NaN()),
    x(std::numeric_limits<float>::signaling_NaN()),
    y(std::numeric_limits<float>::signaling_NaN()),
    w2(std::numeric_limits<float>::signaling_NaN()),
    det(kUNKNOWN),
    mode(kOther),
    generator(kUnknownGenerator),
    genVersion(),
    prim()
  {  }


} // end namespace caf
////////////////////////////////////////////////////////////////////////
