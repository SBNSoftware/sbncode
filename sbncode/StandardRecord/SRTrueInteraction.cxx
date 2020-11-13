////////////////////////////////////////////////////////////////////////
// \file    SRTrueInteraction.cxx
// \brief   SRTrueInteraction holds true interaction info.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "SRTrueInteraction.h"


namespace caf
{

  SRTrueInteraction::SRTrueInteraction():
    initpdg(-1),
    pdg(-1),
    inttype(0),
    index(-1),
    targetPDG(std::numeric_limits<int>::signaling_NaN()),
    genie_intcode(std::numeric_limits<int>::signaling_NaN()),
    parentPDG(std::numeric_limits<int>::signaling_NaN()),
    parentDecayMode(std::numeric_limits<int>::signaling_NaN()),
    isnc(std::numeric_limits<bool>::signaling_NaN()),
    iscc(std::numeric_limits<bool>::signaling_NaN()),
    isvtxcont(std::numeric_limits<bool>::signaling_NaN()),
    is_numucc_primary(std::numeric_limits<bool>::signaling_NaN()),

    E(std::numeric_limits<float>::signaling_NaN()),
    visE(std::numeric_limits<float>::signaling_NaN()),
    visEinslc(std::numeric_limits<float>::signaling_NaN()),
    visEcosmic(std::numeric_limits<float>::signaling_NaN()),
    time(std::numeric_limits<float>::signaling_NaN()),
    bjorkenX(std::numeric_limits<float>::signaling_NaN()),
    inelasticityY(std::numeric_limits<float>::signaling_NaN()),
    Q2(std::numeric_limits<float>::signaling_NaN()),
    q0(std::numeric_limits<float>::signaling_NaN()),
    modq(std::numeric_limits<float>::signaling_NaN()),
    q0_lab(std::numeric_limits<float>::signaling_NaN()),
    modq_lab(std::numeric_limits<float>::signaling_NaN()),
    w(std::numeric_limits<float>::signaling_NaN()),
    t(std::numeric_limits<float>::signaling_NaN()),
    eccqe(std::numeric_limits<float>::signaling_NaN()),
    baseline(std::numeric_limits<float>::signaling_NaN()),
    det(kUNKNOWN),
    mode(kOther),
    generator(kUnknownGenerator),
    genVersion(),
    prim()
  {  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
