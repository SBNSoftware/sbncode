////////////////////////////////////////////////////////////////////////
// \file    SRTrueInteraction.cxx
// \brief   SRTrueInteraction holds true interaction info.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "sbncode/StandardRecord/SRTrueInteraction.h"


namespace caf
{

  SRTrueInteraction::SRTrueInteraction():
    initpdg(-1),
    pdg(-1),
    inttype(0),
    index(-1),
    targetPDG(-999),
    genie_intcode(-999),
    parentPDG(-999),
    parentDecayMode(-999),
    isnc(false),
    iscc(false),
    isvtxcont(false),
    is_numucc_primary(false),

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
