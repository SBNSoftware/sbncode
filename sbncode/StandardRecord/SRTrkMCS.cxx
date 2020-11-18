////////////////////////////////////////////////////////////////////////
// \file    SRTrkMCS.cxx
// \brief   An SRTrkMCS contains reco momentum from MCS track fit
////////////////////////////////////////////////////////////////////////
#include "sbncode/StandardRecord/SRTrkMCS.h"

#include <limits>

namespace caf
{

  SRTrkMCS::SRTrkMCS():
    fwdP_muon(std::numeric_limits<double>::signaling_NaN()),
    fwdP_pion(std::numeric_limits<double>::signaling_NaN()),
    fwdP_kaon(std::numeric_limits<double>::signaling_NaN()),
    fwdP_proton(std::numeric_limits<double>::signaling_NaN()),
    fwdP_err_muon(std::numeric_limits<double>::signaling_NaN()),
    fwdP_err_pion(std::numeric_limits<double>::signaling_NaN()),
    fwdP_err_kaon(std::numeric_limits<double>::signaling_NaN()),
    fwdP_err_proton(std::numeric_limits<double>::signaling_NaN()),
    bwdP_muon(std::numeric_limits<double>::signaling_NaN()),
    bwdP_pion(std::numeric_limits<double>::signaling_NaN()),
    bwdP_kaon(std::numeric_limits<double>::signaling_NaN()),
    bwdP_proton(std::numeric_limits<double>::signaling_NaN()),
    bwdP_err_muon(std::numeric_limits<double>::signaling_NaN()),
    bwdP_err_pion(std::numeric_limits<double>::signaling_NaN()),
    bwdP_err_kaon(std::numeric_limits<double>::signaling_NaN()),
    bwdP_err_proton(std::numeric_limits<double>::signaling_NaN()),
    is_bwd_muon(std::numeric_limits<bool>::signaling_NaN()),
    is_bwd_pion(std::numeric_limits<bool>::signaling_NaN()),
    is_bwd_kaon(std::numeric_limits<bool>::signaling_NaN()),
    is_bwd_proton(std::numeric_limits<bool>::signaling_NaN())
  {  }

  SRTrkMCS::~SRTrkMCS(){  }

  void SRTrkMCS::setDefault()
  {
    fwdP_muon        = -5.0;
    fwdP_pion        = -5.0;
    fwdP_kaon        = -5.0;
    fwdP_proton      = -5.0;
    fwdP_err_muon    = -5.0;
    fwdP_err_pion    = -5.0;
    fwdP_err_kaon    = -5.0;
    fwdP_err_proton  = -5.0;
    bwdP_muon        = -5.0;
    bwdP_pion        = -5.0;
    bwdP_kaon        = -5.0;
    bwdP_proton      = -5.0;
    bwdP_err_muon    = -5.0;
    bwdP_err_pion    = -5.0;
    bwdP_err_kaon    = -5.0;
    bwdP_err_proton  = -5.0;
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
