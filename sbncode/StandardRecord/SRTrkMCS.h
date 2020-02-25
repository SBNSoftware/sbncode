////////////////////////////////////////////////////////////////////////
// \file    SRTrkMCS.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRKMCS_H
#define SRTRKMCS_H

namespace caf
{
  /// Representation of the reco momentum and PID a recob::Track for 
  /// muon, pion, kaon, and proton assumptions 
  class SRTrkMCS
    {
    public:

      SRTrkMCS();
      virtual ~SRTrkMCS();

      float fwdP_muon; //!< Momentum from start->end fit for muon [GeV/c]
      float fwdP_pion; //!< Momentum from start->end fit for pion [GeV/c]
      float fwdP_kaon; //!< Momentum from start->end fit for kaon [GeV/c]
      float fwdP_proton; //!< Momentum from start->end fit for proton [GeV/c]

      float fwdP_err_muon; //!< Error on momentum from start->end fit for muon [GeV/c]
      float fwdP_err_pion; //!< Error on momentum from start->end fit for pion [GeV/c]
      float fwdP_err_kaon; //!< Error on momentum from start->end fit for kaon [GeV/c]
      float fwdP_err_proton; //!< Error on momentum from start->end fit for proton [GeV/c]

      float bwdP_muon; //!< Momentum result from end->start fit for muon [Ge\V/c]
      float bwdP_pion; //!< Momentum result from end->start fit for pion [Ge\V/c]
      float bwdP_kaon; //!< Momentum result from end->start fit for kaon [Ge\V/c]
      float bwdP_proton; //!< Momentum result from end->start fit for proton [Ge\V/c]

      float bwdP_err_muon; //!< Error on momentum from end->start fit for muon [GeV/c]
      float bwdP_err_pion; //!< Error on momentum from end->start fit for pion [GeV/c]
      float bwdP_err_kaon; //!< Error on momentum from end->start fit for kaon [GeV/c]
      float bwdP_err_proton; //!< Error on momentum from end->start fit for proton [GeV/c]

      bool  is_bwd_muon; //!< MCS fit thinks the track is backwards for muon assumption
      bool  is_bwd_pion; //!< MCS fit thinks the track is backwards for pion assumption
      bool  is_bwd_kaon; //!< MCS fit thinks the track is backwards for kaon assumption
      bool  is_bwd_proton; //!< MCS fit thinks the track is backwards for proton assumption

      void setDefault();
    };

} // end namespace

#endif // SRTRKMCS_H
//////////////////////////////////////////////////////////////////////////////
