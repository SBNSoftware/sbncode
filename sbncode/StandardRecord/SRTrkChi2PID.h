////////////////////////////////////////////////////////////////////////
// \file    SRTrkChi2PID.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRKCHI2PID_H
#define SRTRKCHI2PID_H

namespace caf
{
  /// Track PID from dE/dx v. residual range Chi2
  class SRTrkChi2PID
    {
    public:

      SRTrkChi2PID();
      virtual ~SRTrkChi2PID();

      int            pdg;          ///< Min Chi2 pdg   
      int            pid_ndof;    ///< Number of degress of freedom in Chi2 PID fit 
      float          chi2_muon;   ///< dE/dx v. residual range Chi2 (muon hypothesis)
      float          chi2_pion;   ///< dE/dx v. residual range Chi2 (pion hypothesis)
      float          chi2_kaon;   ///< dE/dx v. residual range Chi2 (kaon hypothesis)
      float          chi2_proton; ///< dE/dx v. residual range Chi2 (proton hypothesis)

      float          pida;        ///< PIDA

      void setDefault();
    };

} // end namespace

#endif // SRTRKCHI2PID_H
//////////////////////////////////////////////////////////////////////////////
