////////////////////////////////////////////////////////////////////////
// \file    SRTrkRange.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRKRANGE_H
#define SRTRKRANGE_H

namespace caf
{
  /// Representation of the reco momentum and PID a rb::Track from range
  class SRTrkRange
    {
    public:

      SRTrkRange();
      virtual ~SRTrkRange();

      float    p_muon;   ///< momentum estimate from trk range (muon hypothesis)
      float    p_proton; ///< momentum estimate from trk range (proton hypothesis)

      void setDefault();
    };

} // end namespace

#endif // SRTRKRANGE_H
//////////////////////////////////////////////////////////////////////////////
