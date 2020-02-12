////////////////////////////////////////////////////////////////////////
// \file    SRTrack.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRACK_H
#define SRTRACK_H

/* #include "SRVector3D.h" */

#include "SRTrackTruth.h"
#include "SRTrkChi2PID.h"
#include <TVector3.h>
#include <vector>

namespace caf
{
  /// Representation of a rb::Track, knows energy and direction, but not a list
  /// of hits.
  class SRTrack
    {
    public:


      SRTrack();
      ~SRTrack(){  };

      unsigned short npts;         ///< number of points (recob Track.NPoints)
      float          len;          ///< track length [cm]
      float          costh;       ///< Costh of start direction of track
      TVector3       start;       ///< Start point of track
      TVector3       end;         ///< End point of track
      int            ID;          ///< ID of this track (taken from the pandora particle "ID" of this track)

      SRTrkChi2PID   chi2pid;     ///< larana Chi2 Particle PID

      SRTrackTruth   truth;        ///< truth information

      // TO DO: Move the following into SRObjects

      float          range_momentum_muon; ///< Momentum from range (muon hypothesis)
      float          range_momentum_proton; ///< Momentum from range (proton hypothesis)
      
      /// Information on the result of MCS fitting for momentum
      class MCSFitResult {
        public:
        float fwd_momentum; //!< Momentum result of fitting the track from start -> end [GeV/c]
        float fwd_momentum_err; //!< Error on momentum result of fitting the track from start -> end [GeV/c]
        float bwd_momentum; //!< Momentum result of fitting the track from end -> start [GeV/c]
        float bwd_momentum_err; //!< Error on momentum result of fitting the track from end -> start [GeV/c]
        bool is_bwd; //!< Whether the MCS fit thinks the track is backwards
      };

      /* struct CRTMatch { */
        /* struct Track { */
        /*   bool present; //!< Whether a CRT track match exists */
        /*   float time; //!< time of the CRT track [mus -- t=0 is spill time] */
        /*   float angle; //!< Angle of the match between the TPC track and the CRT track [rad] */
        /* }; */
 
        /* struct Hit { */
        /*   bool present; //!< Whether a CRT hit match exists */
        /*   float distance; //!< Distance of closest approach between CRT hit and projected TPC track [cm] */
        /*   float time; //!< Time of CRT hit [mus -- t=0 is spill time] */
        /* }; */

      /*   Track track; //!< CRT track match */
      /*   Hit   hit;   //!< CRT hit match */
      /* }; */
      MCSFitResult   mcs_muon;     ///< MCS fit result using muon particle hypothesis 
      MCSFitResult   mcs_pion;     ///< MCS fit result using pion particle hypothesis 
      MCSFitResult   mcs_kaon;     ///< MCS fit result using kaon particle hypothesis 
      MCSFitResult   mcs_proton;   ///< MCS fit result using proton particle hypothesis 

      //      CRTMatch       crt_match;   ///< Matching to CRT information
      std::vector<int> daughters; ///< ID's of daughters of this track


    };

} // end namespace

#endif // SRTRACK_H
//////////////////////////////////////////////////////////////////////////////
