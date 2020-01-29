////////////////////////////////////////////////////////////////////////
// \file    SRTrack.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRACK_H
#define SRTRACK_H

/* #include "SRVector3D.h" */

#include <TVector3.h>
#include <vector>

namespace caf
{
  /// Representation of a rb::Track, knows energy and direction, but not a list
  /// of hits.
  class SRTrack
    {
    public:
      /// Truth information associated with a reconstructed track
      struct Truth {
        float total_deposited_energy; //!< True total deposited energy associated with this Track [GeV]

        /// Match from a reconstructed track to a true particle 
        struct ParticleMatch {
          int G4ID; //!< ID of the particle match, taken from G4
          float energy; //!< Total energy matching between reco track and true particle [GeV]
        }; 

        std::vector<ParticleMatch> matches; //!< List of particle matches, sorted by most energy matched

        /// Get the G4ID of the primary (most energy) particle match to this track (returns -1 if no match)
        int GetPrimaryMatchID() const;

        /// Get the purity of the energy associated with this track relative to the best-matched particle (0 if no match) 
        float Purity() const;

      };

      /// Information on the result of MCS fitting for momentum
      struct MCSFitResult {
        float fwd_momentum; //!< Momentum result of fitting the track from start -> end [GeV/c]
        float fwd_momentum_err; //!< Error on momentum result of fitting the track from start -> end [GeV/c]
        float bwd_momentum; //!< Momentum result of fitting the track from end -> start [GeV/c]
        float bwd_momentum_err; //!< Error on momentum result of fitting the track from end -> start [GeV/c]
        bool is_bwd; //!< Whether the MCS fit thinks the track is backwards
      };

      struct CRTMatch {
        struct Track {
          bool present; //!< Whether a CRT track match exists
          float time; //!< time of the CRT track [mus -- t=0 is spill time]
          float angle; //!< Angle of the match between the TPC track and the CRT track [rad]
        };
 
        struct Hit {
          bool present; //!< Whether a CRT hit match exists
          float distance; //!< Distance of closest approach between CRT hit and projected TPC track [cm]
          float time; //!< Time of CRT hit [mus -- t=0 is spill time]
        };

        Track track; //!< CRT track match
        Hit   hit;   //!< CRT hit match
      };

      SRTrack();
      ~SRTrack(){  };

      unsigned short npts;         ///< number of points (recob Track.NPoints)
      //      SRVector3D     start;        ///< Shower start point in detector coordinates. [cm]
      float          len;          ///< track length [cm]
      Truth          truth;        ///< truth information
      MCSFitResult   mcs_muon;     ///< MCS fit result using muon particle hypothesis 
      MCSFitResult   mcs_pion;     ///< MCS fit result using pion particle hypothesis 
      MCSFitResult   mcs_kaon;     ///< MCS fit result using kaon particle hypothesis 
      MCSFitResult   mcs_proton;   ///< MCS fit result using proton particle hypothesis 
      float          range_momentum_muon; ///< Momentum from range calculation using muon particle hypothesis
      float          range_momentum_proton; ///< Momentum from range calculation using proton particle hypothesis
      
      float          chi2_muon;   ///< Reduced Chi2 value of dE/dx v. residual range profile compared to muon particle hypothesis
      float          chi2_pion;   ///< Reduced Chi2 value of dE/dx v. residual range profile compared to pion particle hypothesis
      float          chi2_kaon;   ///< Reduced Chi2 value of dE/dx v. residual range profile compared to kaon particle hypothesis
      float          chi2_proton; ///< Reduced Chi2 value of dE/dx v. residual range profile compared to proton particle hypothesis
      int            pid_ndof;    ///< Number of degress of freedom in Chi2 PID fit 

      float          costh;       ///< Costh of start direction of track
      TVector3       start;       ///< Start point of track
      TVector3       end;         ///< End point of track
      int            ID;          ///< ID of this track (taken from the pandora particle "ID" of this track)
      CRTMatch       crt_match;   ///< Matching to CRT information
      std::vector<int> track_daughters; ///< ID's of track daughters of this track

    };

} // end namespace

#endif // SRTRACK_H
//////////////////////////////////////////////////////////////////////////////
