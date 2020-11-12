////////////////////////////////////////////////////////////////////////
// \file    SRTrack.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRACK_H
#define SRTRACK_H

/* #include "SRVector3D.h" */

#include "SRTrackTruth.h"
#include "SRTrkChi2PID.h"
#include "SRTrkMCS.h"
#include "SRTrkRange.h"
#include "SRCRTHitMatch.h"
#include "SRTrackCalo.h"

#include "SRVector3D.h"
#include "SREnums.h"
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

      unsigned producer;    ///< Index of the producer that produced this object. In ICARUS, this is the same as the cryostat.
      unsigned short npts;         ///< number of points (recob Track.NPoints)
      float          len;          ///< track length [cm]
      float          costh;       ///< Costh of start direction of track
      float          phi;         ///< Angle of the start direction of the track in the x-y plane
      SRVector3D     start;       ///< Start point of track
      SRVector3D     end;         ///< End point of track
      int            ID;          ///< ID of this track (taken from the pandora particle "ID" of this track)

      int            nchi2pid;
      std::vector<SRTrkChi2PID> chi2pid; ///< 3-item list of larana Chi2 Particle PID on each plane ordered (1st ind., 2nd ind., coll)
      int            ncalo;
      std::vector<SRTrackCalo> calo; ///< 3-item list of Calorimetry information on each plane ordered (1st ind., 2nd ind., coll)
      Plane_t            bestplane;   ///< Plane index with the most hits. -1 if no calorimetry

      SRTrkMCS       mcsP;
      SRTrkRange     rangeP;

      SRTrackTruth   truth;        ///< truth information
      SRCRTHitMatch  crthit;       ///< CRT Hit match

      // TO DO: Move the following into SRObjects      

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

      //      CRTMatch       crt_match;   ///< Matching to CRT information
      std::vector<int> daughters; ///< ID's of daughters of this track
      int parent;                 ///< ID of parent particle of this track
      bool parent_is_primary;

      int slcID;


    };

} // end namespace

#endif // SRTRACK_H
//////////////////////////////////////////////////////////////////////////////
