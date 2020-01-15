////////////////////////////////////////////////////////////////////////
// \file    SRTrack.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRACK_H
#define SRTRACK_H

/* #include "SRVector3D.h" */


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
      //      SRVector3D     start;        ///< Shower start point in detector coordinates. [cm]
      float          len;          ///< track length [cm]

    };

} // end namespace

#endif // SRTRACK_H
//////////////////////////////////////////////////////////////////////////////
