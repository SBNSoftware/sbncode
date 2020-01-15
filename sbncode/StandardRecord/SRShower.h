////////////////////////////////////////////////////////////////////////
// \file    SRShower.h
////////////////////////////////////////////////////////////////////////
#ifndef SRSHOWER_H
#define SRSHOWER_H

/* #include "SRVector3D.h" */


namespace caf
{
  /// Representation of a rb::Shower, knows energy and direction, but not a list
  /// of hits.
  class SRShower
    {
    public:
      SRShower();
      ~SRShower(){  };

      //      SRVector3D     start;        ///< Shower start point in detector coordinates. [cm]
      float          len;          ///< shower length [cm]

    };

} // end namespace

#endif // SRSHOWER_H
//////////////////////////////////////////////////////////////////////////////
