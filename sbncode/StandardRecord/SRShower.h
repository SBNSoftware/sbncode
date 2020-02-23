////////////////////////////////////////////////////////////////////////
// \file    SRShower.h
////////////////////////////////////////////////////////////////////////
#ifndef SRSHOWER_H
#define SRSHOWER_H

#include "SRVector3D.h"


namespace caf
{
  /// Representation of a rb::Shower, knows energy and direction, but not a list
  /// of hits.
  class SRShower
    {
    public:
      SRShower();
      ~SRShower(){  };
      int bestplane;     ///< shower best plane
      float dEdx;        ///< shower calculated dEdx at best plane
      float energy;      ///< shower calculated energy at best plane
      float len;         ///< shower length [cm]
      float openAngle;   ///< shower opening angle [rad]
      SRVector3D dir;    ///< direction cosines at the start of the shower
      SRVector3D start;  ///< shower start point in detector coordinates [cm]
    };

  // TO DO:
  // * find dEdx and energy units
  // * description of best plane

} // end namespace

#endif // SRSHOWER_H
//////////////////////////////////////////////////////////////////////////////
