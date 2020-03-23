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
      int bestplane;             ///< shower best reconstructed plane
      double bestplane_dEdx;     ///< shower dEdx at best plane
      double bestplane_energy;   ///< shower energy at best plane
      float conversion_gap;      ///< shower start and vertex position difference [cm]
      float density;             ///< shower density [MeV/cm]
      float len;                 ///< shower length [cm]
      float open_angle;          ///< shower opening angle [rad]
      std::vector<double> dEdx;     ///< shower calculated dEdx at best plane []
      std::vector<double> energy;   ///< shower calculated energy at best plane [MeV]
      SRVector3D dir;               ///< direction cosines at the start of the shower
      SRVector3D start;             ///< shower start point in detector coordinates [cm]
    };

  // TO DO:
  // * find dEdx units

} // end namespace

#endif // SRSHOWER_H
//////////////////////////////////////////////////////////////////////////////
