////////////////////////////////////////////////////////////////////////
// \file    SRShower.h
////////////////////////////////////////////////////////////////////////
#ifndef SRSHOWER_H
#define SRSHOWER_H

#include "SRVector3D.h"
#include "SRShowerSelection.h"
#include "SRTrackTruth.h"

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
      double bestplane_dEdx;     ///< shower dEdx at best plane [MeV/cm]
      double bestplane_energy;   ///< shower energy at best plane [MeV]
      float conversion_gap;      ///< shower start and vertex position difference [cm]
      float density;             ///< shower density [MeV/cm]
      float len;                 ///< shower length [cm]
      float open_angle;          ///< shower opening angle [rad]
      std::vector<double> dEdx;     ///< shower calculated dEdx at best plane [MeV/cm]
      std::vector<double> energy;   ///< shower calculated energy at best plane [MeV]
      SRVector3D dir;               ///< direction cosines at the start of the shower
      SRVector3D start;             ///< shower start point in detector coordinates [cm]

      int            ID;          ///< ID of this shower (taken from the pandora particle "ID" of this PFP)
      std::vector<int> daughters; ///< ID's of daughters of this shower
      int parent;                 ///< ID of parent particle of this shower

      SRShowerSelection selVars;
      SRTrackTruth   truth;        ///< truth information TODO: make seperate showe info class

      bool parent_is_primary;
      
      int slcID;            ///< ID of the slice that this shower is in
      unsigned producer;    ///< Index of the producer that produced this object. In ICARUS, this is the same as the cryostat.

    };

} // end namespace

#endif // SRSHOWER_H
//////////////////////////////////////////////////////////////////////////////
