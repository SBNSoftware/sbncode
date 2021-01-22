////////////////////////////////////////////////////////////////////////
// \file    SRShower.h
////////////////////////////////////////////////////////////////////////
#ifndef SRSHOWER_H
#define SRSHOWER_H

#include "sbncode/StandardRecord/SRVector3D.h"
#include "sbncode/StandardRecord/SRShowerSelection.h"
#include "sbncode/StandardRecord/SRTrackTruth.h"

namespace caf
{
  /// Representation of a rb::Shower, knows energy and direction, but not a list
  /// of hits.
  class SRShower
    {
    public:
      SRShower();
      ~SRShower(){  }
      int bestplane;                ///< shower best reconstructed plane
      float bestplane_dEdx;         ///< shower dEdx at best plane [MeV/cm]
      float bestplane_energy;       ///< shower energy at best plane [MeV]
      float conversion_gap;         ///< shower start and vertex position difference [cm]
      float density;                ///< shower density [MeV/cm]
      float len;                    ///< shower length [cm]
      float open_angle;             ///< shower opening angle [rad]
      float dEdx_plane0;            ///< shower calculated dEdx at each plane [MeV/cm]
      float dEdx_plane1;            ///< shower calculated dEdx at each plane [MeV/cm]
      float dEdx_plane2;            ///< shower calculated dEdx at each plane [MeV/cm]
      float energy_plane0;          ///< shower calculated energy at each plane [MeV]
      float energy_plane1;          ///< shower calculated energy at each plane [MeV]
      float energy_plane2;          ///< shower calculated energy at each plane [MeV]
      unsigned int nHits_plane0;    ///< Number of hits associated to the shower
      unsigned int nHits_plane1;    ///< Number of hits associated to the shower
      unsigned int nHits_plane2;    ///< Number of hits associated to the shower
      float wirePitch_plane0;       ///< Wire pitch corrected for the angle of the shower [cm]
      float wirePitch_plane1;       ///< Wire pitch corrected for the angle of the shower [cm]
      float wirePitch_plane2;       ///< Wire pitch corrected for the angle of the shower [cm]
      std::vector<float> dEdx;      ///< shower calculated dEdx at each plane [MeV/cm]
      std::vector<float> energy;    ///< shower calculated energy at each plane [MeV]
      SRVector3D dir;               ///< direction cosines at the start of the shower
      SRVector3D start;             ///< shower start point in detector coordinates [cm]
      SRVector3D end;               ///< shower end point (start+len*dir) in detector coordinates [cm]

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
