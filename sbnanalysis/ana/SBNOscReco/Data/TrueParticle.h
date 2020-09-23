#ifndef _sbnumurecodata_TrueParticle_hh
#define _sbnumurecodata_TrueParticle_hh

#include <vector>

#include "TVector3.h"

namespace numu {

enum Wall {
  wNone=0,
  wTop=1,
  wBottom=2,
  wLeft=3,
  wRight=4,
  wFront=5,
  wBack=6
};

/**
* Information on true particles
*/
struct TrueParticle {
  TVector3 start_momentum; //!< Particle directional momentum for first trajectory point inside TPC AV [GeV] 
  TVector3 end_momentum; //!< Particle directional momentum for last trajectory point inside TPC AV [GeV]
  float start_energy; //!< Particle energy for first point inside TPC AV [GeV]
  float end_energy; //!< Particle energy for last point inside TPC AV [GeV]
  float deposited_energy; //!< Total particle energy depositive in TPC AV [GeV]
  
  TVector3 start; //!< start position of track
  TVector3 end; //!< end position of track
  float start_time; //!< start time of track
  float end_time; //!< end time of track 

  Wall wall_enter; //!< the face of the TPC that the particle crosses on enter
  Wall wall_exit; //!< the face of the TPC that the particle crosses on exit
  bool contained_in_cryo; //!< is it contained a single cryostat?
  bool contained_in_tpc; //!< is it contained in a single TPC?
  bool crosses_tpc; //!< does it cross a tpc?
  bool is_contained; //!< is it contained in a single cryostat active volume

  float length; //!< Length of track contained in any TPC active volume [cm]
  int pdgid; //!< Particle ID code

  bool is_cosmic; //!< Whether this particle is of cosmic origin

  int ID; //!< ID/index of this particle (taken from MCParticle ID)
  
  TrueParticle():
    start_momentum(-1, -1, -1),
    end_momentum(-1, -1, -1),
    start_energy(-1),
    end_energy(-1),
    start(-999, -999, -999),
    end(-999, -999, -999),
    start_time(-1),
    end_time(-1),
    wall_enter(numu::wNone),
    wall_exit(numu::wNone),
    contained_in_cryo(false),
    contained_in_tpc(false),
    crosses_tpc(false),
    is_contained(false),
    length(-1),
    pdgid(-1),
    ID(-1)
    {}
};
} // namespace numu
#endif
