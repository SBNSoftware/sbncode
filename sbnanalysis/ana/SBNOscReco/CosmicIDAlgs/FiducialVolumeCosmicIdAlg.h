#ifndef FIDUCIALVOLUMECOSMICIDALG_H_SEEN
#define FIDUCIALVOLUMECOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// FiducialVolumeCosmicIdAlg.h
//
// Functions for fiducial volume cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

// sbncode
#include "core/ProviderManager.hh"

// c++
#include <vector>


namespace ana{

  class FiducialVolumeCosmicIdAlg {
  public:

    struct Fiducial {
      using Name = fhicl::Name;

      fhicl::Atom<double> MinX { Name("MinX") };
      fhicl::Atom<double> MinY { Name("MinY") };
      fhicl::Atom<double> MinZ { Name("MinZ") };
      fhicl::Atom<double> MaxX { Name("MaxX") };
      fhicl::Atom<double> MaxY { Name("MaxY") };
      fhicl::Atom<double> MaxZ { Name("MaxZ") };

    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<Fiducial> FiducialCuts {
        Name("FiducialCuts"),
        Comment("Fiducial volume cuts (cm)")
      };

    };

    FiducialVolumeCosmicIdAlg(const core::ProviderManager &manager, const Config& config);

    FiducialVolumeCosmicIdAlg(const core::ProviderManager &manager, const fhicl::ParameterSet& pset) :
      FiducialVolumeCosmicIdAlg(manager, fhicl::Table<Config>(pset, {})()) {}

    FiducialVolumeCosmicIdAlg();

    ~FiducialVolumeCosmicIdAlg();

    void reconfigure(const core::ProviderManager &manager, const Config& config);

    // Check if point in fiducial volume used by this algorithm
    bool InFiducial(geo::Point_t point);

    // Check both start and end points of track are in fiducial volume
    bool FiducialVolumeCosmicId(recob::Track track);

  private:
    std::vector<geo::BoxBoundedGeo> fFiducialVolumes;

    double fMinX;
    double fMinY;
    double fMinZ;
    double fMaxX;
    double fMaxY;
    double fMaxZ;
  };

}

#endif
