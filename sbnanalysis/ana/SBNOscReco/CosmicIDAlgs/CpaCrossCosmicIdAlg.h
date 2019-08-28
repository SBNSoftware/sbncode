#ifndef CPACROSSCOSMICIDALG_H_SEEN
#define CPACROSSCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// CpaCrossCosmicIdAlg.h
//
// Functions for CPA stitching cosmic tagger
// Only useful if you forgot to run track stitching
// in your tracking
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// sbncode
#include "core/ProviderManager.hh"

// framework
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"

// c++
#include <vector>
#include <utility>


namespace ana{

  class CpaCrossCosmicIdAlg {
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

    struct BeamTime {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> BeamTimeMin {
        Name("BeamTimeMin"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeMax {
        Name("BeamTimeMax"),
        Comment("")
      };

    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> CpaStitchDistance {
        Name("CpaStitchDistance"),
        Comment("")
      };

      fhicl::Atom<double> CpaStitchAngle {
        Name("CpaStitchAngle"),
        Comment("")
      };

      fhicl::Atom<double> CpaXDifference {
        Name("CpaXDifference"),
        Comment("")
      };

      fhicl::Table<Fiducial> FiducialCuts {
        Name("FiducialCuts"),
        Comment("")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    };

    CpaCrossCosmicIdAlg(const core::ProviderManager &manager, const Config& config);

    CpaCrossCosmicIdAlg(const core::ProviderManager &manager, const fhicl::ParameterSet& pset) :
      CpaCrossCosmicIdAlg(manager, fhicl::Table<Config>(pset, {})()) {}

    CpaCrossCosmicIdAlg();

    ~CpaCrossCosmicIdAlg();

    void reconfigure(const Config& config);

    // Calculate the time by stitching tracks across the CPA
    std::pair<double, bool> T0FromCpaStitching(recob::Track t1, std::vector<recob::Track> tracks);

    // Tag tracks as cosmics from CPA stitching t0
    bool CpaCrossCosmicId(recob::Track track, std::vector<recob::Track> tracks, art::FindManyP<recob::Hit> hitAssoc);

  private:

    double fCpaStitchDistance;
    double fCpaStitchAngle;
    double fCpaXDifference;
    double fMinX;
    double fMinY;
    double fMinZ;
    double fMaxX;
    double fMaxY;
    double fMaxZ;
    double fBeamTimeMin;
    double fBeamTimeMax;

    detinfo::DetectorPropertiesStandard const* fDetectorProperties;
    TPCGeoAlg fTpcGeo;

  };

}

#endif
