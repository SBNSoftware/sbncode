#ifndef APACROSSCOSMICIDALG_H_SEEN
#define APACROSSCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// ApaCrossCosmicIdAlg.h
//
// Functions for APA crossing cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// sbncode
#include "core/ProviderManager.hh"
#include "Util.h"

// framework
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Persistency/Common/Ptr.h" 

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "larcorealg/Geometry/GeometryCore.h"

// c++
#include <vector>
#include <utility>
#include <map>

namespace ana{

  class ApaCrossCosmicIdAlg {
  public:

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

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("")
      };

      fhicl::Atom<double> MaxApaDistance {
        Name("MaxApaDistance"),
        Comment("")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    };

    ApaCrossCosmicIdAlg(const core::ProviderManager &manager, const Config& config);

    ApaCrossCosmicIdAlg(const core::ProviderManager &manager, const fhicl::ParameterSet& pset) :
      ApaCrossCosmicIdAlg(manager, fhicl::Table<Config>(pset, {})()) {}

    ApaCrossCosmicIdAlg();

    ~ApaCrossCosmicIdAlg();

    void reconfigure(const core::ProviderManager &manager, const Config& config);

    // Get the minimum distance from track to APA for different times
    std::pair<double, double> MinApaDistance(const recob::Track &track, std::vector<double> &t0List, geo::TPCID &tpcid);

    // Get time by matching tracks which cross the APA
    double T0FromApaCross(const recob::Track &track, std::vector<art::Ptr<recob::Hit>> hits, std::map<geo::CryostatID, std::vector<double>> &t_zeros);

    // Get the distance from track to APA at fixed time
    double ApaDistance(recob::Track track, double t0, std::vector<art::Ptr<recob::Hit>> hits);

    // Tag tracks with times outside the beam
    bool ApaCrossCosmicId(const recob::Track &track, std::vector<art::Ptr<recob::Hit>> &hits, std::map<geo::CryostatID, std::vector<double>> &t_zeros);

  private:

    double fDistanceLimit;
    double fMaxApaDistance;
    double fBeamTimeMin;
    double fBeamTimeMax;

    geo::GeometryCore const* fGeometry;
    detinfo::DetectorPropertiesStandard const* fDetectorProperties;
  };

}

#endif
