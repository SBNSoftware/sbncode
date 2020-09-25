#ifndef GEOMETRYCOSMICIDALG_H_SEEN
#define GEOMETRYCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// GeometryCosmicIdAlg.h
//
// Functions for fiducial volume cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Persistency/Common/Ptr.h" 

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcorealg/Geometry/GeometryCore.h"

// sbncode
#include "core/ProviderManager.hh"
#include "Util.h"

// c++
#include <vector>
#include <map>


namespace ana{

  class GeometryCosmicIdAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

    };

    GeometryCosmicIdAlg(const core::ProviderManager &manager, const Config& config);

    GeometryCosmicIdAlg(const core::ProviderManager &manager, const fhicl::ParameterSet& pset) :
      GeometryCosmicIdAlg(manager, fhicl::Table<Config>(pset, {})()) {}

    GeometryCosmicIdAlg();

    ~GeometryCosmicIdAlg();

    void reconfigure(const Config& config);

    // Remove any tracks in different TPC to beam activity
    bool GeometryCosmicId(recob::Track &track, std::vector<art::Ptr<recob::Hit>> &hits, std::map<geo::TPCID, bool> &tpc_flashes);

  private:

    geo::GeometryCore const* fGeometry;
  };

}

#endif
