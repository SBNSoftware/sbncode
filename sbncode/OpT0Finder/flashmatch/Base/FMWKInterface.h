#ifndef __OPT0FINDERFMWKINTERFACE_H__
#define __OPT0FINDERFMWKINTERFACE_H__

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#include <unordered_map>
#include <boost/functional/hash.hpp>

#if USING_LARSOFT == 0
#include "FMWKTools/ConfigManager.h"
#include "flashmatch/GeoAlgo/GeoAABox.h"
#include "FMWKTools/PhotonVoxels.h"
namespace flashmatch {
  /// Configuration object
  using Config_t = flashmatch::PSet;
}
#else
// #include "sbncode/OpT0Finder/flashmatch/GeoAlgo/GeoAABox.h"
#include "larcorealg/GeoAlgo/GeoAABox.h"
#include "OpT0FinderException.h"
#include "OpT0FinderLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
namespace flashmatch{
  /// Configuration object
  using Config_t = fhicl::ParameterSet;
}
#endif


namespace flashmatch {

  class DetectorSpecs {

  public:
    DetectorSpecs(std::string filename="specs.cfg");
    ~DetectorSpecs(){}

    inline static DetectorSpecs& GetME(std::string filename="detector_specs.cfg")
    {
      if(!_me) _me = new DetectorSpecs(filename);
      return *_me;
    }

    /// PMT XYZ position filler
    inline const geoalgo::Point_t& PMTPosition(size_t opch) { return _pmt_v.at(opch); }

    /// Detector active volume
    inline const geoalgo::AABox& ActiveVolume() const { return _bbox; }

    /// Detector active volume given cryo and tpc
    const geoalgo::AABox& ActiveVolume(int tpc, int cryo=0) const;

    /// # of PMTs
    inline size_t NOpDets() const { return _pmt_v.size(); }

    /// Drift velocity
    inline double DriftVelocity() const { return _drift_velocity; }

    /// Visibility
    float GetVisibility(double x, double y, double z, unsigned int opch) const;

    /// Visibility Reflected
    float GetVisibilityReflected(double x, double y, double z, unsigned int opch) const;

    /// Photon Library data access FIXME
    const std::vector<std::vector<float > >& GetPhotonLibraryData() const;

    /// Voxel definition
    #if USING_LARSOFT == 0
    inline const sim::PhotonVoxelDef& GetVoxelDef() const { return _voxel_def; }
    #endif

  private:
    static DetectorSpecs* _me;
    std::vector<geoalgo::Point_t> _pmt_v;
    geoalgo::AABox _bbox;
    std::unordered_map<std::pair<int, int>, geoalgo::AABox, boost::hash<std::pair<int, int>>> _bbox_map; ///< A bbox map (cryo,tpc) -> bbox
    double _drift_velocity;
    #if USING_LARSOFT == 0
    sim::PhotonVoxelDef _voxel_def;
    #endif
  };

}
#endif
