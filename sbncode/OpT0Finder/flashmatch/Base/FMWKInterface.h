#ifndef __OPT0FINDERFMWKINTERFACE_H__
#define __OPT0FINDERFMWKINTERFACE_H__

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#include "OpT0FinderException.h"
#include "LoggerFeature.h"

#if USING_LARSOFT == 0
#include "FMWKTools/ConfigManager.h"
#include "FMWKTools/PhotonVoxels.h"
#include "flashmatch/GeoAlgo/GeoAABox.h"
#include "flashmatch/Base/FMWKTools/PSetUtils.h"
#include "flashmatch/Base/FMWKTools/PhotonVisibilityService.h"
namespace flashmatch {
  /// Configuration object
  using Config_t = flashmatch::PSet;
}
#else
#include "lardataobj/Utilities/LazyVector.h"
#include "sbncode/OpT0Finder/flashmatch/GeoAlgo/GeoAABox.h"
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

  class DetectorSpecs : public LoggerFeature {

  public:
    DetectorSpecs(const Config_t& cfg);
    ~DetectorSpecs(){}

    inline static DetectorSpecs& GetME()
    {
      if(_me) return *_me;
      std::cerr << "DetectorSpecs::GetME() without an argument must be called only after it is called with Config_t argument!" 
		<< std::endl << std::endl << std::flush;
      throw OpT0FinderException();
    }

    inline static DetectorSpecs& GetME(const Config_t& cfg) 
    {
      if(!_me) _me = new DetectorSpecs(cfg);
      return *_me;
    }

    /// Info dump
    void DumpInfo() const;

    /// PMT XYZ position filler
    inline const geoalgo::Point_t& PMTPosition(size_t opch) { return _pmt_v.at(opch); }

    /// Detector active volume
    inline const geoalgo::AABox& ActiveVolume() const { return _bbox; }

    /// Photon library volume
    inline const geoalgo::AABox& PhotonLibraryVolume() const { return _photon_bbox; } // TODO find better name

    /// Detector active volume given cryo and tpc
    const geoalgo::AABox& ActiveVolume(int tpc, int cryo=0) const;

    /// # of PMTs
    inline size_t NOpDets() const { return _pmt_v.size(); }

    /// Drift velocity
    inline double DriftVelocity() const { return _drift_velocity; }

    /// Light yield
    inline double LightYield() const { return _light_yield; }

    /// MIP dE/dx
    inline double MIPdEdx() const { return _MIPdEdx; }
    
    /// Visibility
    float GetVisibility(double x, double y, double z, unsigned int opch) const;

    /// Visibility Reflected
    float GetVisibilityReflected(double x, double y, double z, unsigned int opch) const;

    float GetVisibility(int vox_id, unsigned int opch) const;

    float GetVisibilityReflected(int vox_id, unsigned int opch) const;

    #if USING_LARSOFT == 0
    /// Photon Library data access
    const std::vector<float>& GetLibraryEntries(int vox_id) const;
    /// For non-larsoft option, configure via filename
    inline static DetectorSpecs& GetME(std::string filename)
    {
      assert(!filename.empty());
      if(filename.find("/") != 0)
        filename = std::string(getenv("FMATCH_DATADIR")) + "/" + filename;

      auto cfg = CreatePSetFromFile(filename,"cfg");
      auto const& p = cfg.get<::flashmatch::Config_t>("DetectorSpecs");

      if(!_me) _me = new DetectorSpecs(p);
      return *_me;
    }
    #else
    /// Photon Library data access
    phot::IPhotonLibrary::Counts_t GetLibraryEntries(int vox_id, bool reflWanted=false) const;
    /// Set which cryostats to use
    void EnableCryostats(std::vector<size_t> cryo_id_v);
    #endif
    
    /// Voxel definition
    const sim::PhotonVoxelDef& GetVoxelDef() const;

  private:
    static DetectorSpecs* _me;
    std::vector<geoalgo::Point_t> _pmt_v;
    geoalgo::AABox _bbox, _photon_bbox;
    std::map<std::pair<int, int>, geoalgo::AABox> _bbox_map; ///< A bbox map (cryo,tpc) -> bbox
    double _drift_velocity;
    double _light_yield;
    double _MIPdEdx;
    std::vector<int> _cryo_id_v;
  };

}
#endif
