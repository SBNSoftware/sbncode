#ifndef CRTGEOALG_H_SEEN
#define CRTGEOALG_H_SEEN


///////////////////////////////////////////////
// CRTGeoAlg.h
//
// Wrapper for some awkward CRT geometry things
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// c++
#include <vector>

// ROOT
#include "TVector3.h"
#include "TGeoManager.h"

namespace sbn{

  struct CRTSipmGeo{
    uint32_t channel;
    double x;
    double y;
    double z;
    std::string strip;
    bool null;
  };

  // CRT strip geometry struct contains dimensions and mother module
  struct CRTStripGeo{
    std::string name;
    int sensitiveVolumeID;
    double minX;
    double maxX;
    double minY;
    double maxY;
    double minZ;
    double maxZ;
    geo::Vector_t normal;
    double width;
    std::string module;
    std::pair<int, int> sipms;
    bool null;
  };

  // CRT module geometry struct contains dimensions, daughter strips and mother tagger
  struct CRTModuleGeo{
    std::string name;
    int auxDetID;
    double minX;
    double maxX;
    double minY;
    double maxY;
    double minZ;
    double maxZ;
    geo::Vector_t normal;
    size_t planeID;
    bool top;
    std::string tagger;
    std::map<std::string, CRTStripGeo> strips;
    bool null;
  };

  // CRT tagger geometry struct contains dimensions and daughter modules
  struct CRTTaggerGeo{
    std::string name;
    double minX;
    double maxX;
    double minY;
    double maxY;
    double minZ;
    double maxZ;
    // change to ref?
    std::map<std::string, CRTModuleGeo> modules;
    bool null;
  };


  class CRTGeoAlg {
  public:

    CRTGeoAlg(geo::GeometryCore const *geometry, geo::AuxDetGeometryCore const *auxdet_geometry);
    CRTGeoAlg();

    ~CRTGeoAlg();

    // Return the volume enclosed by the whole CRT system
    std::vector<double> CRTLimits() const;

    // Get the number of taggers in the geometry
    size_t NumTaggers() const;

    // Get the total number of modules in the geometry
    size_t NumModules() const;
    // Get the number of modules in a tagger by name
    size_t NumModules(std::string taggerName) const;
    // Get the number of modules in a tagger by index
    size_t NumModules(size_t tagger_i) const;

    // Get the total number of strips in the geometry
    size_t NumStrips() const;
    // Get the number of strips in module by name
    size_t NumStrips(std::string moduleName) const;
    // Get the number of strips in  module by global index
    size_t NumStrips(size_t module_i) const;
    // Get the number of strips in module by tagger index and local module index
    size_t NumStrips(size_t tagger_i, size_t module_i) const;

    // Get the tagger geometry object by name
    CRTTaggerGeo GetTagger(std::string taggerName) const;
    // Get the tagger geometry object by index
    CRTTaggerGeo GetTagger(size_t tagger_i) const;

    // Get the module geometry object by name
    CRTModuleGeo GetModule(std::string moduleName) const;
    // Get the module geometry object by global index
    CRTModuleGeo GetModule(size_t module_i) const;
    // Get the module geometry object by tagger index and local module index
    CRTModuleGeo GetModule(size_t tagger_i, size_t module_i) const;

    // Get the strip geometry object by name
    CRTStripGeo GetStrip(std::string stripName) const;
    // Get the strip geometry object by global index
    CRTStripGeo GetStrip(size_t strip_i) const;
    // Get the strip geometry object by global module index and local strip index
    CRTStripGeo GetStrip(size_t module_i, size_t strip_i) const;
    // Get the strip geometry object by tagger index, local module index and local strip index
    CRTStripGeo GetStrip(size_t tagger_i, size_t module_i, size_t strip_i) const;

    // Get tagger name from strip or module name
    std::string GetTaggerName(std::string name) const;

    // Get the name of the strip from the SiPM channel ID
    std::string ChannelToStripName(size_t channel) const;

    // Get the world position of Sipm from the channel ID
    geo::Point_t ChannelToSipmPosition(size_t channel) const;
    
    // Get the sipm channels on a strip
    std::pair<int, int> GetStripSipmChannels(std::string stripName) const;

    // Recalculate strip limits including charge sharing
    std::vector<double> StripLimitsWithChargeSharing(std::string stripName, double x, double ex);

    // Return the distance to a sipm in the plane of the sipms
    double DistanceBetweenSipms(geo::Point_t position, size_t channel) const;
    // Returns max distance from sipms in strip
    double DistanceBetweenSipms(geo::Point_t position, std::string stripName) const;
    // Return the distance along the strip (from sipm end)
    double DistanceDownStrip(geo::Point_t position, std::string stripName) const;

    // Determine if a point is inside CRT volume
    bool IsInsideCRT(TVector3 point);
    bool IsInsideCRT(geo::Point_t point);
    // Determine if a point is inside a tagger by name
    bool IsInsideTagger(std::string taggerName, geo::Point_t point);
    bool IsInsideTagger(const CRTTaggerGeo& tagger, geo::Point_t point);
    // Determine if a point is inside a module by name
    bool IsInsideModule(std::string moduleName, geo::Point_t point);
    bool IsInsideModule(const CRTModuleGeo& module, geo::Point_t point);
    // Determine if a point is inside a strip by name
    bool IsInsideStrip(std::string stripName, geo::Point_t point);
    bool IsInsideStrip(const CRTStripGeo& strip, geo::Point_t point);

    // Check if two modules overlap in 2D
    bool CheckOverlap(const CRTModuleGeo& module1, const CRTModuleGeo& module2);
    // Check is a module overlaps with a perpendicual module in the same tagger
    bool HasOverlap(const CRTModuleGeo& module);
    bool StripHasOverlap(std::string stripName);
    std::vector<double> StripOverlap(std::string strip1Name, std::string strip2Name);

    // Find the average of the tagger entry and exit points of a true particle trajectory
    geo::Point_t TaggerCrossingPoint(std::string taggerName, const simb::MCParticle& particle);
    geo::Point_t TaggerCrossingPoint(const CRTTaggerGeo& tagger, const simb::MCParticle& particle);
    bool CrossesTagger(const CRTTaggerGeo& tagger, const simb::MCParticle& particle);
    // Find the average of the module entry and exit points of a true particle trajectory
    geo::Point_t ModuleCrossingPoint(std::string moduleName, const simb::MCParticle& particle);
    geo::Point_t ModuleCrossingPoint(const CRTModuleGeo& module, const simb::MCParticle& particle);
    bool CrossesModule(const CRTModuleGeo& module, const simb::MCParticle& particle);
    // Find the average of the strip entry and exit points of a true particle trajectory
    geo::Point_t StripCrossingPoint(std::string stripName, const simb::MCParticle& particle);
    geo::Point_t StripCrossingPoint(const CRTStripGeo& strip, const simb::MCParticle& particle);
    bool CrossesStrip(const CRTStripGeo& strip, const simb::MCParticle& particle);

    // Work out which strips the true particle crosses
    std::vector<std::string> CrossesStrips(const simb::MCParticle& particle);

    // Find the angle of true particle trajectory to tagger
    double AngleToTagger(std::string taggerName, const simb::MCParticle& particle);

    // Check if a particle enters the CRT volume
    bool EntersVolume(const simb::MCParticle& particle);
    // Check if a particle crosses the CRT volume
    bool CrossesVolume(const simb::MCParticle& particle);

    // Determine if a particle would be able to produce a hit in a tagger
    bool ValidCrossingPoint(std::string taggerName, const simb::MCParticle& particle);


  private:

    std::map<std::string, CRTTaggerGeo> fTaggers;
    std::map<std::string, CRTModuleGeo> fModules;
    std::map<std::string, CRTStripGeo> fStrips;
    std::map<int, CRTSipmGeo> fSipms;

    geo::GeometryCore const* fGeometryService;
    const geo::AuxDetGeometryCore* fAuxDetGeoCore;

  };

}

#endif

