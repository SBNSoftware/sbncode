#include "CRTGeoAlg.h"

namespace sbn{

// Constructor - get values from the auxdet geometry service
CRTGeoAlg::CRTGeoAlg():
  CRTGeoAlg::CRTGeoAlg(lar::providerFrom<geo::Geometry>(), art::ServiceHandle<geo::AuxDetGeometry>()->GetProviderPtr())
{}

CRTGeoAlg::CRTGeoAlg(geo::GeometryCore const *geometry, geo::AuxDetGeometryCore const *auxdet_geometry) {
  fGeometryService = geometry;
  fAuxDetGeoCore = auxdet_geometry;

  // Keep track of which objects we've recorded
  std::vector<std::string> usedTaggers;
  std::vector<std::string> usedModules;
  std::vector<std::string> usedStrips;

  // Get the auxdets (strip arrays for some reason)
  const std::vector<geo::AuxDetGeo>& auxDets = fAuxDetGeoCore->AuxDetGeoVec();

  int ad_i = 0;
  // Loop over them
  for(auto const& auxDet : auxDets){

    int sv_i = 0;
    // Loop over the strips in the arrays
    for(size_t i = 0; i < auxDet.NSensitiveVolume(); i++){

      // Get the geometry object for the strip
      geo::AuxDetSensitiveGeo const& auxDetSensitive = auxDet.SensitiveVolume(i);
      std::set<std::string> volNames = {auxDetSensitive.TotalVolume()->GetName()};
      std::vector<std::vector<TGeoNode const*>> paths = fGeometryService->FindAllVolumePaths(volNames);

      // Build up a path that ROOT understands
      std::string path = "";
      for (size_t inode = 0; inode < paths.at(0).size(); inode++){
        path += paths.at(0).at(inode)->GetName();
        if(inode < paths.at(0).size() - 1){
          path += "/";
        }
      }

      // Access the path
      TGeoManager* manager = fGeometryService->ROOTGeoManager();
      manager->cd(path.c_str());

      // Get all of the things we're intersted in in the path
      TGeoNode* nodeStrip = manager->GetCurrentNode();
      TGeoNode* nodeModule = manager->GetMother(2);
      TGeoNode* nodeTagger = manager->GetMother(3);
      TGeoNode* nodeDet = manager->GetMother(4);

      // Fill the tagger information
      std::string taggerName = nodeTagger->GetName();
      if(std::find(usedTaggers.begin(), usedTaggers.end(), taggerName) == usedTaggers.end()){
        usedTaggers.push_back(taggerName);

        // Get the limits in local coords
        double halfWidth = ((TGeoBBox*)nodeTagger->GetVolume()->GetShape())->GetDX();
        double halfHeight = ((TGeoBBox*)nodeTagger->GetVolume()->GetShape())->GetDY();
        double halfLength = ((TGeoBBox*)nodeTagger->GetVolume()->GetShape())->GetDZ()/2;

        // Transform those limits to world coords
        double limits[3] = {halfWidth, halfHeight, halfLength};
        double limitsDet[3];
        nodeTagger->LocalToMaster(limits, limitsDet);
        double limitsWorld[3];
        nodeDet->LocalToMaster(limitsDet, limitsWorld);

        double limits2[3] = {-halfWidth, -halfHeight, -halfLength};
        double limitsDet2[3];
        nodeTagger->LocalToMaster(limits2, limitsDet2);
        double limitsWorld2[3];
        nodeDet->LocalToMaster(limitsDet2, limitsWorld2);

        // Create a tagger geometry object
        CRTTaggerGeo tagger = {};
        tagger.name = taggerName;
        tagger.minX = std::min(limitsWorld[0], limitsWorld2[0]);
        tagger.maxX = std::max(limitsWorld[0], limitsWorld2[0]);
        tagger.minY = std::min(limitsWorld[1], limitsWorld2[1]);
        tagger.maxY = std::max(limitsWorld[1], limitsWorld2[1]);
        tagger.minZ = std::min(limitsWorld[2], limitsWorld2[2]);
        tagger.maxZ = std::max(limitsWorld[2], limitsWorld2[2]);
        tagger.null = false;
        fTaggers[taggerName] = tagger;
      }

      // Fill the module information
      std::string moduleName = nodeModule->GetName();
      if(std::find(usedModules.begin(), usedModules.end(), moduleName) == usedModules.end()){
        usedModules.push_back(moduleName);

        // Technically the auxdet is the strip array but this is basically the same as the module
        // Get the limits in local coordinates
        double halfWidth = auxDet.HalfWidth1();
        double halfHeight = auxDet.HalfHeight();
        double halfLength = auxDet.Length()/2;

        // Transform to world coordinates
        geo::AuxDetGeo::LocalPoint_t const limits{halfWidth, halfHeight, halfLength};
        geo::AuxDetGeo::LocalPoint_t const limits2{-halfWidth, -halfHeight, -halfLength};
        auto const limitsWorld = auxDet.toWorldCoords(limits);
        auto const limitsWorld2 = auxDet.toWorldCoords(limits2);

        // Determine which plane the module is in in the tagger (XY configuration)
        double origin[3] = {0, 0, 0};
        double modulePosMother[3];
        nodeModule->LocalToMaster(origin, modulePosMother);
        size_t planeID = (modulePosMother[2] > 0);
        // Are the sipms at the top or bottom
        bool top = (planeID == 1) ? (modulePosMother[1] > 0) : (modulePosMother[0] < 0);

        // Create a module geometry object
        CRTModuleGeo module = {};
        module.name = moduleName;
        module.auxDetID = ad_i;
        module.minX = std::min(limitsWorld.X(), limitsWorld2.X());
        module.maxX = std::max(limitsWorld.X(), limitsWorld2.X());
        module.minY = std::min(limitsWorld.Y(), limitsWorld2.Y());
        module.maxY = std::max(limitsWorld.Y(), limitsWorld2.Y());
        module.minZ = std::min(limitsWorld.Z(), limitsWorld2.Z());
        module.maxZ = std::max(limitsWorld.Z(), limitsWorld2.Z());
        module.normal = auxDet.GetNormalVector();
        module.null = false;
        module.planeID = planeID;
        module.top = top;
        module.tagger = taggerName;
        fModules[moduleName] = module;
      }

      // Fill the strip information
      std::string stripName = nodeStrip->GetName();
      if(std::find(usedStrips.begin(), usedStrips.end(), stripName) == usedStrips.end()){
        usedStrips.push_back(stripName);

        // Get the limits in local coordinates
        double halfWidth = auxDetSensitive.HalfWidth1();
        double halfHeight = auxDetSensitive.HalfHeight();
        double halfLength = auxDetSensitive.Length()/2;

        // Transform to world coordinates
        geo::AuxDetSensitiveGeo::LocalPoint_t const limits{halfWidth, halfHeight, halfLength};
        geo::AuxDetSensitiveGeo::LocalPoint_t const limits2{-halfWidth, -halfHeight, -halfLength};
        auto const limitsWorld = auxDetSensitive.toWorldCoords(limits);
        auto const limitsWorld2 = auxDetSensitive.toWorldCoords(limits2);

        // Create a strip geometry object
        CRTStripGeo strip = {};
        strip.name = stripName;
        strip.sensitiveVolumeID = sv_i;
        strip.minX = std::min(limitsWorld.X(), limitsWorld2.X());
        strip.maxX = std::max(limitsWorld.X(), limitsWorld2.X());
        strip.minY = std::min(limitsWorld.Y(), limitsWorld2.Y());
        strip.maxY = std::max(limitsWorld.Y(), limitsWorld2.Y());
        strip.minZ = std::min(limitsWorld.Z(), limitsWorld2.Z());
        strip.maxZ = std::max(limitsWorld.Z(), limitsWorld2.Z());
        strip.normal = auxDetSensitive.GetNormalVector();
        strip.width = halfWidth * 2.;
        strip.null = false;
        strip.module = moduleName;

        // Create a sipm geometry object
        uint32_t channel0 = 32 * ad_i + 2 * sv_i + 0;
        uint32_t channel1 = 32 * ad_i + 2 * sv_i + 1;
        // Sipm0 is on the left in local coords
        double sipm0X = -halfWidth;
        double sipm1X = halfWidth;
        // In local coordinates the Y position is at half height (top if top) (bottom if not)
        double sipmY = halfHeight;
        if(!fModules[moduleName].top) sipmY = - halfHeight;
        geo::AuxDetSensitiveGeo::LocalPoint_t const sipm0XYZ{sipm0X, sipmY, 0};
        auto const sipm0XYZWorld = auxDetSensitive.toWorldCoords(sipm0XYZ);

        CRTSipmGeo sipm0 = {};
        sipm0.channel = channel0;
        sipm0.x = sipm0XYZWorld.X();
        sipm0.y = sipm0XYZWorld.Y();
        sipm0.z = sipm0XYZWorld.Z();
        sipm0.strip = stripName;
        sipm0.null = false;
        fSipms[channel0] = sipm0;

        geo::AuxDetSensitiveGeo::LocalPoint_t const sipm1XYZ{sipm1X, sipmY, 0};
        auto const sipm1XYZWorld = auxDetSensitive.toWorldCoords(sipm1XYZ);

        CRTSipmGeo sipm1 = {};
        sipm1.channel = channel1;
        sipm1.x = sipm1XYZWorld.X();
        sipm1.y = sipm1XYZWorld.Y();
        sipm1.z = sipm1XYZWorld.Z();
        sipm1.strip = stripName;
        sipm1.null = false;
        fSipms[channel1] = sipm1;

        strip.sipms = std::make_pair(channel0, channel1);
        fStrips[stripName] = strip;
        // Add the strip to the relevant module
        fModules[moduleName].strips[stripName] = strip;
      }
      sv_i++;
    }
    ad_i++;
  }

  // Need to fill the tagger modules after all the strip associations have been made
  // Loop over the modules
  for(auto const& module : fModules){
    // Fill the tagger map
    std::string taggerName = module.second.tagger;
    std::string moduleName = module.second.name;
    fTaggers[taggerName].modules[moduleName] = module.second;
  }

}


CRTGeoAlg::~CRTGeoAlg(){

}

// ----------------------------------------------------------------------------------
// Return the volume enclosed by the whole CRT system
std::vector<double> CRTGeoAlg::CRTLimits() const {
  std::vector<double> limits;

  std::vector<double> minXs;
  std::vector<double> minYs;
  std::vector<double> minZs;
  std::vector<double> maxXs;
  std::vector<double> maxYs;
  std::vector<double> maxZs;
  for(auto const& tagger : fTaggers){
    minXs.push_back(tagger.second.minX);
    minYs.push_back(tagger.second.minY);
    minZs.push_back(tagger.second.minZ);
    maxXs.push_back(tagger.second.maxX);
    maxYs.push_back(tagger.second.maxY);
    maxZs.push_back(tagger.second.maxZ);
  }
  limits.push_back(*std::min_element(minXs.begin(), minXs.end()));
  limits.push_back(*std::min_element(minYs.begin(), minYs.end()));
  limits.push_back(*std::min_element(minZs.begin(), minZs.end()));
  limits.push_back(*std::max_element(maxXs.begin(), maxXs.end()));
  limits.push_back(*std::max_element(maxYs.begin(), maxYs.end()));
  limits.push_back(*std::max_element(maxZs.begin(), maxZs.end()));

  return limits;
}

// ----------------------------------------------------------------------------------
// Get the number of taggers in the geometry
size_t CRTGeoAlg::NumTaggers() const{
  return fTaggers.size();
}


// ----------------------------------------------------------------------------------
// Get the total number of modules in the geometry
size_t CRTGeoAlg::NumModules() const{
  return fModules.size();
}

// Get the number of modules in a tagger by name
size_t CRTGeoAlg::NumModules(std::string taggerName) const{
  CRTTaggerGeo tagger = GetTagger(taggerName);
  if(!tagger.null) return tagger.modules.size();
  return 0;
}

// Get the number of modules in a tagger by index
size_t CRTGeoAlg::NumModules(size_t tagger_i) const{
  CRTTaggerGeo tagger = GetTagger(tagger_i);
  if(!tagger.null) return tagger.modules.size();
  return 0;
}


// ----------------------------------------------------------------------------------
// Get the total number of strips in the geometry
size_t CRTGeoAlg::NumStrips() const{
  return fStrips.size();
}

// Get the number of strips in module by name
size_t CRTGeoAlg::NumStrips(std::string moduleName) const{
  CRTModuleGeo module = GetModule(moduleName);
  if(!module.null) return module.strips.size();
  return 0;
}

// Get the number of strips in  module by global index
size_t CRTGeoAlg::NumStrips(size_t module_i) const{
  CRTModuleGeo module = GetModule(module_i);
  if(!module.null) return module.strips.size();
  return 0;
}

// Get the number of strips in module by tagger index and local module index
size_t CRTGeoAlg::NumStrips(size_t tagger_i, size_t module_i) const{
  CRTModuleGeo module = GetModule(tagger_i, module_i);
  if(!module.null) return module.strips.size();
  return 0;
}

// ----------------------------------------------------------------------------------
// Get the tagger geometry object by name
CRTTaggerGeo CRTGeoAlg::GetTagger(std::string taggerName) const{
  for(auto const& tagger : fTaggers){
    if(taggerName == tagger.first) return tagger.second;
  }
  CRTTaggerGeo nullTagger = {};
  nullTagger.null = true;
  return nullTagger;
}

// Get the tagger geometry object by index
CRTTaggerGeo CRTGeoAlg::GetTagger(size_t tagger_i) const{
  size_t index = 0;
  for(auto const& tagger : fTaggers){
    if(tagger_i == index) return tagger.second;
    index++;
  }
  CRTTaggerGeo nullTagger = {};
  nullTagger.null = true;
  return nullTagger;
}


// ----------------------------------------------------------------------------------
// Get the module geometry object by name
CRTModuleGeo CRTGeoAlg::GetModule(std::string moduleName) const{
  for(auto const& module : fModules){
    if(moduleName == module.first) return module.second;
  }
  CRTModuleGeo nullModule = {};
  nullModule.null = true;
  return nullModule;
}

// Get the module geometry object by global index
CRTModuleGeo CRTGeoAlg::GetModule(size_t module_i) const{
  size_t index = 0;
  for(auto const& module : fModules){
    if(module_i == index) return module.second;
    index++;
  }
  CRTModuleGeo nullModule = {};
  nullModule.null = true;
  return nullModule;
}

// Get the module geometry object by tagger index and local module index
CRTModuleGeo CRTGeoAlg::GetModule(size_t tagger_i, size_t module_i) const{
  CRTModuleGeo nullModule = {};
  nullModule.null = true;

  CRTTaggerGeo tagger = GetTagger(tagger_i);
  if(tagger.null) return nullModule;

  size_t index = 0;
  for(auto const& module : tagger.modules){
    if(module_i == index) return module.second;
    index++;
  }
  return nullModule;
}


// ----------------------------------------------------------------------------------
// Get the strip geometry object by name
CRTStripGeo CRTGeoAlg::GetStrip(std::string stripName) const{
  for(auto const& strip : fStrips){
    if(stripName == strip.first) return strip.second;
  }
  CRTStripGeo nullStrip = {};
  nullStrip.null = true;
  return nullStrip;
}

// Get the strip geometry object by global index
CRTStripGeo CRTGeoAlg::GetStrip(size_t strip_i) const{
  size_t index = 0;
  for(auto const& strip : fStrips){
    if(strip_i == index) return strip.second;
    index++;
  }
  CRTStripGeo nullStrip = {};
  nullStrip.null = true;
  return nullStrip;
}

// Get the strip geometry object by global module index and local strip index
CRTStripGeo CRTGeoAlg::GetStrip(size_t module_i, size_t strip_i) const{
  CRTStripGeo nullStrip = {};
  nullStrip.null = true;

  CRTModuleGeo module = GetModule(module_i);
  if(module.null) return nullStrip;

  size_t index = 0;
  for(auto const& strip : module.strips){
    if(strip_i == index) return strip.second;
    index++;
  }
  return nullStrip;
}

// Get the strip geometry object by tagger index, local module index and local strip index
CRTStripGeo CRTGeoAlg::GetStrip(size_t tagger_i, size_t module_i, size_t strip_i) const{
  CRTStripGeo nullStrip = {};
  nullStrip.null = true;

  CRTModuleGeo module = GetModule(tagger_i, module_i);
  if(module.null) return nullStrip;

  size_t index = 0;
  for(auto const& strip : module.strips){
    if(strip_i == index) return strip.second;
    index++;
  }
  return nullStrip;
}


// Get the tagger name from strip or module name
std::string CRTGeoAlg::GetTaggerName(std::string name) const{
  if(fModules.find(name) != fModules.end()){
    return fModules.at(name).tagger;
  }
  if(fStrips.find(name) != fStrips.end()){
    return fModules.at(fStrips.at(name).module).tagger;
  }
  return "";
}

// Get the name of the strip from the SiPM channel ID
std::string CRTGeoAlg::ChannelToStripName(size_t channel) const{
  for(auto const& sipm : fSipms){
    if(sipm.second.channel == channel){
      return sipm.second.strip;
    }
  }
  return "";
}


// Recalculate strip limits including charge sharing
std::vector<double> CRTGeoAlg::StripLimitsWithChargeSharing(std::string stripName, double x, double ex){
  int module = fModules.at(fStrips.at(stripName).module).auxDetID;
  std::string moduleName = fAuxDetGeoCore->AuxDet(module).TotalVolume()->GetName();
  auto const& sensitiveGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(moduleName,
                                                                      2*fStrips.at(stripName).sensitiveVolumeID);

  double halfWidth = sensitiveGeo.HalfWidth1();
  double halfHeight = sensitiveGeo.HalfHeight();
  double halfLength = sensitiveGeo.HalfLength();

  // Get the maximum strip limits in world coordinates
  geo::AuxDetSensitiveGeo::LocalPoint_t const l1{-halfWidth + x + ex, halfHeight, halfLength};
  geo::AuxDetSensitiveGeo::LocalPoint_t const l2{-halfWidth + x - ex, -halfHeight, -halfLength};
  auto const w1 = sensitiveGeo.toWorldCoords(l1);
  auto const w2 = sensitiveGeo.toWorldCoords(l2);

  // Use this to get the limits in the two variable directions
  std::vector<double> limits = {std::min(w1.X(),w2.X()), std::max(w1.X(),w2.X()),
                                std::min(w1.Y(),w2.Y()), std::max(w1.Y(),w2.Y()),
                                std::min(w1.Z(),w2.Z()), std::max(w1.Z(),w2.Z())};
  return limits;
}

// Get the world position of Sipm from the channel ID
geo::Point_t CRTGeoAlg::ChannelToSipmPosition(size_t channel) const{
  for(auto const& sipm : fSipms){
    if(sipm.second.channel == channel){
      geo::Point_t position {sipm.second.x, sipm.second.y, sipm.second.z};
      return position;
    }
  }
  geo::Point_t null {-99999, -99999, -99999};
  return null;
}

// Get the sipm channels on a strip
std::pair<int, int> CRTGeoAlg::GetStripSipmChannels(std::string stripName) const{
  for(auto const& strip : fStrips){
    if(stripName == strip.first) return strip.second.sipms;
  }
  return std::make_pair(-99999, -99999);
}

// Return the distance to a sipm in the plane of the sipms
double CRTGeoAlg::DistanceBetweenSipms(geo::Point_t position, size_t channel) const{
  double distance = -99999;

  for(auto const& sipm : fSipms){
    if(sipm.second.channel == channel){
      geo::Point_t pos {sipm.second.x, sipm.second.y, sipm.second.z};
      // Get the other sipm
      int otherChannel = channel + 1;
      if(channel % 2) otherChannel = channel - 1;
      // Work out which coordinate is different
      if(fSipms.at(otherChannel).x != pos.X()) distance = position.X() - pos.X();
      if(fSipms.at(otherChannel).y != pos.Y()) distance = position.Y() - pos.Y();
      if(fSipms.at(otherChannel).z != pos.Z()) distance = position.Z() - pos.Z();
      // Return distance in that coordinate
      return distance;
    }
  }
  return distance;
}

// Returns max distance from sipms in strip
double CRTGeoAlg::DistanceBetweenSipms(geo::Point_t position, std::string stripName) const{
  std::pair<int, int> sipms = GetStripSipmChannels(stripName);
  double sipmDist = std::max(DistanceBetweenSipms(position, sipms.first), DistanceBetweenSipms(position, sipms.second));
  return sipmDist;
}

// Return the distance along the strip (from sipm end)
double CRTGeoAlg::DistanceDownStrip(geo::Point_t position, std::string stripName) const{
  double distance = -99999;
  for(auto const& strip : fStrips){
    if(stripName == strip.first){
      geo::Point_t pos = ChannelToSipmPosition(strip.second.sipms.first);
      // Work out the longest dimension of strip
      double xdiff = std::abs(strip.second.maxX-strip.second.minX);
      double ydiff = std::abs(strip.second.maxY-strip.second.minY);
      double zdiff = std::abs(strip.second.maxZ-strip.second.minZ);
      if(xdiff > ydiff && xdiff > zdiff) distance = position.X() - pos.X();
      if(ydiff > xdiff && ydiff > zdiff) distance = position.Y() - pos.Y();
      if(zdiff > xdiff && zdiff > ydiff) distance = position.Z() - pos.Z();
      return std::abs(distance);
    }
  }
  return distance;
}

// ----------------------------------------------------------------------------------
// Determine if a point is inside CRT volume
bool CRTGeoAlg::IsInsideCRT(TVector3 point){
  geo::Point_t pt {point.X(), point.Y(), point.Z()};
  return IsInsideCRT(pt);
}

bool CRTGeoAlg::IsInsideCRT(geo::Point_t point){
  std::vector<double> limits = CRTLimits();
  if(point.X() > limits[0] && point.Y() > limits[1] && point.Z() > limits[2]
     && point.X() < limits[3] && point.Y() < limits[4] && point.Z() < limits[5]){
    return true;
  }
  return false;
}

// ----------------------------------------------------------------------------------
// Determine if a point is inside a tagger by name
bool CRTGeoAlg::IsInsideTagger(std::string taggerName, geo::Point_t point){
  CRTTaggerGeo tagger = GetTagger(taggerName);
  return IsInsideTagger(tagger, point);
}

bool CRTGeoAlg::IsInsideTagger(const CRTTaggerGeo& tagger, geo::Point_t point){
  if(tagger.null) return false;
  double x = point.X();
  double y = point.Y();
  double z = point.Z();
  double xmin = tagger.minX;
  double ymin = tagger.minY;
  double zmin = tagger.minZ;
  double xmax = tagger.maxX;
  double ymax = tagger.maxY;
  double zmax = tagger.maxZ;
  if(x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax) return true;
  return false;
}

// ----------------------------------------------------------------------------------
// Determine if a point is inside a module by name
bool CRTGeoAlg::IsInsideModule(std::string moduleName, geo::Point_t point){
  CRTModuleGeo module = GetModule(moduleName);
  return IsInsideModule(module, point);
}

bool CRTGeoAlg::IsInsideModule(const CRTModuleGeo& module, geo::Point_t point){
  if(module.null) return false;
  double x = point.X();
  double y = point.Y();
  double z = point.Z();
  double xmin = module.minX;
  double ymin = module.minY;
  double zmin = module.minZ;
  double xmax = module.maxX;
  double ymax = module.maxY;
  double zmax = module.maxZ;
  // Make the width limits a bit more generous to account for steps in true trajectory
  /*if(std::abs(xmax-xmin)==1){ xmin-=1; xmax+=1; }
  if(std::abs(ymax-ymin)==1){ ymin-=1; ymax+=1; }
  if(std::abs(zmax-zmin)==1){ zmin-=1; zmax+=1; }*/
  if(x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax) return true;
  return false;
}

// ----------------------------------------------------------------------------------
// Determine if a point is inside a strip by name
bool CRTGeoAlg::IsInsideStrip(std::string stripName, geo::Point_t point){
  CRTStripGeo strip = GetStrip(stripName);
  return IsInsideStrip(strip, point);
}

bool CRTGeoAlg::IsInsideStrip(const CRTStripGeo& strip, geo::Point_t point){
  if(strip.null) return false;
  double x = point.X();
  double y = point.Y();
  double z = point.Z();
  double xmin = strip.minX;
  double ymin = strip.minY;
  double zmin = strip.minZ;
  double xmax = strip.maxX;
  double ymax = strip.maxY;
  double zmax = strip.maxZ;
  // Make the width limits a bit more generous to account for steps in true trajectory
  /*if(std::abs(xmax-xmin)==1){ xmin-=0.5; xmax+=0.5; }
  if(std::abs(ymax-ymin)==1){ ymin-=0.5; ymax+=0.5; }
  if(std::abs(zmax-zmin)==1){ zmin-=0.5; zmax+=0.5; }*/
  if(x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax) return true;
  return false;
}

// ----------------------------------------------------------------------------------
// Check if two modules overlap in 2D
bool CRTGeoAlg::CheckOverlap(const CRTModuleGeo& module1, const CRTModuleGeo& module2){
  // Get the minimum and maximum X, Y, Z coordinates
  double minX = std::max(module1.minX, module2.minX);
  double maxX = std::min(module1.maxX, module2.maxX);
  double minY = std::max(module1.minY, module2.minY);
  double maxY = std::min(module1.maxY, module2.maxY);
  double minZ = std::max(module1.minZ, module2.minZ);
  double maxZ = std::min(module1.maxZ, module2.maxZ);

  // If the two strips overlap in 2 dimensions then return true
  if ((minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ)) return true;
  // Otherwise return a "null" value
  return false;
}

// ----------------------------------------------------------------------------------
// Check is a module overlaps with a perpendicual module in the same tagger
bool CRTGeoAlg::HasOverlap(const CRTModuleGeo& module){
  // Record plane of mother module
  size_t planeID = module.planeID;
  // Get mother tagger of module
  std::string taggerName = module.tagger;
  // Loop over other modules in tagger
  for(auto const& module2 : fTaggers[taggerName].modules){
    // If in other plane loop over strips
    if(module2.second.planeID == planeID) continue;
    // Check for overlaps
    if(CheckOverlap(module, module2.second)) return true;
  }
  return false;
}

bool CRTGeoAlg::StripHasOverlap(std::string stripName){
  return HasOverlap(fModules.at(fStrips.at(stripName).module));
}

std::vector<double> CRTGeoAlg::StripOverlap(std::string strip1Name, std::string strip2Name){
  auto const& strip1 = fStrips.at(strip1Name);
  auto const& strip2 = fStrips.at(strip2Name);

  double minX = std::max(strip1.minX, strip2.minX);
  double maxX = std::min(strip1.maxX, strip2.maxY);
  double minY = std::max(strip1.minY, strip2.minY);
  double maxY = std::min(strip1.maxY, strip2.maxY);
  double minZ = std::max(strip1.minZ, strip2.minZ);
  double maxZ = std::min(strip1.maxZ, strip2.maxZ);

  std::vector<double> null = {-99999, -99999, -99999, -99999, -99999, -99999};
  std::vector<double> overlap = {minX, maxX, minY, maxY, minZ, maxZ};

  // If the two strips overlap in 2 dimensions then return the overlap
  if ((minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ)) return overlap;
  // Otherwise return a "null" value
  return null;
}

// ----------------------------------------------------------------------------------
// Find the average of the tagger entry and exit points of a true particle trajectory
geo::Point_t CRTGeoAlg::TaggerCrossingPoint(std::string taggerName, const simb::MCParticle& particle){
  const CRTTaggerGeo& tagger = fTaggers.at(taggerName);
  return TaggerCrossingPoint(tagger, particle);
}

geo::Point_t CRTGeoAlg::TaggerCrossingPoint(const CRTTaggerGeo& tagger, const simb::MCParticle& particle){
  geo::Point_t entry {-99999, -99999, -99999};
  geo::Point_t exit {-99999, -99999, -99999};

  bool first = true;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(!IsInsideTagger(tagger, point)) continue;
    if(first){
      entry = point;
      first = false;
    }
    exit = point;
  }

  geo::Point_t cross {(entry.X()+exit.X())/2., (entry.Y()+exit.Y())/2., (entry.Z()+exit.Z())/2};
  return cross;
}

bool CRTGeoAlg::CrossesTagger(const CRTTaggerGeo& tagger, const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideTagger(tagger, point)) return true;
  }
  return false;
}

// ----------------------------------------------------------------------------------
// Find the average of the module entry and exit points of a true particle trajectory
geo::Point_t CRTGeoAlg::ModuleCrossingPoint(std::string moduleName, const simb::MCParticle& particle){
  const CRTModuleGeo& module = fModules.at(moduleName);
  return ModuleCrossingPoint(module, particle);
}

geo::Point_t CRTGeoAlg::ModuleCrossingPoint(const CRTModuleGeo& module, const simb::MCParticle& particle){
  geo::Point_t entry {-99999, -99999, -99999};
  geo::Point_t exit {-99999, -99999, -99999};

  bool first = true;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(!IsInsideModule(module, point)){
      // Look at the mid point
      if(i == particle.NumberTrajectoryPoints()-1) continue;
      geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
      geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
      if(!IsInsideModule(module, mid)) continue;
      point = mid;
    }
    if(first){
      entry = point;
      first = false;
    }
    exit = point;

    if(i == particle.NumberTrajectoryPoints()-1) continue;
    geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
    geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
    if(!IsInsideModule(module, mid)) continue;
    exit = mid;
  }

  geo::Point_t cross {(entry.X()+exit.X())/2., (entry.Y()+exit.Y())/2., (entry.Z()+exit.Z())/2};
  return cross;
}

bool CRTGeoAlg::CrossesModule(const CRTModuleGeo& module, const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideModule(module, point)) return true;

    // Also look at midpoint with next point
    if(i == particle.NumberTrajectoryPoints()-1) continue;
    geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
    geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
    if(IsInsideModule(module, mid)) return true;
  }
  return false;
}

// ----------------------------------------------------------------------------------
// Find the average of the strip entry and exit points of a true particle trajectory
geo::Point_t CRTGeoAlg::StripCrossingPoint(std::string stripName, const simb::MCParticle& particle){
  const CRTStripGeo& strip = fStrips.at(stripName);
  return StripCrossingPoint(strip, particle);
}

geo::Point_t CRTGeoAlg::StripCrossingPoint(const CRTStripGeo& strip, const simb::MCParticle& particle){
  geo::Point_t entry {-99999, -99999, -99999};
  geo::Point_t exit {-99999, -99999, -99999};

  bool first = true;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(!IsInsideStrip(strip, point)){
      if(i == particle.NumberTrajectoryPoints()-1) continue;
      geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
      geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
      if(!IsInsideStrip(strip, mid)) continue;
      point = mid;
    }
    if(first){
      entry = point;
      first = false;
    }
    exit = point;
    if(i == particle.NumberTrajectoryPoints()-1) continue;
    geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
    geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
    if(!IsInsideStrip(strip, mid)) continue;
    exit = mid;
  }

  geo::Point_t cross {(entry.X()+exit.X())/2., (entry.Y()+exit.Y())/2., (entry.Z()+exit.Z())/2};
  return cross;
}

bool CRTGeoAlg::CrossesStrip(const CRTStripGeo& strip, const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideStrip(strip, point)) return true;

    // Also look at midpoint with next point
    if(i == particle.NumberTrajectoryPoints()-1) continue;
    geo::Point_t next {particle.Vx(i+1), particle.Vy(i+1), particle.Vz(i+1)};
    geo::Point_t mid {(point.X()+next.X())/2, (point.Y()+next.Y())/2, (point.Z()+next.Z())/2};
    if(IsInsideStrip(strip, mid)) return true;
  }
  return false;
}


// ----------------------------------------------------------------------------------
// Work out which strips the true particle crosses
std::vector<std::string> CRTGeoAlg::CrossesStrips(const simb::MCParticle& particle){
  std::vector<std::string> stripNames;
  for(auto const& tagger : fTaggers){
    if(!CrossesTagger(tagger.second, particle)) continue;
    for(auto const& module : tagger.second.modules){
      if(!CrossesModule(module.second, particle)) continue;
      for(auto const& strip : module.second.strips){
        if(!CrossesStrip(strip.second, particle)) continue;
        if(std::find(stripNames.begin(), stripNames.end(), strip.first) != stripNames.end()) continue;
        stripNames.push_back(strip.first);
      }
    }
  }
  return stripNames;
}


// ----------------------------------------------------------------------------------
// Find the angle of true particle trajectory to tagger
double CRTGeoAlg::AngleToTagger(std::string taggerName, const simb::MCParticle& particle){
  // Get normal to tagger using the top modules
  TVector3 normal (0,0,0);
  for(auto const& module : GetTagger(taggerName).modules){
    normal.SetXYZ(module.second.normal.X(), module.second.normal.Y(), module.second.normal.Z());
    break;
  }
  //FIXME this is pretty horrible
  if(normal.X()<0.5 && normal.X()>-0.5) normal.SetX(0);
  if(normal.Y()<0.5 && normal.Y()>-0.5) normal.SetY(0);
  if(normal.Z()<0.5 && normal.Z()>-0.5) normal.SetZ(0);

  if(std::abs(normal.X())==1 && GetTagger(taggerName).minX < 0) normal.SetX(-1);
  if(std::abs(normal.X())==1 && GetTagger(taggerName).minX > 0) normal.SetX(1);
  if(std::abs(normal.Y())==1 && GetTagger(taggerName).minY < 0) normal.SetY(-1);
  if(std::abs(normal.Y())==1 && GetTagger(taggerName).minY > 0) normal.SetY(1);
  if(std::abs(normal.Z())==1 && GetTagger(taggerName).minZ < 0) normal.SetZ(-1);
  if(std::abs(normal.Z())==1 && GetTagger(taggerName).minZ > 0) normal.SetZ(1);

  TVector3 start (particle.Vx(), particle.Vy(), particle.Vz());
  TVector3 end (particle.EndX(), particle.EndY(), particle.EndZ());
  TVector3 diff = end - start;

  if(normal.X() == 1 && start.X() > end.X()) diff = start - end;
  if(normal.X() == -1 && start.X() < end.X()) diff = start - end;
  if(normal.Y() == 1 && start.Y() > end.Y()) diff = start - end;
  if(normal.Y() == -1 && start.Y() < end.Y()) diff = start - end;
  if(normal.Z() == 1 && start.Z() > end.Z()) diff = start - end;
  if(normal.Z() == -1 && start.Z() < end.Z()) diff = start - end;

  return normal.Angle(diff);
}


// ----------------------------------------------------------------------------------
// Check if a particle enters the CRT volume
bool CRTGeoAlg::EntersVolume(const simb::MCParticle& particle){
  bool enters = false;
  bool startOutside = false;
  bool endOutside = false;
  std::vector<double> limits = CRTLimits();
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideCRT(point)){
      enters = true;
    }
    else if(i == 0) startOutside = true;
    else if(i == particle.NumberTrajectoryPoints()-1) endOutside = true;
  }
  if(enters && (startOutside || endOutside)) return true;
  return false;
}

// ----------------------------------------------------------------------------------
// Check if a particle crosses the CRT volume
bool CRTGeoAlg::CrossesVolume(const simb::MCParticle& particle){
  bool enters = false;
  bool startOutside = false;
  bool endOutside = false;
  std::vector<double> limits = CRTLimits();
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    geo::Point_t point {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    if(IsInsideCRT(point)){
      enters = true;
    }
    else if(i == 0) startOutside = true;
    else if(i == particle.NumberTrajectoryPoints()-1) endOutside = true;
  }
  if(startOutside && enters && endOutside) return true;
  return false;
}


// ----------------------------------------------------------------------------------
// Determine if a particle would be able to produce a hit in a tagger
bool CRTGeoAlg::ValidCrossingPoint(std::string taggerName, const simb::MCParticle& particle){

  // Get all the crossed strips in the tagger
  std::vector<std::string> crossedModules;
  for(auto const& module : fTaggers[taggerName].modules){
    geo::Point_t crossPoint = ModuleCrossingPoint(module.second.name, particle);
    if(crossPoint.X() != -99999) crossedModules.push_back(module.second.name);
  }

  // Check if the strip has a possible overlap, return true if not
  for(size_t i = 0; i < crossedModules.size(); i++){
    if(!HasOverlap(fModules[crossedModules[i]])) return true;
    // Check if any of the crossed strips overlap, return true if they do
    for(size_t j = i; j < crossedModules.size(); j++){
      if(CheckOverlap(fModules[crossedModules[i]], fModules[crossedModules[j]])) return true;
    }
  }
  return false;
}

}
