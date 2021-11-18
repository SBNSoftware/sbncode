
#include "TPCGeoAlg.h"

namespace sbn{

// Constructor - get values from the geometry service
TPCGeoAlg::TPCGeoAlg(){

  fMinX = 99999;
  fMinY = 99999;
  fMinZ = 99999;
  fMaxX = -99999;
  fMaxY = -99999;
  fMaxZ = -99999;
  fCpaWidth = 0;

  fGeometryService = lar::providerFrom<geo::Geometry>();

  for(size_t cryo_i = 0; cryo_i < fGeometryService->Ncryostats(); cryo_i++){
    const geo::CryostatGeo& cryostat = fGeometryService->Cryostat(cryo_i);

    for (size_t tpc_i = 0; tpc_i < cryostat.NTPC(); tpc_i++)
    {
      const geo::TPCGeo& tpcg = cryostat.TPC(tpc_i);
      if (tpcg.MinX() < fMinX) fMinX = tpcg.MinX();
      if (tpcg.MaxX() > fMaxX) fMaxX = tpcg.MaxX();
      if (tpcg.MinY() < fMinY) fMinY = tpcg.MinY();
      if (tpcg.MaxY() > fMaxY) fMaxY = tpcg.MaxY();
      if (tpcg.MinZ() < fMinZ) fMinZ = tpcg.MinZ();
      if (tpcg.MaxZ() > fMaxZ) fMaxZ = tpcg.MaxZ();
      fCpaWidth = std::min(std::abs(tpcg.MinX()), std::abs(tpcg.MaxX()));
    }
  }
}


TPCGeoAlg::~TPCGeoAlg(){

}

// ----------------------------------------------------------------------------------
// Getters
double TPCGeoAlg::MinX() const{
  return fMinX;
}

double TPCGeoAlg::MinY() const{
  return fMinY;
}

double TPCGeoAlg::MinZ() const{
  return fMinZ;
}

double TPCGeoAlg::MaxX() const{
  return fMaxX;
}

double TPCGeoAlg::MaxY() const{
  return fMaxY;
}

double TPCGeoAlg::MaxZ() const{
  return fMaxZ;
}

double TPCGeoAlg::CpaWidth() const{
  return fCpaWidth;
}

// ----------------------------------------------------------------------------------
// Functions for applying fiducial volume cuts to total volume
bool TPCGeoAlg::InFiducial(geo::Point_t point, double fiducial){
  return InFiducial(point, fiducial, fiducial);
}

bool TPCGeoAlg::InFiducial(geo::Point_t point, double fiducial, double fiducialTop){
  return InFiducial(point, fiducial, fiducial, fiducial, fiducial, fiducialTop, fiducial);
}

bool TPCGeoAlg::InFiducial(geo::Point_t point, double minXCut, double minYCut, double minZCut, 
                           double maxXCut, double maxYCut, double maxZCut){
  
  double xmin = fMinX + minXCut;
  double xmax = fMaxX - maxXCut;
  double ymin = fMinY + minYCut;
  double ymax = fMaxY - maxYCut;
  double zmin = fMinZ + minZCut;
  double zmax = fMaxZ - maxZCut;

  double x = point.X();
  double y = point.Y();
  double z = point.Z();
  if(x>xmin && x<xmax && y>ymin && y<ymax && z>zmin && z<zmax) return true;

  return false;
}

// ----------------------------------------------------------------------------------
// Determine which TPC a collection of hits is detected in (-1 if multiple) 
int TPCGeoAlg::DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits){
  // Return tpc of hit collection or -1 if in multiple
  if(hits.size() == 0) return -1;
  int tpc = hits[0]->WireID().TPC;
  for(size_t i = 0; i < hits.size(); i++){
    if((int)hits[i]->WireID().TPC != tpc) return -1;
  }
  return tpc;
}

// Determine the drift direction for a collection of hits (-1, 0 or 1 assuming drift in X)
int TPCGeoAlg::DriftDirectionFromHits(std::vector<art::Ptr<recob::Hit>> hits){
  // If there are no hits then return 0
  if(hits.size() == 0) return 0;
  
  // If the track is stitched (in multiple TPCs) return 0
  if(DetectedInTPC(hits) == -1) return 0;

  // Work out the drift direction
  geo::TPCID tpcID = hits[0]->WireID().asTPCID();
  const geo::TPCGeo& tpcGeo = fGeometryService->GetElement(tpcID);
  double driftDirection = tpcGeo.DetectDriftDirection();
  if(std::abs(driftDirection) != 1) driftDirection = 0;
  return driftDirection;
}

// Work out the drift limits for a collection of hits
std::pair<double, double> TPCGeoAlg::XLimitsFromHits(std::vector<art::Ptr<recob::Hit>> hits){
  // If there are no hits then return 0
  if(hits.size() == 0) return std::make_pair(0, 0);
  
  // If the track is stitched (in multiple TPCs) return 0
  if(DetectedInTPC(hits) == -1) return std::make_pair(0, 0);

  // Work out the drift direction
  geo::TPCID tpcID = hits[0]->WireID().asTPCID();
  const geo::TPCGeo& tpcGeo = fGeometryService->GetElement(tpcID);
  return std::make_pair(tpcGeo.MinX(), tpcGeo.MaxX());
}

// Is point inside given TPC
bool TPCGeoAlg::InsideTPC(geo::Point_t point, const geo::TPCGeo& tpc, double buffer){
  if(point.X() < (tpc.MinX()-buffer) || point.X() > (tpc.MaxX()+buffer)
      || point.Y() < (tpc.MinY()-buffer) || point.Y() > (tpc.MaxY()+buffer)
      || point.Z() < (tpc.MinZ()-buffer) || point.Z() > (tpc.MaxZ()+buffer)) return false;
  return true;
}

// Minimum distance to a TPC wall
double TPCGeoAlg::MinDistToWall(geo::Point_t point){

  std::vector<double> dists;
  
  dists.push_back(std::abs(point.X() - fMinX));
  dists.push_back(std::abs(point.X() - fMaxX));
  dists.push_back(std::abs(point.Y() - fMinY));
  dists.push_back(std::abs(point.Y() - fMaxY));
  dists.push_back(std::abs(point.Z() - fMinZ));
  dists.push_back(std::abs(point.Z() - fMaxZ));

  std::sort(dists.begin(), dists.end());

  return dists[0];

}

// ----------------------------------------------------------------------------------
// Determine if a true particle is ever inside the TPC volume
bool TPCGeoAlg::InVolume(const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    if(x > fMinX && y > fMinY && z > fMinZ && x < fMaxX && y < fMaxY && z < fMaxZ){
      return true;
    }
  }
  return false;
}

// ----------------------------------------------------------------------------------
// Determine if a true particle is contained inside the TPC volume
bool TPCGeoAlg::IsContained(const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    if(x < fMinX || y < fMinY || z < fMinZ || x > fMaxX || y > fMaxY || z > fMaxZ){
      return false;
    }
  }
  return true;
}

// ----------------------------------------------------------------------------------
// Determine if a true particle enters the TPC volume
bool TPCGeoAlg::EntersVolume(const simb::MCParticle& particle){
  bool enters = false;
  bool startOutside = false;
  bool endOutside = false;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    if(x > fMinX && y > fMinY && z > fMinZ && x < fMaxX && y < fMaxY && z < fMaxZ){
      enters = true;
    }
    else if(i == 0) startOutside = true;
    else if(i == particle.NumberTrajectoryPoints()-1) endOutside = true;
  }
  if(enters && (startOutside || endOutside)) return true;
  return false;
}

// ----------------------------------------------------------------------------------
// Determine if a true particle crosses the TPC volume
bool TPCGeoAlg::CrossesVolume(const simb::MCParticle& particle){
  bool enters = false;
  bool startOutside = false;
  bool endOutside = false;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    if(x > fMinX && y > fMinY && z > fMinZ && x < fMaxX && y < fMaxY && z < fMaxZ){
      enters = true;
    }
    else if(i == 0) startOutside = true;
    else if(i == particle.NumberTrajectoryPoints()-1) endOutside = true;
  }
  if(startOutside && enters && endOutside) return true;
  return false;
}

// ----------------------------------------------------------------------------------
// Determine if a true particle crosses either APA
bool TPCGeoAlg::CrossesApa(const simb::MCParticle& particle){
  for(size_t i = 0; i < particle.NumberTrajectoryPoints()-1; i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    double x1 = particle.Vx(i+1); 
    double y1 = particle.Vy(i+1);
    double z1 = particle.Vz(i+1);
    if(y >= fMinY && z >= fMinZ && y <= fMaxY && z <= fMaxZ
       && y1 >= fMinY && z1 >= fMinZ && y1 <= fMaxY && z1 <= fMaxZ){
      if(x <= fMinX && x1 >= fMinX) return true;
      if(x >= fMinX && x1 <= fMinX) return true;
      if(x <= fMaxX && x1 >= fMaxX) return true;
      if(x >= fMaxX && x1 <= fMaxX) return true;
    }
  }
  return false;
}

std::pair<TVector3, TVector3> TPCGeoAlg::CrossingPoints(const simb::MCParticle& particle){
  bool first = true;
  TVector3 start (-99999, -99999, -99999);
  TVector3 end (-99999, -99999, -99999);
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    if(x > fMinX && y > fMinY && z > fMinZ && x < fMaxX && y < fMaxY && z < fMaxZ){
      if(first){
        first = false;
        start.SetXYZ(x, y, z);
      }
      end.SetXYZ(x, y, z);
    }
  }
  return std::make_pair(start, end);
}

double TPCGeoAlg::TpcLength(const simb::MCParticle& particle){
  bool first = true;
  double length = 0;
  TVector3 point, disp;
  for(size_t i = 0; i < particle.NumberTrajectoryPoints(); i++){
    double x = particle.Vx(i); 
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    if(x > fMinX && y > fMinY && z > fMinZ && x < fMaxX && y < fMaxY && z < fMaxZ){
      point.SetXYZ(x, y, z);
      if(!first){
        disp -= point;
        length += disp.Mag();
      }
      first = false;
      disp = point;
    }
  }
  return length;
}


}


