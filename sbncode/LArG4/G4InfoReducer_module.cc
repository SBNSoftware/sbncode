////////////////////////////////////////////////////////////////////////
// Class:       G4InfoReducer
// Plugin Type: producer (Unknown Unknown)
// File:        G4InfoReducer_module.cc
//
// Generated at Thu Jul 14 11:04:27 2022 by Laura Domine using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimEnergyDepositLite.h"
#include "larcore/CoreUtils/ServiceUtil.h" // for lar::providerFrom
#include "lardata/DetectorInfoServices/DetectorPropertiesServiceStandard.h" // for DetectorClocksService

#include <memory>
#include <set>

class G4InfoReducer;


class G4InfoReducer : public art::EDProducer {
public:
  explicit G4InfoReducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  G4InfoReducer(G4InfoReducer const&) = delete;
  G4InfoReducer(G4InfoReducer&&) = delete;
  G4InfoReducer& operator=(G4InfoReducer const&) = delete;
  G4InfoReducer& operator=(G4InfoReducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fSedLabel; ///< module making the SimEnergyDeposit
  double fMinX, fMinY, fMinZ; ///< bottom left coordinate of union of all TPC active volumes
  double fVoxelSizeX, fVoxelSizeY, fVoxelSizeZ; ///< size of a voxel (cm)
  bool fUseOrigTrackID; //Use orig track ID boolean
  //services
  const geo::GeometryCore& fGeometry;
};


G4InfoReducer::G4InfoReducer(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  , fGeometry(*lar::providerFrom<geo::Geometry>())
{
  fSedLabel = p.get<art::InputTag>("SimEnergyDepositLabel", "largeant:TPCActive");

  // Prepare for voxelization
  fVoxelSizeX = p.get<double>("VoxelSizeX", 0.3);
  fVoxelSizeY = p.get<double>("VoxelSizeY", 0.3);
  fVoxelSizeZ = p.get<double>("VoxelSizeZ", 0.3);

  //Use orig track id
  fUseOrigTrackID = p.get<bool>("useOrigTrackID",true);

  if (fVoxelSizeX <= 0. || fVoxelSizeY <= 0. || fVoxelSizeZ <= 0.) {
    std::cerr << "Voxel size must be strictly greater than zero." << std::endl;
    throw std::exception();
  }

  double min_x = std::numeric_limits<double>::max();
  double min_y = std::numeric_limits<double>::max();
  double min_z = std::numeric_limits<double>::max();
  for(geo::TPCGeo const& tpc: fGeometry.Iterate<geo::TPCGeo>()) {
    auto const& tpcabox = tpc.ActiveBoundingBox();
    min_x = std::min(min_x, tpcabox.MinX());
    min_y = std::min(min_y, tpcabox.MinY());
    min_z = std::min(min_z, tpcabox.MinZ());
  }
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  // Take into account TPC readout window size
  // Ionization electrons travel at ~1.6mm/us at 500 V/cm electric field
  fMinX = min_x + clockData.TriggerOffsetTPC() * 1.6 / 10;
  fMinY = min_y;
  fMinZ = min_z;
  //std::cout << "G4InfoReducer " << fMinX << " " << fMinY << " " << fMinZ << std::endl;
  //std::cout << clockData.TriggerOffsetTPC() << std::endl;

  // Call appropriate produces<>() functions here.
  produces<std::vector<sim::SimEnergyDepositLite>>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
  consumes<std::vector<sim::SimEnergyDeposit>>(fSedLabel);
}

void G4InfoReducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  // Collect SimEnergyDeposit
  auto handle = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fSedLabel);
  if (!handle.isValid()) {
    // throw some error
    std::cerr << "SimEnergyDeposit not found" << std::endl;
    throw std::exception();
  }

  struct comp {
    bool operator() (sim::SimEnergyDepositLite lhs, sim::SimEnergyDepositLite rhs) const {
      if (lhs.X() < rhs.X()) return true;
      if (lhs.X() > rhs.X()) return false;
      if (lhs.Y() < rhs.Y()) return true;
      if (lhs.Y() > rhs.Y()) return false;
      if (lhs.Z() < rhs.Z()) return true;
      if (lhs.Z() > rhs.Z()) return false;
      return lhs.TrackID() < rhs.TrackID();
    }
  };
  std::set<sim::SimEnergyDepositLite, comp> sedlite_v2;
  sim::SimEnergyDepositLite sed_lite;
  sim::SimEnergyDepositLite new_sed_lite;


  auto const& sed_v = *handle;

  /*double total_edep = 0.;
  for (size_t idx = 0; idx < sed_v.size(); ++idx) {
    total_edep += sed_v[idx].E();
  }
  std::cout << "total edep = " << total_edep << " count = " << sed_v.size() << std::endl;*/

  for (size_t idx = 0; idx < sed_v.size(); ++idx) {
    auto const& sed = sed_v[idx];
    std::cout << "G4info - orig track ID : " << sed.OrigTrackID() 
    << " trackID : " << sed.TrackID()
    << " useOrigID : " << fUseOrigTrackID
    << std::endl;
    // Voxelize coordinates to closest coordinate using fVoxelSize
    double x = sed.X() - std::fmod(sed.X() - fMinX, fVoxelSizeX);
    double y = sed.Y() - std::fmod(sed.Y() - fMinY, fVoxelSizeY);
    double z = sed.Z() - std::fmod(sed.Z() - fMinZ, fVoxelSizeZ);

    // Copy info to SimEnergyDepositLite
    if (fUseOrigTrackID){
      sed_lite = sim::SimEnergyDepositLite(sed.E(), geo::Point_t(x, y, z), sed.T(), sed.OrigTrackID());
    }
    else{
      sed_lite = sim::SimEnergyDepositLite(sed.E(), geo::Point_t(x, y, z), sed.T(), sed.TrackID());
    }
    auto it = sedlite_v2.find(sed_lite);
    // Attempt to insert, if already exist we sum up energy deposit
    
    if (it!=sedlite_v2.end()) {
      double new_energy = sed_lite.E() + it->E();
      double new_time = std::min(sed_lite.T(), it->T());
      sedlite_v2.erase(it);
      if (fUseOrigTrackID){
        new_sed_lite = sim::SimEnergyDepositLite(new_energy, geo::Point_t(x, y, z), new_time, sed.OrigTrackID());
      }
      else{
        new_sed_lite = sim::SimEnergyDepositLite(new_energy, geo::Point_t(x, y, z), new_time, sed.TrackID());
      }
      sedlite_v2.insert(new_sed_lite);
      
    }
    else {
      sedlite_v2.insert(std::move(sed_lite));
    }
  }

  /*double new_total_edep = 0.;
  int counts = 0;

  for (auto it : sedlite_v2) {
    new_total_edep += it.E();
    counts += 1;
  }
  std::cout << "new total edep = " << new_total_edep << " with counts " << counts << " " << sedlite_v2.size() << std::endl;
  */

  // Create vector for SEDLite
  std::unique_ptr<std::vector<sim::SimEnergyDepositLite>> sedlite_v(new std::vector<sim::SimEnergyDepositLite>(sedlite_v2.begin(), sedlite_v2.end()));

  // Store SEDLite in event
  e.put(std::move(sedlite_v));

}

DEFINE_ART_MODULE(G4InfoReducer)
