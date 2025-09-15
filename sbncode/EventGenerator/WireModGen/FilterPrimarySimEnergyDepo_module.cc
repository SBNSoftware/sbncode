////////////////////////////////////////////////////////////////////////
// Class:       FilterPrimarySimEnergyDepo
// Plugin Type: filter (art v3_05_01)
// File:        FilterPrimarySimEnergyDepo_module.cc
//
// Module for selecting events based on the presence of data products.
//
// Useful for processing data events, which sometimes have missing 
// data products.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include <memory>

namespace sbn {
    class FilterPrimarySimEnergyDepo;
}

class sbn::FilterPrimarySimEnergyDepo : public art::EDFilter {
public:
  explicit FilterPrimarySimEnergyDepo(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilterPrimarySimEnergyDepo(FilterPrimarySimEnergyDepo const&) = delete;
  FilterPrimarySimEnergyDepo(FilterPrimarySimEnergyDepo&&) = delete;
  FilterPrimarySimEnergyDepo& operator=(FilterPrimarySimEnergyDepo const&) = delete;
  FilterPrimarySimEnergyDepo& operator=(FilterPrimarySimEnergyDepo&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fSimEnergyDepoLabel;
  std::vector<int> fTrackIDs;
  bool fVerbose;
};


sbn::FilterPrimarySimEnergyDepo::FilterPrimarySimEnergyDepo(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fSimEnergyDepoLabel(p.get<art::InputTag>("SimEnergyDepoLabel")),
    fTrackIDs(p.get<std::vector<int>>("TrackIDs", {1})),
    fVerbose(p.get<bool>("Verbose", true))
{}

bool sbn::FilterPrimarySimEnergyDepo::filter(art::Event& e) {
  art::ValidHandle<std::vector<sim::SimEnergyDeposit>> simdepo = e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fSimEnergyDepoLabel);

  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  auto const &clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const &dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  for (const sim::SimEnergyDeposit &d: *simdepo) {
    if (fVerbose) std::cout << "depo at time: " << d.T() << "x: " << d.X() << " y: " << d.Y() << " z: " << d.Z() << " nelec: " << d.NumElectrons() << " id: " << d.TrackID() << std::endl;

    // check if has correct ID, then check if fiducial and in-time
    int match_index = -1;
    for (unsigned i = 0; i < fTrackIDs.size(); i++) {
      if (d.TrackID() == fTrackIDs[i]) {
        match_index = (int)i;
      }
    }
    if (match_index == -1) continue;

    for (auto const &cryo: geometry->Iterate<geo::CryostatGeo>()) {
      for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryo.ID())) {
        if (TPC.ContainsPosition(d.MidPoint())) {
          // check drifted point in time
  
          geo::PlaneID plane(TPC.ID(), 0); // check against front plane

          // get the tick of the point. Subtract out the trigger time, which gets put into t_tick below
          double tick = dprop.ConvertXToTicks(d.X(), plane);

          // convert the gen time to a TDC
          double t_tdc = clock_data.TPCTick2TDC(clock_data.Time2Tick(d.T()/1e3));

          // add in the offset from the drift
          double x_tick = t_tdc + tick;

          if (fVerbose) std::cout << "Gen time: " << d.T() << " is at TPC tdc: " << t_tdc << ". After drift from X=" << d.X() << " arrival tick is: " << x_tick << ". TPC window is 0 to " << dprop.NumberTimeSamples() << std::endl;

          if ((x_tick > 0) && (x_tick < dprop.NumberTimeSamples())) {
            if (fVerbose) std::cout << "Found matching deposit for track ID: " << fTrackIDs[match_index] << std::endl;
            return true;
          }

        }
      }
    }
  }

  if (fVerbose) std::cout << "No matching deposit." << std::endl;

  return false;
}

DEFINE_ART_MODULE(sbn::FilterPrimarySimEnergyDepo)
