////////////////////////////////////////////////////////////////////////
// Class:       FilterPrimarySimChannel
// Plugin Type: filter (art v3_05_01)
// File:        FilterPrimarySimChannel_module.cc
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
#include "lardataobj/Simulation/SimChannel.h"

#include <memory>

namespace sbn {
    class FilterPrimarySimChannel;
}

class sbn::FilterPrimarySimChannel : public art::EDFilter {
public:
  explicit FilterPrimarySimChannel(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilterPrimarySimChannel(FilterPrimarySimChannel const&) = delete;
  FilterPrimarySimChannel(FilterPrimarySimChannel&&) = delete;
  FilterPrimarySimChannel& operator=(FilterPrimarySimChannel const&) = delete;
  FilterPrimarySimChannel& operator=(FilterPrimarySimChannel&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fSimChannelLabel;
  std::vector<int> fTrackIDs;
  bool fVerbose;
};


sbn::FilterPrimarySimChannel::FilterPrimarySimChannel(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fSimChannelLabel(p.get<art::InputTag>("SimChannelLabel")),
    fTrackIDs(p.get<std::vector<int>>("TrackIDs", {1})),
    fVerbose(p.get<bool>("Verbose", true))
{}

bool sbn::FilterPrimarySimChannel::filter(art::Event& e) {
  art::ValidHandle<std::vector<sim::SimChannel>> simchan = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel);

  for (const sim::SimChannel &c: *simchan) {
    const sim::SimChannel::TDCIDEs_t &tides = c.TDCIDEMap();

    for (const sim::TDCIDE &tide: tides) {
      for (const sim::IDE &ide: tide.second) {
        if (fVerbose) std::cout << "ide at tdc: " << tide.first << " id: " << ide.trackID << " x: " << ide.x << " y: " << ide.y << " z: " << ide.z << " nelec: " << ide.numElectrons << std::endl;
        for (int tid: fTrackIDs) {
          if (ide.trackID == tid) {
            if (fVerbose) std::cout << "Found matching sim channel for track ID: " << tid << std::endl;
            return true; 
          }
        }
      }
    }
  }

  if (fVerbose) std::cout << "No matching sim channel." << std::endl;

  return false;
}

DEFINE_ART_MODULE(sbn::FilterPrimarySimChannel)
