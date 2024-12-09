////////////////////////////////////////////////////////////////////////
// Class:       TrackHitDumper
// Plugin Type: analyzer (art v3_03_01)
// File:        TrackHitDumper_module.cc
//
// Generated at Tue Apr 14 11:12:49 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <memory>

namespace sbn {
  class TrackHitDumper;
}


class sbn::TrackHitDumper : public art::EDAnalyzer {
public:
  explicit TrackHitDumper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackHitDumper(TrackHitDumper const&) = delete;
  TrackHitDumper(TrackHitDumper&&) = delete;
  TrackHitDumper& operator=(TrackHitDumper const&) = delete;
  TrackHitDumper& operator=(TrackHitDumper&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  std::vector<std::string> fHitLabels;
};

sbn::TrackHitDumper::TrackHitDumper(fhicl::ParameterSet const& p)
  : art::EDAnalyzer{p}
{
  fHitLabels = p.get<std::vector<std::string>>("HitLabels");
  std::cout << "SBN Track Hit Dumper!\n";
  for (const std::string &label: fHitLabels) {
    std::cout << "Processing hits with label: " << label << std::endl;
  }
  
}

void sbn::TrackHitDumper::analyze(art::Event const& e)
{

  std::cout << "New Event!\n";

  // iterate over the hit tags
  for (const std::string &label: fHitLabels) {
    std::cout << "Processing hits with label: " << label << std::endl;

    auto const& hits = e.getProduct<std::vector<recob::Hit>>(label);
    for (recob::Hit const& hit: hits) {
      std::cout << "Hit on wire: " << hit.WireID() << " channel: " << hit.Channel() << " view: " << hit.View() << std::endl;
      std::cout << "Hit amp: " << hit.PeakAmplitude() << " RMS: " << hit.RMS() << " width: " << hit.SigmaPeakAmplitude();
      std::cout << " from: " << hit.StartTick() << " to: " << hit.EndTick() <<  " center: " << hit.PeakTime() << std::endl;
      std::cout << "area: " << hit.Integral() << " sim ADC: " << hit.ROISummedADC() << std::endl;
    }
  }

}

DEFINE_ART_MODULE(sbn::TrackHitDumper)
