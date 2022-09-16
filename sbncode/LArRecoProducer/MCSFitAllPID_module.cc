////////////////////////////////////////////////////////////////////////
// Class:       MCSFitAllPID
// Plugin Type: producer (art v3_02_06)
// File:        MCSFitAllPID_module.cc
//
// Generated at Wed Feb 19 17:38:15 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_07_02.
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
#include "lardata/Utilities/AssociationUtil.h"

#include "LArReco/TrajectoryMCSFitter.h"

#include <memory>

namespace sbn {
  class MCSFitAllPID;
}


class sbn::MCSFitAllPID : public art::EDProducer {
public:
  explicit MCSFitAllPID(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCSFitAllPID(MCSFitAllPID const&) = delete;
  MCSFitAllPID(MCSFitAllPID&&) = delete;
  MCSFitAllPID& operator=(MCSFitAllPID const&) = delete;
  MCSFitAllPID& operator=(MCSFitAllPID&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  trkf::sbn::TrajectoryMCSFitter fMCSCalculator;
  art::InputTag fTrackLabel;
  float fMinTrackLength;
};

const static std::vector<int> PIDs {13, 211, 321, 2212};
const static std::vector<std::string> names {"muon", "pion", "kaon", "proton"};

sbn::MCSFitAllPID::MCSFitAllPID(fhicl::ParameterSet const& p)
  : EDProducer{p},
    // fMCSCalculator(p.get<fhicl::Table<trkf::TrajectoryMCSFitter::Config>>("MCS")),
    fMCSCalculator(p.get<fhicl::ParameterSet>("MCS")),
    fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraTrack")),
    fMinTrackLength(p.get<float>("MinTrackLength", 10.))
{
  for (unsigned i = 0; i < names.size(); i++) {
    produces<std::vector<recob::MCSFitResult>>(names[i]);
    produces<art::Assns<recob::Track, recob::MCSFitResult>>(names[i]); 
  }
}

void sbn::MCSFitAllPID::produce(art::Event& e)
{

  // Implementation of required member function here.
  art::Handle<std::vector<recob::Track>> track_handle;
  e.getByLabel(fTrackLabel, track_handle);

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, track_handle);

  for (unsigned i = 0; i < PIDs.size(); i++) {
    std::unique_ptr<std::vector<recob::MCSFitResult>> mcscol(new std::vector<recob::MCSFitResult>);
    std::unique_ptr<art::Assns<recob::Track, recob::MCSFitResult>> assn(new art::Assns<recob::Track, recob::MCSFitResult>);

    for (const art::Ptr<recob::Track> track: tracks) {
      if (fMinTrackLength > 0. && track->Length() < fMinTrackLength) continue;

      mcscol->push_back(fMCSCalculator.fitMcs(*track, PIDs[i]));        
      util::CreateAssn(*this, e, *mcscol, track, *assn, names[i]);
    }

    e.put(std::move(mcscol), names[i]);
    e.put(std::move(assn), names[i]);
  }

}

DEFINE_ART_MODULE(sbn::MCSFitAllPID)
