////////////////////////////////////////////////////////////////////////
// Class:       RangePAllPID
// Plugin Type: producer (art v3_02_06)
// File:        RangePAllPID_module.cc
//
// Generated at Wed Feb 19 17:38:21 2020 by Gray Putnam using cetskelgen
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

#include "TDatabasePDG.h"

#include "LArReco/TrackMomentumCalculator.h"
#include "sbnobj/Common/Reco/RangeP.h"

#include <memory>

namespace sbn {
  class RangePAllPID;
}


class sbn::RangePAllPID : public art::EDProducer {
public:
  explicit RangePAllPID(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RangePAllPID(RangePAllPID const&) = delete;
  RangePAllPID(RangePAllPID&&) = delete;
  RangePAllPID& operator=(RangePAllPID const&) = delete;
  RangePAllPID& operator=(RangePAllPID&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  trkf::TrackMomentumCalculator fRangeCalculator;
  art::InputTag fTrackLabel;
  double fMuonMass;
  double fPionMass;

};

const static std::vector<int> PIDs {13, 211, 2212};
const static std::vector<std::string> names {"muon", "pion", "proton"};

sbn::RangePAllPID::RangePAllPID(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fRangeCalculator(p.get<float>("MinTrackLength", 10.)),
    fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraTrack"))
{
  for (unsigned i = 0; i < names.size(); i++) {
    produces<std::vector<sbn::RangeP>>(names[i]);
    produces<art::Assns<recob::Track, sbn::RangeP>>(names[i]); 
  }

  const TDatabasePDG *PDGTable = TDatabasePDG::Instance();
  fMuonMass = PDGTable->GetParticle(13)->Mass();
  fPionMass = PDGTable->GetParticle(211)->Mass();
}

void sbn::RangePAllPID::produce(art::Event& e)
{
  art::Handle<std::vector<recob::Track>> track_handle;
  e.getByLabel(fTrackLabel, track_handle);

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, track_handle);

  for (unsigned i = 0; i < PIDs.size(); i++) {
    std::unique_ptr<std::vector<sbn::RangeP>> rangecol(new std::vector<sbn::RangeP>);
    std::unique_ptr<art::Assns<recob::Track, sbn::RangeP>> assn(new art::Assns<recob::Track, sbn::RangeP>);

    for (const art::Ptr<recob::Track> track: tracks) {
      sbn::RangeP rangep;
      // Rescale the input and output as described by https://inspirehep.net/literature/1766384 (eq. 6.2)
      if (PIDs[i] == 211) {
        rangep.range_p = fRangeCalculator.GetTrackMomentum(track->Length() * (fMuonMass / fPionMass), 13) * (fPionMass / fMuonMass); 
      }
      else {
        rangep.range_p = fRangeCalculator.GetTrackMomentum(track->Length(), PIDs[i]);
      }
      rangep.trackID = track->ID();
      rangecol->push_back(rangep);
      util::CreateAssn(*this, e, *rangecol, track, *assn, names[i]);
    }

    e.put(std::move(rangecol), names[i]);
    e.put(std::move(assn), names[i]);
  }

}

DEFINE_ART_MODULE(sbn::RangePAllPID)
