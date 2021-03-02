////////////////////////////////////////////////////////////////////////
// Class:       TrackScatterDCA
// Plugin Type: producer (art v3_06_03)
// File:        TrackScatterDCA_module.cc
//
// Generated at Mon Mar  1 10:30:45 2021 by Edward Tyley using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Track.h"
#include "sbnobj/Common/Reco/ScatterDCA.h"

#include <memory>

namespace sbn {
class TrackScatterDCA;
}

class sbn::TrackScatterDCA : public art::EDProducer {
  public:
  explicit TrackScatterDCA(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackScatterDCA(TrackScatterDCA const&) = delete;
  TrackScatterDCA(TrackScatterDCA&&) = delete;
  TrackScatterDCA& operator=(TrackScatterDCA const&) = delete;
  TrackScatterDCA& operator=(TrackScatterDCA&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  private:
  // Declare member data here.
  const art::InputTag fTrackLabel;
  const float fMinTrackLength;

  sbn::ScatterDCA CalculateDCA(const recob::Track& track) const;
};

sbn::TrackScatterDCA::TrackScatterDCA(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fTrackLabel(p.get<std::string>("TrackLabel"))
    , fMinTrackLength(p.get<float>("MinTrackLength"))
{
  produces<std::vector<sbn::ScatterDCA>>();
  produces<art::Assns<recob::Track, sbn::ScatterDCA>>();
}

void sbn::TrackScatterDCA::produce(art::Event& e)
{
  // Implementation of required member function here.
  auto const trackHandle(e.getValidHandle<std::vector<recob::Track>>(fTrackLabel));

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, trackHandle);

  auto dcaVec = std::make_unique<std::vector<sbn::ScatterDCA>>();
  auto trackAssns = std::make_unique<art::Assns<recob::Track, sbn::ScatterDCA>>();

  for (auto const& track : tracks) {

    if (track->Length() < fMinTrackLength)
      continue;

    sbn::ScatterDCA dca(this->CalculateDCA(*track));

    dcaVec->push_back(dca);
    util::CreateAssn(*this, e, *dcaVec, track, *trackAssns);
  }

  e.put(std::move(dcaVec));
  e.put(std::move(trackAssns));
}

sbn::ScatterDCA sbn::TrackScatterDCA::CalculateDCA(const recob::Track& track) const
{
  const TVector3 start(track.Start<TVector3>());
  const TVector3 dir((start - track.End<TVector3>()).Unit());

  float sumDCA(0), maxDCA(0);
  unsigned int counter(0);
  for (size_t i = 0; i < track.NumberTrajectoryPoints(); i++) {
    if (!track.HasValidPoint(i))
      continue;
    counter++;

    const TVector3 pos(track.LocationAtPoint<TVector3>(i));
    const TVector3 disp(pos - start);
    const float proj(disp.Dot(dir));
    const float thisDCA((disp - proj * dir).Mag());

    sumDCA += thisDCA;
    maxDCA = std::max(maxDCA, thisDCA);
  }

  if (!counter)
    return sbn::ScatterDCA();

  const float meanDCA(sumDCA / counter);

  float sumStdDev(0);
  for (size_t i = 0; i < track.NumberTrajectoryPoints(); i++) {
    if (!track.HasValidPoint(i))
      continue;

    const TVector3 pos(track.LocationAtPoint<TVector3>(i));
    const TVector3 disp(pos - start);
    const float proj(disp.Dot(dir));
    const float thisDCA((disp - proj * dir).Mag());
    const float thisDCADev(thisDCA - meanDCA);

    sumStdDev += thisDCADev * thisDCADev;
  }

  const float stdDevDCA(std::sqrt(sumStdDev / counter));

  return sbn::ScatterDCA(meanDCA, stdDevDCA, maxDCA);
}

DEFINE_ART_MODULE(sbn::TrackScatterDCA)
