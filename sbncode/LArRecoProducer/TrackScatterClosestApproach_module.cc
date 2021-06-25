////////////////////////////////////////////////////////////////////////
// Class:       TrackScatterClosestApproach
// Plugin Type: producer (art v3_06_03)
// File:        TrackScatterClosestApproach_module.cc
//
// Generated at Mon Mar  1 10:30:45 2021 by Edward Tyley using cetskelgen
// from cetlib version v3_11_01.
//
// Producer to look at the average distance of closest approach from
// the centroid of a track to quantify the scattering
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
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"

#include <memory>

namespace sbn {
class TrackScatterClosestApproach : public art::EDProducer {
  public:
  explicit TrackScatterClosestApproach(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackScatterClosestApproach(TrackScatterClosestApproach const&) = delete;
  TrackScatterClosestApproach(TrackScatterClosestApproach&&) = delete;
  TrackScatterClosestApproach& operator=(TrackScatterClosestApproach const&) = delete;
  TrackScatterClosestApproach& operator=(TrackScatterClosestApproach&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  private:
  // Declare member data here.
  const art::InputTag fTrackLabel;
  const float fMinTrackLength;

  ScatterClosestApproach CalculateClosestApproach(const recob::Track& track) const;
};

TrackScatterClosestApproach::TrackScatterClosestApproach(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fTrackLabel(p.get<std::string>("TrackLabel"))
    , fMinTrackLength(p.get<float>("MinTrackLength"))
{
  produces<std::vector<ScatterClosestApproach>>();
  produces<art::Assns<recob::Track, ScatterClosestApproach>>();
}

void TrackScatterClosestApproach::produce(art::Event& e)
{
  // Implementation of required member function here.
  auto const trackHandle(e.getValidHandle<std::vector<recob::Track>>(fTrackLabel));

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, trackHandle);

  auto closestapproachVec = std::make_unique<std::vector<ScatterClosestApproach>>();
  auto trackAssns = std::make_unique<art::Assns<recob::Track, ScatterClosestApproach>>();

  for (auto const& track : tracks) {

    if (track->Length() < fMinTrackLength)
      continue;

    ScatterClosestApproach closestapproach(this->CalculateClosestApproach(*track));

    closestapproachVec->push_back(closestapproach);
    util::CreateAssn(*this, e, *closestapproachVec, track, *trackAssns);
  }

  e.put(std::move(closestapproachVec));
  e.put(std::move(trackAssns));
}

ScatterClosestApproach TrackScatterClosestApproach::CalculateClosestApproach(const recob::Track& track) const
{
  // Interpolate the start and end of the track to find the centroid axis
  const TVector3 start(track.Start<TVector3>());
  const TVector3 dir((start - track.End<TVector3>()).Unit());

  // Calculate the perpendicular distance from the centroid to each traj point
  float sumClosestApproach(0), maxClosestApproach(0);
  unsigned int counter(0);
  for (size_t i = 0; i < track.NumberTrajectoryPoints(); i++) {
    if (!track.HasValidPoint(i))
      continue;
    counter++;

    const TVector3 pos(track.LocationAtPoint<TVector3>(i));
    const TVector3 disp(pos - start);
    const float proj(disp.Dot(dir));
    const float thisClosestApproach((disp - proj * dir).Mag());

    sumClosestApproach += thisClosestApproach;
    maxClosestApproach = std::max(maxClosestApproach, thisClosestApproach);
  }

  if (!counter)
    return ScatterClosestApproach();

  const float meanClosestApproach(sumClosestApproach / counter);

  // Calculate the spread in ClosestApproach around the mean value
  float sumStdDev(0);
  for (size_t i = 0; i < track.NumberTrajectoryPoints(); i++) {
    if (!track.HasValidPoint(i))
      continue;

    const TVector3 pos(track.LocationAtPoint<TVector3>(i));
    const TVector3 disp(pos - start);
    const float proj(disp.Dot(dir));
    const float thisClosestApproach((disp - proj * dir).Mag());
    const float thisClosestApproachDev(thisClosestApproach - meanClosestApproach);

    sumStdDev += thisClosestApproachDev * thisClosestApproachDev;
  }

  const float stdDevClosestApproach(std::sqrt(sumStdDev / counter));

  return ScatterClosestApproach(meanClosestApproach, stdDevClosestApproach, maxClosestApproach);
}
}

DEFINE_ART_MODULE(sbn::TrackScatterClosestApproach)
