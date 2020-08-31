////////////////////////////////////////////////////////////////////////
// Class:       TrackHitFilter
// Plugin Type: producer (art v3_02_06)
// File:        TrackHitFilter_module.cc
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
#include "lardataobj/Utilities/sparse_vector.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcorealg/Geometry/Exceptions.h"

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
  class TrackHitFilter;
}


class sbn::TrackHitFilter : public art::EDProducer {
public:
  explicit TrackHitFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackHitFilter(TrackHitFilter const&) = delete;
  TrackHitFilter(TrackHitFilter&&) = delete;
  TrackHitFilter& operator=(TrackHitFilter const&) = delete;
  TrackHitFilter& operator=(TrackHitFilter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fTrackLabel;
  art::InputTag fHitLabel;
  bool fPassBadHits;

};

sbn::TrackHitFilter::TrackHitFilter(fhicl::ParameterSet const& p)
  : EDProducer{p},
//    fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraKalmanTrackICARUSCryo0")),
   
fTrackLabel(p.get<art::InputTag>("TrackLabel", "pandoraTrackGausCryo0")),
//    fHitLabel(p.get<art::InputTag>("HitLabel", "icarusHitTPC0")),
    fHitLabel(p.get<art::InputTag>("HitLabel", "gaushitTPC0")),
    fPassBadHits(p.get<bool>("PassBadHits", false))
{
  produces<std::vector<recob::Hit>>();
  produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
  // produces<art::Assns<recob::Wire, recob::Hit>>();
}

void sbn::TrackHitFilter::produce(art::Event& e)
{
  // output data products
  std::unique_ptr<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>> assn(new art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>);
  std::unique_ptr<std::vector<recob::Hit>> outHits(new std::vector<recob::Hit>);
  // std::unique_ptr<art::Assns<recob::Wire, recob::Hit>> wireAssn(new art::Assns<recob::Wire, recob::Hit>);

  art::PtrMaker<recob::Hit> hitPtrMaker{e};

  // input data
  art::Handle<std::vector<recob::Track>> track_handle;
  e.getByLabel(fTrackLabel, track_handle);

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, track_handle);

  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmHits(tracks, e, fTrackLabel);

  for (unsigned i = 0; i < tracks.size(); i++) {
    const recob::Track &track = *tracks[i];

    // get the input hits
    const std::vector<art::Ptr<recob::Hit>> &trkHits = fmHits.at(i);
    const std::vector<const recob::TrackHitMeta*> &trkHitMetas = fmHits.data(i);

    // art::FindManyP<recob::Wire> hitWires(trkHits, e, fHitLabel);

    for (unsigned i_hit = 0; i_hit < trkHits.size(); i_hit++) {
      const recob::Hit &hit = *trkHits[i_hit];
      const recob::TrackHitMeta &meta = *trkHitMetas[i_hit];

      // figure out if bad hit -- copy of what Calorimetry module does
      bool badhit = (meta.Index() == std::numeric_limits<int>::max()) ||
                    (!track.HasValidPoint(meta.Index()));
      if (!badhit || fPassBadHits) {
        // save to output data
        outHits->push_back(hit);
        art::Ptr<recob::Hit> thisHitPtr = hitPtrMaker(outHits->size()-1);
        // add to the association
        assn->addSingle(tracks[i], thisHitPtr, meta);
        // wireAssn->addSingle(hitWires.at(i_hit).at(0), thisHitPtr);
      } 
    }
  }

  e.put(std::move(outHits));
  e.put(std::move(assn));
  // e.put(std::move(wireAssn));

}

DEFINE_ART_MODULE(sbn::TrackHitFilter)
