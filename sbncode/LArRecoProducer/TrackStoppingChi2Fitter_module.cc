////////////////////////////////////////////////////////////////////////
// Class:       TrackStoppingChi2Fitter
// Plugin Type: producer (art v3_06_03)
// File:        TrackStoppingChi2Fitter_module.cc
//
// Generated at Mon Mar  1 10:30:45 2021 by Edward Tyley using cetskelgen
// from cetlib version v3_11_01.
//
// Producer to identify Bragg peaks from stopping tracks
// Based on the StoppingParticleCosmicIDAlg by Tom Brooks in sbndcode
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
#include "sbncode/LArRecoProducer/TrackStoppingChi2Alg.h"

#include <memory>

namespace sbn {
class TrackStoppingChi2Fitter : public art::EDProducer {
  public:
  explicit TrackStoppingChi2Fitter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackStoppingChi2Fitter(TrackStoppingChi2Fitter const&) = delete;
  TrackStoppingChi2Fitter(TrackStoppingChi2Fitter&&) = delete;
  TrackStoppingChi2Fitter& operator=(TrackStoppingChi2Fitter const&) = delete;
  TrackStoppingChi2Fitter& operator=(TrackStoppingChi2Fitter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  private:
  // Declare member data here.
  const art::InputTag fTrackLabel, fCaloLabel;
  const float fMinTrackLength;

  sbn::TrackStoppingChi2Alg fTrackStoppingChi2Alg;
};

TrackStoppingChi2Fitter::TrackStoppingChi2Fitter(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fTrackLabel(p.get<std::string>("TrackLabel"))
    , fCaloLabel(p.get<std::string>("CaloLabel"))
    , fMinTrackLength(p.get<float>("MinTrackLength"))
    , fTrackStoppingChi2Alg(p)
{
  produces<std::vector<StoppingChi2Fit>>();
  produces<art::Assns<recob::Track, StoppingChi2Fit>>();
  produces<art::Assns<anab::Calorimetry, StoppingChi2Fit>>();
}

void TrackStoppingChi2Fitter::produce(art::Event& e)
{
  // Implementation of required member function here.
  auto const trackHandle(e.getValidHandle<std::vector<recob::Track>>(fTrackLabel));

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, trackHandle);

  art::FindManyP<anab::Calorimetry> fmTrackCalo(trackHandle, e, fCaloLabel);

  auto fitVec = std::make_unique<std::vector<StoppingChi2Fit>>();
  auto trackAssns = std::make_unique<art::Assns<recob::Track, StoppingChi2Fit>>();
  auto caloAssns = std::make_unique<art::Assns<anab::Calorimetry, StoppingChi2Fit>>();

  for (auto const& track : tracks) {

    if (track->Length() < fMinTrackLength)
      continue;

    const std::vector<art::Ptr<anab::Calorimetry>> caloVec(fmTrackCalo.at(track.key()));

    if (caloVec.size() != 3)
      continue;

    // Find the plane with the most hits: prefer collection > 1st induction > 2nd induction if multiple planes have the same number
    const unsigned int maxHits(std::max({ caloVec[0]->dEdx().size(), caloVec[1]->dEdx().size(), caloVec[2]->dEdx().size() }));
    const int bestPlane((caloVec[2]->dEdx().size() == maxHits) ? 2 : (caloVec[0]->dEdx().size() == maxHits) ? 0 : (caloVec[1]->dEdx().size() == maxHits) ? 1 : -1);

    if (bestPlane == -1)
      continue;

    StoppingChi2Fit thisFit(fTrackStoppingChi2Alg.RunFit(*caloVec.at(bestPlane)));

    if (thisFit.pol0Chi2 < 0.f || thisFit.expChi2 < 0.f)
      continue;

    fitVec->push_back(thisFit);
    util::CreateAssn(*this, e, *fitVec, track, *trackAssns);
    util::CreateAssn(*this, e, *fitVec, caloVec[bestPlane], *caloAssns);
  }

  e.put(std::move(fitVec));
  e.put(std::move(trackAssns));
  e.put(std::move(caloAssns));
}
}

DEFINE_ART_MODULE(sbn::TrackStoppingChi2Fitter)
