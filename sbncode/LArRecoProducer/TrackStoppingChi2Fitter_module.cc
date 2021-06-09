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

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"

#include "TF1.h"
#include "TGraph.h"

#include <memory>

namespace sbn {
class TrackStoppingChi2Fitter;
}

class sbn::TrackStoppingChi2Fitter : public art::EDProducer {
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
  const float fMinTrackLength, fMaxTrackLength, fMaxdEdx;
  const unsigned int fMinHits;

  sbn::StoppingChi2Fit RunFit(const anab::Calorimetry& calo) const;
};

sbn::TrackStoppingChi2Fitter::TrackStoppingChi2Fitter(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fTrackLabel(p.get<std::string>("TrackLabel"))
    , fCaloLabel(p.get<std::string>("CaloLabel"))
    , fMinTrackLength(p.get<float>("MinTrackLength"))
    , fMaxTrackLength(p.get<float>("MaxTrackLength"))
    , fMaxdEdx(p.get<float>("MaxdEdx"))
    , fMinHits(p.get<unsigned int>("MinHits"))
{
  produces<std::vector<sbn::StoppingChi2Fit>>();
  produces<art::Assns<recob::Track, sbn::StoppingChi2Fit>>();
  produces<art::Assns<anab::Calorimetry, sbn::StoppingChi2Fit>>();
}

void sbn::TrackStoppingChi2Fitter::produce(art::Event& e)
{
  // Implementation of required member function here.
  auto const trackHandle(e.getValidHandle<std::vector<recob::Track>>(fTrackLabel));

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, trackHandle);

  art::FindManyP<anab::Calorimetry> fmTrackCalo(trackHandle, e, fCaloLabel);

  auto fitVec = std::make_unique<std::vector<sbn::StoppingChi2Fit>>();
  auto trackAssns = std::make_unique<art::Assns<recob::Track, sbn::StoppingChi2Fit>>();
  auto caloAssns = std::make_unique<art::Assns<anab::Calorimetry, sbn::StoppingChi2Fit>>();

  for (auto const& track : tracks) {

    if (track->Length() < fMinTrackLength)
      continue;

    const std::vector<art::Ptr<anab::Calorimetry>> caloVec(fmTrackCalo.at(track.key()));

    if (caloVec.size() != 3)
      continue;

    // Find the plane with the most hits: prefer collection > 1st induction > 2nd induction if multiple planes have the same number
    const unsigned int maxHits(std::max(caloVec[0]->dEdx().size(), std::max(caloVec[1]->dEdx().size(), caloVec[2]->dEdx().size())));
    const int bestPlane((caloVec[2]->dEdx().size() == maxHits) ? 2 : (caloVec[0]->dEdx().size() == maxHits) ? 0 : (caloVec[1]->dEdx().size() == maxHits) ? 1 : -1);

    if (bestPlane == -1)
      continue;

    sbn::StoppingChi2Fit thisFit(this->RunFit(*caloVec.at(bestPlane)));

    if (thisFit.mPol0Chi2 < 0.f || thisFit.mExpChi2 < 0.f)
      continue;

    fitVec->push_back(thisFit);
    util::CreateAssn(*this, e, *fitVec, track, *trackAssns);
    util::CreateAssn(*this, e, *fitVec, caloVec[bestPlane], *caloAssns);
  }

  e.put(std::move(fitVec));
  e.put(std::move(trackAssns));
  e.put(std::move(caloAssns));
}

sbn::StoppingChi2Fit sbn::TrackStoppingChi2Fitter::RunFit(const anab::Calorimetry& calo) const
{

  std::vector<float> dEdxVec, resRangeVec;
  // Fill the dEdx vs res range vectors, ignoring the first/last points
  for (size_t i = 1; i < calo.dEdx().size() - 1; i++) {
    const float thisdEdx(calo.dEdx()[i]);
    const float thisResRange(calo.ResidualRange()[i]);
    if (thisResRange > fMaxTrackLength || thisdEdx > fMaxdEdx)
      continue;

    dEdxVec.push_back(thisdEdx);
    resRangeVec.push_back(thisResRange);
  }

  if (dEdxVec.size() != resRangeVec.size())
    throw cet::exception("TrachStoppingChi2Fitter") << "dEdx and Res Range do not have same length: " << dEdxVec.size() << " and " << resRangeVec.size() << std::endl;

  if (dEdxVec.size() < fMinHits)
    return sbn::StoppingChi2Fit();

  const auto graph(std::make_unique<TGraph>(dEdxVec.size(), &resRangeVec[0], &dEdxVec[0]));

  // Try and fit a flat polynomial
  graph->Fit("pol0", "Q");
  const TF1* polFit = graph->GetFunction("pol0");
  const float pol0Chi2(polFit ? polFit->GetChisquare() : -5.f);
  const float pol0Fit(polFit ? polFit->GetParameter(0) : -5.f);

  // Try to fit an exponential
  graph->Fit("expo", "Q");
  const TF1* expFit = graph->GetFunction("expo");
  const float expChi2(expFit ? expFit->GetChisquare() : -5.f);

  return sbn::StoppingChi2Fit(pol0Chi2, expChi2, pol0Fit);
}

DEFINE_ART_MODULE(sbn::TrackStoppingChi2Fitter)
