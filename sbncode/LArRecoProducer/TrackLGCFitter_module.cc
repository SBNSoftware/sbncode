////////////////////////////////////////////////////////////////////////
// Class:       TrackLGCFitter
// Plugin Type: producer (art v3_06_03)
// File:        TrackLGCFitter_module.cc
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

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "sbncode/LArRecoProducer/LArReco/LGfitter.h"
#include "sbnobj/Common/Reco/LGCFit.h"

#include "TF1.h"
#include "TH1F.h"

#include <memory>

namespace sbn {
class TrackLGCFitter;
}

class sbn::TrackLGCFitter : public art::EDProducer {
  public:
  explicit TrackLGCFitter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackLGCFitter(TrackLGCFitter const&) = delete;
  TrackLGCFitter(TrackLGCFitter&&) = delete;
  TrackLGCFitter& operator=(TrackLGCFitter const&) = delete;
  TrackLGCFitter& operator=(TrackLGCFitter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  private:
  // Declare member data here.
  art::InputTag fTrackLabel, fCaloLabel;
  float fMinTrackLength;

  LGfitter::LGfitter lgFitter;

  sbn::LGCFit RunFit(const anab::Calorimetry& calo) const;
};

sbn::TrackLGCFitter::TrackLGCFitter(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fTrackLabel(p.get<std::string>("TrackLabel"))
    , fCaloLabel(p.get<std::string>("CaloLabel"))
    , fMinTrackLength(p.get<float>("MinTrackLength"))
    , lgFitter(p.get<float>("LGFitterNP", 100.f), p.get<float>("LGFitterSC", 8.f))
{
  produces<std::vector<sbn::LGCFit>>();
  produces<art::Assns<recob::Track, sbn::LGCFit>>();
  produces<art::Assns<anab::Calorimetry, sbn::LGCFit>>();
}

void sbn::TrackLGCFitter::produce(art::Event& e)
{
  // Implementation of required member function here.
  auto const trackHandle(e.getValidHandle<std::vector<recob::Track>>(fTrackLabel));

  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, trackHandle);

  art::FindManyP<anab::Calorimetry> fmTrackCalo(trackHandle, e, fCaloLabel);

  auto fitVec = std::make_unique<std::vector<sbn::LGCFit>>();
  auto trackAssns = std::make_unique<art::Assns<recob::Track, sbn::LGCFit>>();
  auto caloAssns = std::make_unique<art::Assns<anab::Calorimetry, sbn::LGCFit>>();

  for (auto const& track : tracks) {

    if (track->Length() < fMinTrackLength)
      continue;

    const std::vector<art::Ptr<anab::Calorimetry>> caloVec(fmTrackCalo.at(track.key()));

    if (caloVec.size() != 3)
      continue;

    const unsigned int maxHits(std::max(caloVec[0]->dEdx().size(), std::max(caloVec[1]->dEdx().size(), caloVec[2]->dEdx().size())));
    const int bestPlane((caloVec[2]->dEdx().size() == maxHits) ? 2 : (caloVec[0]->dEdx().size() == maxHits) ? 0 : (caloVec[1]->dEdx().size() == maxHits) ? 1 : -1);

    if (bestPlane == -1)
      continue;

    sbn::LGCFit thisFit(this->RunFit(*caloVec.at(bestPlane)));

    if (thisFit.mChi2 < 0.f)
      continue;

    fitVec->push_back(thisFit);
    util::CreateAssn(*this, e, *fitVec, track, *trackAssns);
    util::CreateAssn(*this, e, *fitVec, caloVec[bestPlane], *caloAssns);
  }

  e.put(std::move(fitVec));
  e.put(std::move(trackAssns));
  e.put(std::move(caloAssns));
}

sbn::LGCFit sbn::TrackLGCFitter::RunFit(const anab::Calorimetry& calo) const
{

  TH1F* hist = new TH1F("dEdxHist", "dEdxHist", 50, 0, 20);
  // Fill the dEdx into the histogram, ignoring the first and last points
  for (auto dEdxIter = calo.dEdx().cbegin() + 1; dEdxIter != calo.dEdx().cend() - 1; dEdxIter++)
    hist->Fill(*dEdxIter);

  // Guesstimate the best params
  double* bestFitParams = new double[4] { hist->GetRMS() / 2.f, hist->GetBinCenter(hist->GetMaximumBin()), hist->GetMaximum(), hist->GetRMS() / 2.f };
  TF1* lgcFun(lgFitter.Fit(hist, 1, 16, bestFitParams, nullptr, nullptr));

  const auto lgcFit(lgcFun ? sbn::LGCFit(bestFitParams[0], bestFitParams[1], bestFitParams[2], bestFitParams[3], lgcFun->GetChisquare(), lgcFun->GetNDF())
                           : sbn::LGCFit());

  delete[] bestFitParams;
  delete hist;
  delete lgcFun;

  return lgcFit;
}

DEFINE_ART_MODULE(sbn::TrackLGCFitter)
