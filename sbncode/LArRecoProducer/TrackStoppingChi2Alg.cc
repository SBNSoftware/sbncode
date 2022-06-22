#include "sbncode/LArRecoProducer/TrackStoppingChi2Alg.h"

#include "TF1.h"
#include "TGraph.h"


sbn::TrackStoppingChi2Alg::TrackStoppingChi2Alg(fhicl::ParameterSet const& p) :
  fFitRange(p.get<float>("FitRange"))
  , fMaxdEdx(p.get<float>("MaxdEdx"))
  , fMinHits(p.get<unsigned int>("MinHits"))
{
}

sbn::StoppingChi2Fit sbn::TrackStoppingChi2Alg::RunFit(const std::vector<float> &dEdxVec, const std::vector<float> &resRangeVec) const
{
  if (dEdxVec.size() != resRangeVec.size())
    throw cet::exception("TrackStoppingChi2Alg") << "dEdx and Res Range do not have same length: " << dEdxVec.size() << " and " << resRangeVec.size() << std::endl;

  if (dEdxVec.size() < fMinHits)
    return StoppingChi2Fit();

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

  return StoppingChi2Fit(pol0Chi2, expChi2, pol0Fit);
}

sbn::StoppingChi2Fit sbn::TrackStoppingChi2Alg::RunFit(const anab::Calorimetry& calo) const
{

  std::vector<float> dEdxVec, resRangeVec;
  // Fill the dEdx vs res range vectors, ignoring the first/last points
  for (size_t i = 1; i < calo.dEdx().size() - 1; i++) {
    const float thisdEdx(calo.dEdx()[i]);
    const float thisResRange(calo.ResidualRange()[i]);
    if (thisResRange > fFitRange || thisdEdx > fMaxdEdx)
      continue;

    dEdxVec.push_back(thisdEdx);
    resRangeVec.push_back(thisResRange);
  }

  return this->RunFit(dEdxVec, resRangeVec);
}

sbn::StoppingChi2Fit sbn::TrackStoppingChi2Alg::RunFitForCosmicID(const anab::Calorimetry& calo) const
{
  if(calo.XYZ().size() == 0)
    return StoppingChi2Fit();

  geo::Point_t start(calo.XYZ().front());
  geo::Point_t end(calo.XYZ().back());

  std::vector<float> dEdxVec, resRangeVec;
  // Fill the dEdx vs res range vectors, ignoring the first/last points
  for (size_t i = 1; i < calo.dEdx().size() - 1; i++) {
    const float thisdEdx(calo.dEdx()[i]);
    const float thisResRange(calo.ResidualRange()[i]);
    if (thisResRange > fFitRange || thisdEdx > fMaxdEdx)
      continue;

    dEdxVec.push_back(thisdEdx);
    resRangeVec.push_back(thisResRange);
  }

  if(fTpcGeo.MinDistToWall(start) > fTpcGeo.MinDistToWall(end))
    {
      std::reverse(dEdxVec.begin(), dEdxVec.end());
      std::reverse(resRangeVec.begin(), resRangeVec.end());
    }

  return this->RunFit(dEdxVec, resRangeVec);
}
