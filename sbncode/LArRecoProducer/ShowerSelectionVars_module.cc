////////////////////////////////////////////////////////////////////////
// Class:       ShowerSelectionVars
// Plugin Type: producer (art v3_03_01)
// File:        ShowerSelectionVars_module.cc
//
// Generated at Tue Sep 15 08:54:39 2020 by Edward Tyley using cetskelgen
// from cetlib version v3_08_00.
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

#include "sbnobj/Common/Reco/ShowerSelectionVars.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <memory>

#include "TGraph.h"
#include "TF1.h"

namespace sbn{
  class ShowerSelectionVars;
}

class sbn::ShowerSelectionVars : public art::EDProducer {
public:
  explicit ShowerSelectionVars(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerSelectionVars(ShowerSelectionVars const&) = delete;
  ShowerSelectionVars(ShowerSelectionVars&&) = delete;
  ShowerSelectionVars& operator=(ShowerSelectionVars const&) = delete;
  ShowerSelectionVars& operator=(ShowerSelectionVars&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // FCL params
  const art::InputTag fPandoraLabel;
  const art::InputTag fShowerLabel;
  const unsigned int fNSegments;
  const bool fRemoveStartFin;

  sbn::ShowerDensityFit DensityFitter(const recob::Shower& shower,
      const std::vector<art::Ptr<recob::SpacePoint> >& sps) const;

  sbn::ShowerTrackFit TrackFinder(const recob::Shower& shower, const recob::Track& track,
      const std::vector<art::Ptr<recob::SpacePoint> >& sps) const ;

  std::vector<float> ResidualFinder(const recob::Shower& primaryShower, const std::vector<recob::Shower>& sliceShowers) const;

  double SpacePointPerpendicular(const recob::SpacePoint& sp,
      TVector3 const& vertex, TVector3 const& direction) const;

  double SpacePointPerpendicular(const recob::SpacePoint& sp,
      TVector3 const& vertex, TVector3 const& direction, double const proj) const;

  double SpacePointProjection(const recob::SpacePoint& sp,
      TVector3 const& vertex, TVector3 const& direction) const;

  TVector3 SpacePointPosition(recob::SpacePoint const& sp) const;

};


sbn::ShowerSelectionVars::ShowerSelectionVars(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fPandoraLabel(p.get<art::InputTag>("PandoraLabel"))
  , fShowerLabel(p.get<art::InputTag>("ShowerLabel"))
  , fNSegments(p.get<unsigned int>("NSegments"))
  , fRemoveStartFin(p.get<bool>("RemoveStartFin"))
{
  if (fNSegments==0)
    throw cet::exception("ShowerSelectionVars") << "Number of Segements must be non-zero" << std::endl;

  produces<std::vector<float> >();
  produces<std::vector<sbn::ShowerTrackFit> >();
  produces<std::vector<sbn::ShowerDensityFit> >();

  produces<art::Assns<recob::Shower, float> >();
  produces<art::Assns<recob::Shower, sbn::ShowerTrackFit> >();
  produces<art::Assns<recob::Shower, sbn::ShowerDensityFit> >();
}

void sbn::ShowerSelectionVars::produce(art::Event& e)
{
  //Get the showers
  auto const showerHandle = e.getValidHandle<std::vector<recob::Shower> >(fShowerLabel);
  auto const sliceHandle = e.getValidHandle<std::vector<recob::Slice> >(fPandoraLabel);
  auto const pfpHandle = e.getValidHandle<std::vector<recob::PFParticle> >(fPandoraLabel);
  // auto const spHandle = e.getValidHandle<std::vector<recob::SpacePoint> >(fPandoraLabel);

  std::vector<art::Ptr<recob::Shower> > showers;
  art::fill_ptr_vector(showers, showerHandle);

  std::vector<art::Ptr<recob::Slice> > slices;
  art::fill_ptr_vector(slices, sliceHandle);

  std::vector<art::Ptr<recob::PFParticle> > pfps;
  art::fill_ptr_vector(pfps, pfpHandle);

  art::FindManyP<recob::PFParticle> fmSlicePFP(sliceHandle, e, fPandoraLabel);
  if(!fmSlicePFP.isValid()){
    throw cet::exception("ShowerSelectionVars") << "Slice-PFP association is somehow not valid. Stopping";
    return;
  }
  art::FindManyP<recob::Shower> fmPFPShower(pfpHandle, e, fShowerLabel);
  if(!fmPFPShower.isValid()){
    throw cet::exception("ShowerSelectionVars") << "PFP-Shower association is somehow not valid. Stopping";
    return;
  }
  art::FindManyP<recob::SpacePoint> fmShowerSP(showerHandle, e, fShowerLabel);
  if(!fmShowerSP.isValid()){
    throw cet::exception("ShowerSelectionVars") << "Shower-SP association is somehow not valid. Stopping";
    return;
  }
  art::FindManyP<recob::SpacePoint> fmTrackSP(showerHandle, e, fShowerLabel);
  if(!fmTrackSP.isValid()){
    throw cet::exception("TrackSelectionVars") << "Track-SP association is somehow not valid. Stopping";
    return;
  }
  art::FindManyP<recob::Track> fmShowerTrack(showerHandle, e, fShowerLabel);
  if(!fmShowerTrack.isValid()){
    throw cet::exception("ShowerSelectionVars") << "Shower-Track association is somehow not valid. Stopping";
    return;
  }
  // art::FindManyP<recob::Track> fmSPHit(spHandle, e, fPandoraLabel);
  // if(!fmSPHit.isValid()){
  //   throw cet::exception("ShowerSelectionVars") << "SP-Hit association is somehow not valid. Stopping";
  //   return;
  // }

  // Create a map of PFP ID to PFP object for convinience
  std::map<size_t, art::Ptr<recob::PFParticle> > pfpMap;
  for (auto const& pfp: pfps) {
    pfpMap[pfp->Self()] = pfp;
  }

  std::map<recob::Shower, std::vector<float> > showerResidualMap;

  // Calculate the residuals: Want to calculate the residuals from each primary
  // shower to every other primary shower in the same slice
  for (auto const& slice: slices) {
    std::vector<art::Ptr<recob::PFParticle> >  slicePFPs(fmSlicePFP.at(slice.key()));

    if (slicePFPs.empty())
      continue;

    art::Ptr<recob::PFParticle> pfpNeutrino;
    for (auto const& pfp: slicePFPs) {
      if(std::abs(pfp->PdgCode()) == 12 || std::abs(pfp->PdgCode()) == 14){
        pfpNeutrino = pfp;
        break;
      } // if neutrino
    } // end pfp: slicePFPs

    if (pfpNeutrino.isNull())
      continue;

    std::vector<recob::Shower> primaryShowers;
    for (auto const& daughterID: pfpNeutrino->Daughters()) {
      const auto& pfp(pfpMap[daughterID]);
      if (pfp->PdgCode()==11) {
        primaryShowers.push_back(*fmPFPShower.at(pfp.key()).front());
      } // if showers
    } // end daughterId: pfpNeutrino->Daughters

    if (primaryShowers.size()<2)
      continue;

    for (auto const& shower: primaryShowers){
      showerResidualMap[shower] = ResidualFinder(shower, primaryShowers);
    } // end shower: primaryShowers
  } // end slice: slices

  std::unique_ptr<std::vector<float> >residualCol(std::make_unique<std::vector<float> >());
  std::unique_ptr<std::vector<sbn::ShowerDensityFit> > densityFitCol(std::make_unique<std::vector<sbn::ShowerDensityFit> >());
  std::unique_ptr<std::vector<sbn::ShowerTrackFit> > trackFitCol(std::make_unique<std::vector<sbn::ShowerTrackFit> >());

  std::unique_ptr<art::Assns<recob::Shower, float> >
    residualAssns(std::make_unique<art::Assns<recob::Shower, float> >());
  std::unique_ptr<art::Assns<recob::Shower, sbn::ShowerDensityFit> >
    densityFitAssns(std::make_unique<art::Assns<recob::Shower, sbn::ShowerDensityFit> >());
  std::unique_ptr<art::Assns<recob::Shower, sbn::ShowerTrackFit> >
    trackFitAssns(std::make_unique<art::Assns<recob::Shower, sbn::ShowerTrackFit> >());

  for (auto const& shower: showers){

    for (auto const res: showerResidualMap[*shower]) {
      residualCol->push_back(res);
      util::CreateAssn(*this, e, *residualCol, shower, *residualAssns);
    }

    sbn::ShowerDensityFit showerDensityFit(DensityFitter(*shower, fmShowerSP.at(shower.key())));
    densityFitCol->push_back(showerDensityFit);
    util::CreateAssn(*this, e, *densityFitCol, shower, *densityFitAssns);

    const std::vector<art::Ptr<recob::Track> > showerTrackVec(fmShowerTrack.at(shower.key()));
    sbn::ShowerTrackFit showerTrackFit(showerTrackVec.empty() ? sbn::ShowerTrackFit() :
        TrackFinder(*shower, *showerTrackVec.front(), fmTrackSP.at(showerTrackVec.front().key())));
    trackFitCol->push_back(showerTrackFit);
    util::CreateAssn(*this, e, *trackFitCol, shower, *trackFitAssns);
  }

  e.put(std::move(residualCol));
  e.put(std::move(densityFitCol));
  e.put(std::move(trackFitCol));

  e.put(std::move(residualAssns));
  e.put(std::move(densityFitAssns));
  e.put(std::move(trackFitAssns));
}


sbn::ShowerDensityFit sbn::ShowerSelectionVars::DensityFitter(const recob::Shower& shower,
    const std::vector<art::Ptr<recob::SpacePoint> >& sps) const {

  const TVector3 showerVtx(shower.ShowerStart());
  const TVector3 showerDir(shower.Direction());
  const double OpenAngle(shower.OpenAngle());
  unsigned int totalHits(0);

  if (!shower.has_length() || !shower.has_open_angle() || sps.empty())
    return sbn::ShowerDensityFit();

  std::map<int, std::vector<art::Ptr<recob::SpacePoint> > > segmentSPMap;
  double segmentSize = shower.Length()/fNSegments;

  //Split the the spacepoints into segments.
  for(auto const& sp: sps){

    //Get the position of the spacepoint
    const double projLen(SpacePointProjection(*sp, showerVtx, showerDir));
    const double perpLen(SpacePointPerpendicular(*sp, showerVtx, showerDir, projLen));

    if(perpLen > TMath::Abs(TMath::Tan(OpenAngle)*projLen))
      continue;

    //Get where the sp should be place.
    const int sg_len(round(projLen/segmentSize));
    segmentSPMap[sg_len].push_back(sp);
    ++totalHits;
  }

  TGraph* graph = new TGraph();

  //Calculate the density gradent.
  for(auto& segment: segmentSPMap){
    double sg_len = segment.first;

    if(segment.second.size() < 10){continue;}

    //Calculate the charge in the segement
    double segmentHits = segment.second.size();

    //Calculate the voume
    double lower_dist = sg_len*segmentSize - segmentSize/2;
    double upper_dist = sg_len*segmentSize + segmentSize/2;

    if(fRemoveStartFin){if(sg_len==0 || sg_len==fNSegments){continue;}}

    if(sg_len==0)         {lower_dist = 0;}
    if(sg_len==fNSegments){upper_dist = sg_len*segmentSize;}

    double littlevolume = lower_dist*TMath::Power((TMath::Tan(0.5*OpenAngle)*lower_dist),2)*TMath::Pi()/3;
    double bigvolume    = upper_dist*TMath::Power((TMath::Tan(0.5*OpenAngle)*upper_dist),2)*TMath::Pi()/3;
    double volume       = bigvolume - littlevolume;

    double SegmentDensity = segmentHits/volume;

    double LengthToSegment = (lower_dist+upper_dist)/2;

    graph->SetPoint(graph->GetN(),LengthToSegment,SegmentDensity/totalHits);
  }

  if(graph->GetN() < 3 )
    return sbn::ShowerDensityFit();

  TF1 *fit = new TF1("fit", "[0]/x^[1]");
  fit->SetParLimits(1,1,2);
  fit->SetParLimits(0,0,1);

  graph->Fit(fit,"Q");

  const double grad(fit->GetParameter(0));
  const double pow(fit->GetParameter(1));

  delete fit;
  delete graph;

  return sbn::ShowerDensityFit(grad, pow);
}

sbn::ShowerTrackFit sbn::ShowerSelectionVars::TrackFinder(const recob::Shower& shower, const recob::Track& track,
    const std::vector<art::Ptr<recob::SpacePoint> >& sps) const {

  const unsigned int numSPs(sps.size());

  if (sps.empty())
    return sbn::ShowerTrackFit();

  double perp(0);
  const TVector3 showerVtx(shower.ShowerStart());
  const TVector3 showerDir(shower.Direction());

  for (auto const& sp: sps) {
    perp += SpacePointPerpendicular(*sp, showerVtx, showerDir);
  }

  const double trackLen(track.Length());
  const double trackWidth(perp/numSPs);

  return sbn::ShowerTrackFit(trackLen, trackWidth, numSPs);
}

std::vector<float> sbn::ShowerSelectionVars::ResidualFinder(const recob::Shower& primaryShower, const std::vector<recob::Shower>& sliceShowers) const {

  std::vector<float> residuals;

  const TVector3 primaryVtx(primaryShower.ShowerStart());
  const TVector3 showerDir(primaryShower.Direction());

  for (auto const& shower: sliceShowers) {

    if ((primaryVtx - shower.ShowerStart()).Mag() < std::numeric_limits<float>::epsilon())
      continue;

    const TVector3 showerDisplacement(primaryVtx - shower.ShowerStart());

    const double projection(showerDir.Dot(showerDisplacement));
    const double perpendicualr((showerDisplacement - projection*showerDir).Mag());

    residuals.push_back(perpendicualr);
  }

  return residuals;
}

//TODO move to a helper function (better yet replace with geo::Point_t when updating larsoft

double sbn::ShowerSelectionVars::SpacePointPerpendicular(const recob::SpacePoint& sp,
                TVector3 const& vertex,
                TVector3 const& direction) const {

  const double proj(SpacePointProjection(sp, vertex, direction));
  return SpacePointPerpendicular(sp, vertex, direction, proj);
 }

double sbn::ShowerSelectionVars::SpacePointPerpendicular(const recob::SpacePoint& sp,
                TVector3 const& vertex,
                TVector3 const& direction,
                double const proj) const {

  const TVector3 pos = SpacePointPosition(sp) - vertex;
  return (pos - proj*direction).Mag();
}

double sbn::ShowerSelectionVars::SpacePointProjection(const recob::SpacePoint& sp,
                TVector3 const& vertex,
                TVector3 const& direction) const {

  TVector3 pos = SpacePointPosition(sp) - vertex;
  return pos.Dot(direction);
}

TVector3 sbn::ShowerSelectionVars::SpacePointPosition(recob::SpacePoint const& sp) const {

  const Double32_t* sp_xyz = sp.XYZ();
  TVector3 sp_postiion = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
  return sp_postiion;
}

DEFINE_ART_MODULE(sbn::ShowerSelectionVars)
