////////////////////////////////////////////////////////////////////////
// Class:       Razzle
// Plugin Type: producer (art v3_05_01)
// File:        Razzle_module.cc
//
// Generated at Tue Jan 26 08:37:49 2021 by Edward Tyley using cetskelgen
// from cetlib version v3_10_00.
//
// Module that attempts to classify each shower as either:
// - Electron (11)
// - Photon (22)
// - Other (0)
//
// If a shower is not classified (becuase it is too small) it will be
// assigned a pdg of -5
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/DataViewImpl.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/ShowerSelectionVars.h"

// Root Includes
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"

#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include <iostream>
#include <vector>

namespace sbn {
class Razzle;
}

class sbn::Razzle : public art::EDProducer {
  public:
  explicit Razzle(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Razzle(Razzle const&) = delete;
  Razzle(Razzle&&) = delete;
  Razzle& operator=(Razzle const&) = delete;
  Razzle& operator=(Razzle&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginJob() override;

  private:
  // Declare member data here.
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  art::InputTag fLArGeantLabel, fPFPLabel, fShowerLabel, fShowerSelVarsLabel;
  const float fMinShowerEnergy;
  const bool fMakeTree, fRunMVA;
  const std::string fMethodName, fWeightFile;

  // The metrics actually used in the MVA
  float bestdEdx;      // The dE/dx at the start of the shower (in the best plane)
  float convGap;       // The gap between the shower start and parent vertex
  float openAngle;     // Opening Angle of the shower, defined as atan(width/length)
  float modHitDensity; // The hit density corrected for the pitch (hits/wire)
  // Theory says the length grows like the log of the energy but in practice the sqrt gives a flatter distribution
  // float logEnergyDensity; // The log of the energy divided by the length
  float sqrtEnergyDensity; // Sqrt of the energy divided by the length

  // MVA scores
  float electronScore, photonScore, otherScore, bestScore;
  int bestPDG;

  // Other metrics for filling in tree for analysis
  TMVA::Reader* reader;
  TTree* showerTree;

  int trackHits;
  int truePdg, numHits, bestPlane, recoContained, recoPrimary, numDaughters;

  float length, bestEnergy, bestPlaneHits, bestPitch, logEnergyDensity;
  float trackScore;
  float densityFitGrad, densityFitPow;
  float trackLength, trackWidth;
  float startX, startY, startZ, endX, endY, endZ, trueStartX, trueStartY, trueStartZ, trueEndX, trueEndY, trueEndZ, startDist, endDist, trueP, energyComp, energyPurity;

  std::string trueType, trueEndProcess;

  void ClearTree();
  void FillTrueParticleMetrics(const detinfo::DetectorClocksData& clockData, const recob::Shower& shower, const std::vector<art::Ptr<recob::Hit>>& hits,
      std::vector<art::Ptr<sim::SimChannel>>& simChannels);
  void FillShowerMetrics(const recob::Shower& shower, const std::vector<art::Ptr<recob::Hit>>& hitVec);
  void FillPFPMetrics(const art::Ptr<recob::PFParticle>& pfp, const std::map<size_t, art::Ptr<recob::PFParticle>>& pfpMap,
      const recob::Shower& shower, const art::FindManyP<larpandoraobj::PFParticleMetadata>& fmMeta, const art::FindManyP<recob::Vertex>& fmVertex);
  void FillDensityFitMetrics(const sbn::ShowerDensityFit& densityFit);
  void FillTrackFitMetrics(const sbn::ShowerTrackFit& trackFit);
  sbn::MVAPID RunMVA();

  std::map<size_t, art::Ptr<recob::PFParticle>> GetPFPMap(std::vector<art::Ptr<recob::PFParticle>>& pfps) const;
  float GetPFPTrackScore(const art::Ptr<recob::PFParticle>& pfp, const art::FindManyP<larpandoraobj::PFParticleMetadata>& fmMeta) const;
  bool InFV(const TVector3& pos) const;
  std::string PdgString(const int pdg) const;
};

sbn::Razzle::Razzle(fhicl::ParameterSet const& p)
    : EDProducer { p }
    , fLArGeantLabel(p.get<std::string>("LArGeantLabel"))
    , fPFPLabel(p.get<std::string>("PFPLabel"))
    , fShowerLabel(p.get<std::string>("ShowerLabel"))
    , fShowerSelVarsLabel(p.get<std::string>("ShowerSelVarsLabel"))
    , fMinShowerEnergy(p.get<float>("MinShowerEnergy"))
    , fMakeTree(p.get<bool>("MakeTree"))
    , fRunMVA(p.get<bool>("RunMVA"))
    , fMethodName(p.get<std::string>("MethodName", ""))
    , fWeightFile(p.get<std::string>("WeightFile", ""))
{
  if (!fMakeTree && !fRunMVA)
    throw cet::exception("Razzle") << "Configured to do nothing";

  if (fRunMVA) {
    if (fMethodName == "" || fWeightFile == "")
      throw cet::exception("Razzle") << "Trying to run MVA with inputs not set: MethodName: " << fMethodName << " and WeightFile: " << fWeightFile;

    cet::search_path searchPath("FW_SEARCH_PATH");
    std::string fWeightFileFullPath;
    if (!searchPath.find_file(fWeightFile, fWeightFileFullPath))
      throw cet::exception("Razzle") << "Unable to find weight file: " << fWeightFile << " in FW_SEARCH_PATH: " << searchPath.to_string();

    reader = new TMVA::Reader("V");

    // std::cout << "Adding Variable: recoLen" << std::endl;
    // reader->AddVariable("recoLen", &recoLen);
    reader->AddVariable("bestdEdx", &bestdEdx);
    reader->AddVariable("convGap", &convGap);
    reader->AddVariable("openAngle", &openAngle);
    reader->AddVariable("modHitDensity", &modHitDensity);
    // reader->AddVariable("logEnergyDensity", &logEnergyDensity);
    reader->AddVariable("sqrtEnergyDensity", &sqrtEnergyDensity);

    reader->BookMVA(fMethodName, fWeightFileFullPath);
  }
  // Call appropriate produces<>() functions here.
  produces<std::vector<sbn::MVAPID>>();
  produces<art::Assns<recob::Shower, sbn::MVAPID>>();
}

void sbn::Razzle::beginJob()
{
  if (fMakeTree) {
    showerTree = tfs->make<TTree>("showerTree", "Tree filled per Shower with  PID variables");

    if (fRunMVA) {
      showerTree->Branch("electronScore", &electronScore);
      showerTree->Branch("photonScore", &photonScore);
      showerTree->Branch("otherScore", &otherScore);
      showerTree->Branch("bestScore", &bestScore);
      showerTree->Branch("bestPDG", &bestPDG);
    }

    showerTree->Branch("truePdg", &truePdg);
    showerTree->Branch("trueP", &trueP);
    showerTree->Branch("trueType", &trueType);
    showerTree->Branch("trueEndProcess", &trueEndProcess);
    showerTree->Branch("energyComp", &energyComp);
    showerTree->Branch("energyPurity", &energyPurity);

    showerTree->Branch("startX", &startX);
    showerTree->Branch("startY", &startY);
    showerTree->Branch("startZ", &startZ);

    showerTree->Branch("endX", &endX);
    showerTree->Branch("endY", &endY);
    showerTree->Branch("endZ", &endZ);

    showerTree->Branch("trueStartX", &trueStartX);
    showerTree->Branch("trueStartY", &trueStartY);
    showerTree->Branch("trueStartZ", &trueStartZ);

    showerTree->Branch("trueEndX", &trueEndX);
    showerTree->Branch("trueEndY", &trueEndY);
    showerTree->Branch("trueEndZ", &trueEndZ);

    showerTree->Branch("startDist", &startDist);
    showerTree->Branch("endDist", &endDist);

    showerTree->Branch("numHits", &numHits);
    showerTree->Branch("bestPlane", &bestPlane);

    showerTree->Branch("recoPrimary", &recoPrimary);
    showerTree->Branch("numDaughters", &numDaughters);
    showerTree->Branch("trackScore", &trackScore);
    showerTree->Branch("convGap", &convGap);
    showerTree->Branch("recoContained", &recoContained);

    showerTree->Branch("densityFitGrad", &densityFitGrad);
    showerTree->Branch("densityFitPow", &densityFitPow);

    showerTree->Branch("trackLength", &trackLength);
    showerTree->Branch("trackWidth", &trackWidth);
    showerTree->Branch("trackHits", &trackHits);

    showerTree->Branch("length", &length);
    showerTree->Branch("openAngle", &openAngle);
    showerTree->Branch("bestdEdx", &bestdEdx);
    showerTree->Branch("bestEnergy", &bestEnergy);
    showerTree->Branch("bestPlaneHits", &bestPlaneHits);
    showerTree->Branch("bestPitch", &bestPitch);
    showerTree->Branch("modHitDensity", &modHitDensity);
    showerTree->Branch("sqrtEnergyDensity", &sqrtEnergyDensity);
    showerTree->Branch("logEnergyDensity", &logEnergyDensity);
  }
}

void sbn::Razzle::produce(art::Event& e)
{
  auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));

  auto const pfpHandle(e.getValidHandle<std::vector<recob::PFParticle>>(fPFPLabel));
  auto const simChannelHandle(e.getValidHandle<std::vector<sim::SimChannel>>(fLArGeantLabel));
  auto const showerHandle(e.getValidHandle<std::vector<recob::Shower>>(fShowerLabel));

  std::vector<art::Ptr<recob::PFParticle>> pfps;
  art::fill_ptr_vector(pfps, pfpHandle);

  std::vector<art::Ptr<sim::SimChannel>> simChannels;
  art::fill_ptr_vector(simChannels, simChannelHandle);

  art::FindManyP<recob::Vertex> fmPFPVertex(pfpHandle, e, fPFPLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> fmPFPMeta(pfpHandle, e, fPFPLabel);

  art::FindManyP<recob::Shower> fmPFPShower(pfpHandle, e, fShowerLabel);
  art::FindManyP<recob::Hit> fmShowerHit(showerHandle, e, fShowerLabel);

  art::FindManyP<sbn::ShowerDensityFit> fmShowerDensityFit(showerHandle, e, fShowerSelVarsLabel);
  art::FindManyP<sbn::ShowerTrackFit> fmShowerTrackFit(showerHandle, e, fShowerSelVarsLabel);

  auto mvaPIDVec = std::make_unique<std::vector<sbn::MVAPID>>();
  auto showerAssns = std::make_unique<art::Assns<recob::Shower, sbn::MVAPID>>();

  const std::map<size_t, art::Ptr<recob::PFParticle>> pfpMap(this->GetPFPMap(pfps));

  for (auto const& pfp : pfps) {
    this->ClearTree();

    // Get the shower for the PFP
    std::vector<art::Ptr<recob::Shower>> pfpShowerVec(fmPFPShower.at(pfp.key()));

    // Skip showers and primaries
    if (pfpShowerVec.empty())
      continue;

    // There should be 1-1 PFP to showers
    if (pfpShowerVec.size() > 1)
      throw cet::exception("Razzle") << "Too many showers: " << pfpShowerVec.size();

    art::Ptr<recob::Shower>& pfpShower(pfpShowerVec.front());

    if (pfpShower->best_plane() < 0 || pfpShower->Energy().at(pfpShower->best_plane()) < fMinShowerEnergy)
      continue;

    std::vector<art::Ptr<recob::Hit>> showerHitVec(fmShowerHit.at(pfpShower.key()));

    this->FillShowerMetrics(*pfpShower, showerHitVec);

    this->FillPFPMetrics(pfp, pfpMap, *pfpShower, fmPFPMeta, fmPFPVertex);

    auto const densityFitVec(fmShowerDensityFit.at(pfpShower.key()));
    if (densityFitVec.size() == 1)
      this->FillDensityFitMetrics(*densityFitVec.front());

    auto const trackFitVec(fmShowerTrackFit.at(pfpShower.key()));
    if (trackFitVec.size() == 1)
      this->FillTrackFitMetrics(*trackFitVec.front());

    if (fRunMVA) {
      sbn::MVAPID mvaPID(this->RunMVA());
      mvaPIDVec->push_back(mvaPID);
      util::CreateAssn(*this, e, *mvaPIDVec, pfpShower, *showerAssns);
    }

    // Only fill the truth metrics if we are saving a TTree
    if (fMakeTree) {
      this->FillTrueParticleMetrics(clockData, *pfpShower, showerHitVec, simChannels);
      showerTree->Fill();
    }
  }
  e.put(std::move(mvaPIDVec));
  e.put(std::move(showerAssns));
}

void sbn::Razzle::ClearTree()
{
  truePdg = -5;
  numHits = -5;
  bestPlane = -5;

  recoContained = -5;
  recoPrimary = -5;
  numDaughters = -5;

  trackScore = -5.f;
  convGap = -5.f;

  electronScore = -5.f;
  photonScore = -5.f;
  otherScore = -5.f;
  bestScore = -5.f;
  bestPDG = -5;

  startX = -999.f;
  startY = -999.f;
  startZ = -999.f;
  endX = -999.f;
  endY = -999.f;
  endZ = -999.f;
  trueStartX = -999.f;
  trueStartY = -999.f;
  trueStartZ = -999.f;
  trueEndX = -999.f;
  trueEndY = -999.f;
  trueEndZ = -999.f;
  startDist = -999.f;
  endDist = -999.f;

  trueP = -5.f;
  energyPurity = -5.f;
  energyComp = -5.f;

  densityFitGrad = -5.f;
  densityFitPow = -5.f;

  trackWidth = -5.f;
  trackLength = -5.f;
  trackHits = -5;

  length = -5.f;
  openAngle = -5.f;
  bestdEdx = -5.f;
  bestEnergy = -5.f;
  bestPlaneHits = -5.f;
  bestPitch = -5.f;
  modHitDensity = -5.f;
  logEnergyDensity = -5.f;
  sqrtEnergyDensity = -5.f;

  trueType = "";
  trueEndProcess = "";
}

void sbn::Razzle::FillTrueParticleMetrics(const detinfo::DetectorClocksData& clockData, const recob::Shower& shower, const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<art::Ptr<sim::SimChannel>>& simChannels)
{
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  const int bestMatch(TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, hits, true));

  if (!TruthMatchUtils::Valid(bestMatch))
    return;

  float totalHitEnergy(0.f), totalTrueHitEnergy(0.f);
  for (auto const& hit : hits) {
    const std::vector<sim::TrackIDE> trackIDEs(bt_serv->HitToTrackIDEs(clockData, hit));
    totalHitEnergy = std::accumulate(trackIDEs.cbegin(), trackIDEs.cend(), totalHitEnergy,
        [](float sum, auto const& ide) { return sum + ide.energy; });
    totalTrueHitEnergy = std::accumulate(trackIDEs.cbegin(), trackIDEs.cend(), totalTrueHitEnergy,
        [bestMatch](float sum, auto const& ide) { return (std::abs(ide.trackID) == bestMatch) ? sum + ide.energy : sum; });
  }

  const std::vector<const sim::IDE*> trackIDEs(bt_serv->TrackIdToSimIDEs_Ps(bestMatch));
  float totalTrueEnergy(std::accumulate(trackIDEs.cbegin(), trackIDEs.cend(), 0.f,
      [](float sum, auto const& ide) { return sum + ide->energy; }));

  energyComp = totalTrueHitEnergy / totalTrueEnergy;
  energyPurity = totalTrueHitEnergy / totalHitEnergy;

  const simb::MCParticle* const trueParticle(particleInventory->TrackIdToParticle_P(bestMatch));

  if (!trueParticle)
    return;

  truePdg = trueParticle->PdgCode();
  trueType = this->PdgString(truePdg);
  trueEndProcess = trueParticle->EndProcess();

  trueP = trueParticle->P();

  const TVector3 showerStart(shower.ShowerStart());
  const TVector3 showerEnd(showerStart + shower.Length() * shower.Direction());

  // Select first traj point where the photon loses energy, last be default
  TVector3 PositionTrajStart(trueParticle->Position(0).Vect());

  // For phoyons say the shower started when the photon loses 10% of the energy
  if (trueParticle->PdgCode() == 22) {
    unsigned int trajPoint(1);
    while (trueParticle->E(trajPoint) == trueParticle->E()) {
      trajPoint++;
    }
    PositionTrajStart = trueParticle->Position(trajPoint).Vect();
  }

  const TVector3 trueStart(PositionTrajStart);
  const TVector3 trueEnd(trueParticle->EndPosition().Vect());

  trueStartX = trueStart.X();
  trueStartY = trueStart.Y();
  trueStartZ = trueStart.Z();

  trueEndX = trueEnd.X();
  trueEndY = trueEnd.Y();
  trueEndZ = trueEnd.Z();

  startDist = (trueStart - showerStart).Mag();
  endDist = (trueEnd - showerEnd).Mag();
}

void sbn::Razzle::FillShowerMetrics(const recob::Shower& shower, const std::vector<art::Ptr<recob::Hit>>& hitVec)
{
  const geo::GeometryCore* geom = lar::providerFrom<geo::Geometry>();

  length = shower.Length();
  openAngle = shower.OpenAngle();

  const TVector3 start(shower.ShowerStart());
  const TVector3 end(start + length * shower.Direction());
  recoContained = this->InFV(start) && this->InFV(end);

  numHits = hitVec.size();

  std::array<int, 3> showerPlaneHits { 0, 0, 0 };
  for (auto const& hit : hitVec) {
    showerPlaneHits[hit->WireID().Plane]++;
  }

  std::array<float, 3> showerPlanePitches { -1.f, -1.f, -1.f };
  for (geo::PlaneGeo const& plane : geom->IteratePlanes()) {

    const float angleToVert(geom->WireAngleToVertical(plane.View(), plane.ID()) - 0.5 * M_PI);
    const float cosgamma(std::abs(std::sin(angleToVert) * shower.Direction().Y() + std::cos(angleToVert) * shower.Direction().Z()));

    showerPlanePitches[plane.ID().Plane] = plane.WirePitch() / cosgamma;
  }

  // Fill only for the best plane, defined as the one with the most hits
  // Prefer collection plane > 1st induction > 2nd induction
  bestPlane = shower.best_plane();

  if (bestPlane < 0 || bestPlane > 3)
    throw cet::exception("Razzle") << "Best plane: " << bestPlane;

  bestdEdx = shower.dEdx()[bestPlane];
  bestdEdx = std::min(bestdEdx, 20.f);
  bestdEdx = std::max(bestdEdx, -5.f);

  bestEnergy = shower.Energy()[bestPlane];
  bestPlaneHits = showerPlaneHits[bestPlane];
  bestPitch = showerPlanePitches[bestPlane];

  logEnergyDensity = (length > 0 && bestEnergy > 0) ? std::log(bestEnergy) / length : -5.f;
  sqrtEnergyDensity = (length > 0 && bestEnergy > 0) ? std::sqrt(bestEnergy) / length : -5.f;

  const float wiresHit(bestPitch > std::numeric_limits<float>::epsilon() ? length / bestPitch : -5.f);
  modHitDensity = wiresHit > 1.f ? bestPlaneHits / wiresHit : -5.f;
  modHitDensity = std::min(modHitDensity, 40.f);

  if (!fMakeTree)
    return;

  startX = start.X();
  startY = start.Y();
  startZ = start.Z();

  endX = end.X();
  endY = end.Y();
  endZ = end.Z();
}

bool sbn::Razzle::InFV(const TVector3& pos) const
{
  return (std::abs(pos.X()) < 195 && std::abs(pos.Y()) < 195 && pos.Z() > 5 && pos.Z() < 495);
}

void sbn::Razzle::FillPFPMetrics(const art::Ptr<recob::PFParticle>& pfp, const std::map<size_t, art::Ptr<recob::PFParticle>>& pfpMap,
    const recob::Shower& shower, const art::FindManyP<larpandoraobj::PFParticleMetadata>& fmMeta, const art::FindManyP<recob::Vertex>& fmVertex)
{
  numDaughters = pfp->Daughters().size();
  trackScore = this->GetPFPTrackScore(pfp, fmMeta);

  auto const parentId(pfp->Parent());
  auto const& parentIter(pfpMap.find(parentId));

  if (parentIter == pfpMap.end())
    return;

  recoPrimary = parentIter->second->IsPrimary();

  auto const& parentVertexVec = fmVertex.at(parentIter->second.key());

  if (parentVertexVec.empty())
    return;

  // Need to convert vertex to TVector3 for ease of comparison to shower start
  const auto parentVertexPoint(parentVertexVec.front()->position());
  const TVector3 parentVertex { parentVertexPoint.X(), parentVertexPoint.Y(), parentVertexPoint.Z() };

  convGap = (shower.ShowerStart() - parentVertex).Mag();
  convGap = std::min(convGap, 50.f);
}

void sbn::Razzle::FillDensityFitMetrics(const sbn::ShowerDensityFit& densityFit)
{
  densityFitGrad = densityFit.mDensityGrad;
  densityFitPow = densityFit.mDensityPow;
}
void sbn::Razzle::FillTrackFitMetrics(const sbn::ShowerTrackFit& trackFit)
{
  trackLength = trackFit.mTrackLength;
  trackWidth = trackFit.mTrackWidth;
  trackHits = trackFit.mNumHits;
}

sbn::MVAPID sbn::Razzle::RunMVA()
{
  const std::vector<float> mvaScores(reader->EvaluateMulticlass(fMethodName));

  sbn::MVAPID pidResults;

  pidResults.AddScore(11, mvaScores.at(0));
  pidResults.AddScore(22, mvaScores.at(1));
  pidResults.AddScore(0, mvaScores.at(2));

  if (!fMakeTree)
    return pidResults;

  electronScore = mvaScores.at(0);
  photonScore = mvaScores.at(1);
  otherScore = mvaScores.at(2);

  bestScore = pidResults.BestScore();
  bestPDG = pidResults.BestPDG();

  return pidResults;
}

std::map<size_t, art::Ptr<recob::PFParticle>> sbn::Razzle::GetPFPMap(std::vector<art::Ptr<recob::PFParticle>>& pfps) const
{
  std::map<size_t, art::Ptr<recob::PFParticle>> pfpMap;
  for (auto const& pfp : pfps) {
    pfpMap[pfp->Self()] = pfp;
  }
  return pfpMap;
}

float sbn::Razzle::GetPFPTrackScore(const art::Ptr<recob::PFParticle>& pfp, const art::FindManyP<larpandoraobj::PFParticleMetadata>& fmMeta) const
{
  auto const pfpMetaVec(fmMeta.at(pfp.key()));
  if (pfpMetaVec.size() != 1)
    return -5.f;
  const larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap(pfpMetaVec.front()->GetPropertiesMap());
  auto const& pfpTrackScoreIter = propertiesMap.find("TrackScore");
  return pfpTrackScoreIter == propertiesMap.end() ? -5.f : pfpTrackScoreIter->second;
}

std::string sbn::Razzle::PdgString(const int pdg) const
{
  switch (std::abs(pdg)) {
  case 11:
    return "Electron";
  case 22:
    return "Photon";
  default:
    return "Other";
  }
}

DEFINE_ART_MODULE(sbn::Razzle)
