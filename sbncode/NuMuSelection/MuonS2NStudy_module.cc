////////////////////////////////////////////////////////////////////////
// Class:       MuonS2NStudy
// Plugin Type: analyzer (art v3_03_01)
// File:        MuonS2NStudy_module.cc
//
// Generated at Wed Mar 25 12:10:37 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "art_root_io/TFileService.h"

#include "lardataobj/AnalysisBase/T0.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "sbncode/LArRecoProducer/Products/RangeP.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcorealg/GeoAlgo/GeoAlgo.h"

#include "TH1D.h"
#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"

namespace numu {
  class MuonS2NStudy;
}

class numu::MuonS2NStudy : public art::EDAnalyzer {
public:
  explicit MuonS2NStudy(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MuonS2NStudy(MuonS2NStudy const&) = delete;
  MuonS2NStudy(MuonS2NStudy&&) = delete;
  MuonS2NStudy& operator=(MuonS2NStudy const&) = delete;
  MuonS2NStudy& operator=(MuonS2NStudy&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  std::vector<std::array<TH1D*,3>> fSignal;
  std::array<TH1D*, 3> fNoise;
  std::string fHitLabel;
  std::string fWireLabel;

};

numu::MuonS2NStudy::MuonS2NStudy(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},  // ,
    fHitLabel(p.get<std::string>("HitLabel")),
    fWireLabel(p.get<std::string>("WireLabel"))
  // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs; 

  std::vector<std::string> planes {"U", "V", "Y"};
  std::vector<std::string> angles {"abs(#theta) < 0.2", "#theta < -0.2", "#theta > 0.2"};
  // initilize the signal and noise on each plane
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) {
      fSignal.emplace_back();
      fSignal[j][i] = tfs->make<TH1D>(("Signal " + angles[j] + " " + planes[i]).c_str(), "Signal", 200, 0., 400.);
    }
    fNoise[i] = tfs->make<TH1D>(("Noise " + planes[i]).c_str(), "Noise", 200, 0., 10.);
  }
}

void numu::MuonS2NStudy::analyze(art::Event const& ev)
{
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();

  // Get the Signal Amplitude
  const std::vector<simb::MCParticle> &mcparticle_list = *ev.getValidHandle<std::vector<simb::MCParticle>>("largeant");

  art::Handle<std::vector<recob::Hit>> hit_handle;
  ev.getByLabel(fHitLabel, hit_handle);
  std::vector<art::Ptr<recob::Hit>> hits;
  art::fill_ptr_vector(hits, hit_handle);

  const simb::MCParticle *muon = NULL;

  // get the muon
  for (unsigned i = 0; i < mcparticle_list.size(); i++) {
    if (mcparticle_list[i].PdgCode() == 13 && mcparticle_list[i].Process() == "primary") {
      assert(muon == NULL);
      muon = &mcparticle_list[i];
    }
    break;
  }

  assert(muon != NULL);

  std::cout << "Muon Momentum Theta: " << muon->Momentum().Theta() << std::endl;
  std::cout << "Muon Momentum Phi: " << muon->Momentum().Phi() << std::endl;

  // find the angular bin
  int angle_bin = -1;
  if (abs(muon->Momentum().Theta()) < 0.2) angle_bin = 0;
  else if (muon->Momentum().Phi() < 0.) angle_bin = 1;
  else angle_bin = 2;

  // make sure it's an energetic enough muon
  // if (muon->E() < fMuonEnergyCut) return;

  // make sure the muon is mostly perpindicular to the wire-planes
  // if (muon->Momentum()->CosTheta() < fMuonCosThetaCut) return;

  // get the true matched hits
  std::vector<art::Ptr<recob::Hit>> muon_hits = backtracker->TrackIdToHits_Ps(muon->TrackId(), hits);

  // Only consider the hit with the highest amplitude on each channel
  // (in case of accidental delta rays or split hits)
  std::map<raw::ChannelID_t, float> amplitude_map;

  for (art::Ptr<recob::Hit> hit: muon_hits) {
    if (!amplitude_map.count(hit->Channel()) || amplitude_map.at(hit->Channel()) < hit->PeakAmplitude()) {
      amplitude_map[hit->Channel()] = hit->PeakAmplitude();
    }
  }

  // fill up all of the signals
  for (auto const &pair: amplitude_map) {
    fSignal[angle_bin][geo->View(pair.first)]->Fill(pair.second);
  }

  // now get the background noise RMS

  // figure out which channels do not have any energy depositions
  std::vector<raw::ChannelID_t> allChannels = geo->ChannelsInTPCs();

  std::set<raw::ChannelID_t> channelSet;
  for (raw::ChannelID_t c: allChannels) channelSet.insert(c);

  // ignore wires that are shorter than the max length
  for (raw::ChannelID_t c: allChannels) {
    geo::WireID wireID = geo->ChannelToWire(c).at(0);
    const geo::WireGeo &wireGeo = geo->Wire(wireID);
    float length = wireGeo.Length();
    // TODO: don't hardcode
    if (geo->SignalType(wireID) == geo::kInduction && length < 577.3) {
      channelSet.erase(c);
    }
  }

  art::Handle<std::vector<sim::SimChannel>> simChannels;
  ev.getByLabel("largeant", simChannels);
  for (const sim::SimChannel &c: *simChannels) {
    if (c.TDCIDEMap().size()) {
      channelSet.erase(c.Channel());
    }
  }

  // now get the wires
  art::Handle<std::vector<recob::Wire>> wires;
  ev.getByLabel(fWireLabel, wires);

  for (const recob::Wire &wire: *wires) {
    if (!channelSet.count(wire.Channel())) continue;
    // get the noise RMS
    const recob::Wire::RegionsOfInterest_t &ROI = wire.SignalROI();
    int ncount = 0;
    double rms_sum = 0.;
    for (lar::sparse_vector<float>::datarange_t range: ROI.get_ranges()) {
      for (float v: range) {
        rms_sum += v * v;
        ncount ++;
      }
    }
    float noise = sqrt(rms_sum / ncount);
    // std::cout << "Noise: " << noise;
    fNoise[geo->View(wire.Channel())]->Fill(noise);
  }

}

DEFINE_ART_MODULE(numu::MuonS2NStudy)
