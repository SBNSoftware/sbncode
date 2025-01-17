#ifndef SBN_TrackCaloSkimmer
#define SBN_TrackCaloSkimmer

////////////////////////////////////////////////////////////////////////
// Class:       TrackCaloSkimmer
// Plugin Type: analyzer (art v3_06_03)
//
// Generated at Mon May 17 09:46:34 2021 by Gray Putnam using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////


#include <cstdlib>
#include <iostream>

#include "TTree.h"
#include "TFitter.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "lardataobj/RawData/RawDigit.h"

#include "larcorealg/GeoAlgo/GeoAlgo.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "sbnobj/Common/Calibration/TrackCaloSkimmerObj.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "ITCSSelectionTool.h"

namespace sbn {
  class TrackCaloSkimmer;
  enum EDet {kNOTDEFINED, kSBND, kICARUS}; 
}

class sbn::TrackCaloSkimmer : public art::EDAnalyzer {
public:
  explicit TrackCaloSkimmer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackCaloSkimmer(TrackCaloSkimmer const&) = delete;
  TrackCaloSkimmer(TrackCaloSkimmer&&) = delete;
  TrackCaloSkimmer& operator=(TrackCaloSkimmer const&) = delete;
  TrackCaloSkimmer& operator=(TrackCaloSkimmer&&) = delete;

  ~TrackCaloSkimmer();

  void analyze(art::Event const& e) override;

  void respondToOpenInputFile(const art::FileBlock& fb) override {
    (void) fb;
    fMeta.ifile ++;
  }

private:
  // Internal data struct
  struct GlobalTrackInfo {
    geo::Point_t start;
    geo::Point_t end;
    geo::Vector_t dir;
    geo::Vector_t enddir;
    int ID;
  };


  // Represents a "Snippet" of ADCs shared by a set of hits on a wire
  struct Snippet {
    geo::WireID wire;
    int start_tick;
    int end_tick;

    bool operator==(const Snippet &rhs) const {
      return wire == rhs.wire && start_tick == rhs.start_tick && end_tick == rhs.end_tick;
    }

    bool operator<(const Snippet &rhs) const {
      return wire < rhs.wire || (wire == rhs.wire && start_tick < rhs.start_tick);
    }
  };

  // Fill vars
  void FillTrack(const recob::Track &track, 
    const recob::PFParticle &pfp, float t0, 
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<const recob::TrackHitMeta*> &thms,
    const std::vector<art::Ptr<recob::SpacePoint>> &sps,
    const std::vector<art::Ptr<anab::Calorimetry>> &calo,
    const art::Ptr<recob::MCSFitResult> &mcsUF,
    const art::Ptr<recob::MCSFitResult> &mcsUI,
    const art::Ptr<recob::MCSFitResult> &mcsIF3D,
    const art::Ptr<recob::MCSFitResult> &mcsIF2DI1,
    const art::Ptr<recob::MCSFitResult> &mcsIF2DI2,
    const art::Ptr<recob::MCSFitResult> &mcsIF2DC,
    const art::Ptr<recob::MCSFitResult> &mcsII3D,
    const art::Ptr<recob::MCSFitResult> &mcsII2DI1,
    const art::Ptr<recob::MCSFitResult> &mcsII2DI2,
    const art::Ptr<recob::MCSFitResult> &mcsII2DC,
    const art::Ptr<sbn::RangeP> &range,
    const std::map<geo::WireID, art::Ptr<raw::RawDigit>> &rawdigits,
    const std::vector<GlobalTrackInfo> &tracks,
    const geo::GeometryCore *geo,
    const detinfo::DetectorClocksData &clock_data,
    const cheat::BackTrackerService *bt_serv,
    const sbn::EDet det,
    const int whichcrt,
    const float crtt0,
    const float confidence);

  void FillTrackDaughterRays(const recob::Track &trk,
    const recob::PFParticle &pfp, 
    const std::vector<art::Ptr<recob::PFParticle>> &PFParticleList, 
    const art::FindManyP<recob::SpacePoint> &PFParticleSPs);

  void FillTrackEndHits(const geo::GeometryCore *geometry,
    const detinfo::DetectorPropertiesData &dprop,
    const recob::Track &track,
    const std::vector<art::Ptr<recob::Hit>> &allHits,
    const art::FindManyP<recob::SpacePoint> &allHitSPs);

  void FillTrackTruth(const detinfo::DetectorClocksData &clock_data,
    const std::vector<art::Ptr<recob::Hit>> &trkHits,
    const std::vector<art::Ptr<simb::MCParticle>> &mcparticles,
    const std::vector<geo::BoxBoundedGeo> &active_volumes,
    const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
    const std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map,
    const std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map,
    const detinfo::DetectorPropertiesData &dprop,
    const geo::GeometryCore *geo);

  TrackHitInfo MakeHit(const recob::Hit &hit,
    unsigned hkey,
    const recob::TrackHitMeta &thm,
    const recob::Track &trk,
    const art::Ptr<recob::SpacePoint> &sp,
    const std::vector<art::Ptr<anab::Calorimetry>> &calo,
    const geo::GeometryCore *geo,
    const detinfo::DetectorClocksData &dclock,
    const cheat::BackTrackerService *bt_serv);

  void DoTailFit();

  // config

  // tags
  art::InputTag fPFPproducer;

  art::InputTag fMCSproducerUF;
  art::InputTag fMCSproducerUI;
  art::InputTag fMCSproducerIF3D;
  art::InputTag fMCSproducerIF2DI1;
  art::InputTag fMCSproducerIF2DI2;
  art::InputTag fMCSproducerIF2DC;
  art::InputTag fMCSproducerII3D;
  art::InputTag fMCSproducerII2DI1;
  art::InputTag fMCSproducerII2DI2;
  art::InputTag fMCSproducerII2DC;

  std::vector<art::InputTag> fproducerCRTTag;
  art::InputTag fRangeInputtag;

  std::vector<art::InputTag> fT0producers;
  art::InputTag fCALOproducer;
  art::InputTag fTRKproducer;
  art::InputTag fTRKHMproducer;
  art::InputTag fHITproducer;
  std::vector<art::InputTag> fRawDigitproducers;
  std::string fG4producer;
  std::string fSimChannelproducer;
  bool fRequireT0;
  bool fDoTailFit;
  bool fVerbose;
  bool fSilenceMissingDataProducts;
  double fHitRawDigitsTickCollectWidth;
  double fTailFitResidualRange;
  int fHitRawDigitsWireCollectWidth;
  bool fFillTrackEndHits;
  float fTrackEndHitWireBox;
  float fTrackEndHitTimeBox;

  // tools
  std::vector<std::unique_ptr<sbn::ITCSSelectionTool>> fSelectionTools;

  // persistent info
  MetaInfo fMeta;

  std::map<Snippet, int> fSnippetCount;
  std::map<geo::WireID, std::pair<int, int>> fWiresToSave;

  // Output
  TTree *fTree;
  TrackInfo *fTrack;

  // Fitting info
  TFitter fFitExp;
  TFitter fFitConst;
};
#endif
