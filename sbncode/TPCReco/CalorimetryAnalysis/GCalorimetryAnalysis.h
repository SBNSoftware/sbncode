#ifndef SBN_GCalorimetryAnalysis
#define SBN_GCalorimetryAnalysis

////////////////////////////////////////////////////////////////////////
// Class:       GCalorimetryAnalysis
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

#include "larcorealg/GeoAlgo/GeoAlgo.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "GCalorimetryAnalysisObj.h"
#include "IGCaloSelectionTool.h"

namespace sbn {
  class GCalorimetryAnalysis;
}

class sbn::GCalorimetryAnalysis : public art::EDAnalyzer {
public:
  explicit GCalorimetryAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GCalorimetryAnalysis(GCalorimetryAnalysis const&) = delete;
  GCalorimetryAnalysis(GCalorimetryAnalysis&&) = delete;
  GCalorimetryAnalysis& operator=(GCalorimetryAnalysis const&) = delete;
  GCalorimetryAnalysis& operator=(GCalorimetryAnalysis&&) = delete;

  ~GCalorimetryAnalysis();

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
  };

  // Fill vars
  void FillTrack(const recob::Track &track, 
    const recob::PFParticle &pfp, float t0, 
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<const recob::TrackHitMeta*> &thms,
    const std::vector<art::Ptr<anab::Calorimetry>> &calo,
    const std::vector<GlobalTrackInfo> &tracks);

  HitInfo MakeHit(const recob::Hit &hit,
    unsigned hkey,
    const recob::TrackHitMeta &thm,
    const recob::Track &trk,
    const std::vector<art::Ptr<anab::Calorimetry>> &calo);

  // config

  // tags
  art::InputTag fPFPproducer;
  art::InputTag fT0Producer;
  art::InputTag fCALOproducer;
  art::InputTag fTRKproducer;
  bool fRequireT0;

  // tools
  std::vector<std::unique_ptr<sbn::IGCaloSelectionTool>> fSelectionTools;

  // persistent info
  MetaInfo fMeta;

  // Output
  TTree *fTree;
  TrackInfo *fTrack;

  // Fitting info
  TFitter fFitExp;
  TFitter fFitConst;
};



#endif
