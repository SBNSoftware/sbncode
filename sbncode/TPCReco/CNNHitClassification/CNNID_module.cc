/////////////////////////////////////////////////////////////////////////////////////
// Class:       CNNID
// Module Type: producer
// File:        CNNID_module.cc
//
// apply CNN to hits associated with cluster to get track/shower/noise/endMichel score
// score of a PFP is the average score of all associated hits
// skip clear cosmic PFPs
// largely based on larrecodnn/ImagePatternAlgs/Tensorflow/Modules/EmTrackMichelId_module.cc
//
// M.Jung munjung@uchicago.edu
/////////////////////////////////////////////////////////////////////////////////////

#include "larrecodnn/ImagePatternAlgs/Tensorflow/PointIdAlg/PointIdAlg.h"
#include "lardata/ArtDataHelper/MVAWriter.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "sbnobj/Common/Reco/CNNScore.h"

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sbn {

  class CNNID : public art::EDProducer {
  public:

    typedef std::unordered_map<unsigned int, std::vector<size_t>> view_keymap;
    typedef std::unordered_map<unsigned int, view_keymap> tpc_view_keymap;
    typedef std::unordered_map<unsigned int, tpc_view_keymap> cryo_tpc_view_keymap;

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<nnet::PointIdAlg::Config> PointIdAlg{Name("PointIdAlg")};
      fhicl::Atom<size_t> BatchSize{Name("BatchSize"),
                                    Comment("number of samples processed in one batch")};

      fhicl::Atom<art::InputTag> WireLabel{
        Name("WireLabel"),
        Comment("tag of deconvoluted ADC on wires (recob::Wire)")};

      fhicl::Atom<art::InputTag> HitModuleLabel{
        Name("HitModuleLabel"),
        Comment("tag of hits to be EM/track / Michel tagged")};

      fhicl::Atom<art::InputTag> ClusterModuleLabel{
        Name("ClusterModuleLabel"),
        Comment("tag of clusters")};

      fhicl::Atom<art::InputTag> PFParticleModuleLabel{
        Name("PFParticleModuleLabel"),
        Comment("tag of PFParticles")};

      fhicl::Atom<bool> SkipClearCosmics{
        Name("SkipClearCosmics"),
        Comment("skip clear cosmic PFPs")};

      fhicl::Atom<bool> DoMichel{
        Name("DoMichel"),
        Comment("get Michel Scores")};
    
      fhicl::Sequence<int> MichelRegionSize{
        Name("MichelRegionSize"),
        Comment("size of region around cluster end to get Michel Scores")};

      fhicl::Atom<bool> DoPFP{
        Name("DoPFP"),
        Comment("get PFP Scores")};

      fhicl::Atom<int> SparseLengthCut{
        Name("SparseLengthCut"),
        Comment("sparsely apply CNN if the number of hits is greater than this number")};

      fhicl::Atom<int> SparseRate{
        Name("SparseRate"),
        Comment("rate to sparsely apply CNN")};

      fhicl::Sequence<int> Views{
        Name("Views"),
        Comment("tag clusters in selected views only, or in all views if empty list")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit CNNID(Parameters const& p);

    CNNID(CNNID const&) = delete;
    CNNID(CNNID&&) = delete;
    CNNID& operator=(CNNID const&) = delete;
    CNNID& operator=(CNNID&&) = delete;

  private:
    void produce(art::Event& e) override;
    float getAvgScore(std::vector<float> scores);
    size_t fBatchSize;
    nnet::PointIdAlg fPointIdAlg;
    anab::MVAWriter<4> fMVAWriter; 
    art::InputTag fWireProducerLabel;
    art::InputTag fHitModuleLabel;
    art::InputTag fClusterModuleLabel;
    art::InputTag fPFParticleModuleLabel;
    bool fSkipClearCosmics;
    bool fDoMichel;
    std::vector<int> fMichelRegionSize;
    bool fDoPFP;
    int fSparseLengthCut;
    int fSparseRate;
    std::vector<int> fViews;

  };
  // ------------------------------------------------------

  CNNID::CNNID(CNNID::Parameters const& config)
    : EDProducer{config}
    , fBatchSize(config().BatchSize())
    , fPointIdAlg(config().PointIdAlg())
    , fMVAWriter(producesCollector(), "cnnid")
    , fWireProducerLabel(config().WireLabel())
    , fHitModuleLabel(config().HitModuleLabel())
    , fClusterModuleLabel(config().ClusterModuleLabel())
    , fPFParticleModuleLabel(config().PFParticleModuleLabel())
    , fSkipClearCosmics(config().SkipClearCosmics())
    , fDoMichel(config().DoMichel())
    , fMichelRegionSize(config().MichelRegionSize())
    , fDoPFP(config().DoPFP())
    , fSparseLengthCut(config().SparseLengthCut())
    , fSparseRate(config().SparseRate())
    , fViews(config().Views())
  {
    produces<std::vector<PFPCNNScore>>();
    produces<art::Assns<recob::PFParticle, PFPCNNScore>>();
    fMVAWriter.produces_using<recob::Hit>();
  }
  // ------------------------------------------------------

  void
  CNNID::produce(art::Event& evt)
  {
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    unsigned int cryo, tpc, view;
    CNNID::cryo_tpc_view_keymap cluHitMap;
    CNNID::cryo_tpc_view_keymap endMichelHitMap;

    auto wireHandle = evt.getValidHandle<std::vector<recob::Wire>>(fWireProducerLabel);

    auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);
    std::vector<art::Ptr<recob::Hit>> hitPtrList;
    art::fill_ptr_vector(hitPtrList, hitListHandle);

    auto pfpListHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleModuleLabel);
	std::vector<art::Ptr<recob::PFParticle>> pfpList;
    art::fill_ptr_vector(pfpList, pfpListHandle);

    auto metaHandle = evt.getValidHandle<std::vector<larpandoraobj::PFParticleMetadata>>(fPFParticleModuleLabel);
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metaPtrList;
    art::fill_ptr_vector(metaPtrList, metaHandle);

    auto cluListHandle = evt.getValidHandle<std::vector<recob::Cluster>>(fClusterModuleLabel);
    std::vector<art::Ptr<recob::Cluster>> cluPtrList;
    art::fill_ptr_vector(cluPtrList, cluListHandle);

    art::FindManyP<larpandoraobj::PFParticleMetadata> metaFMpfp(pfpListHandle, evt, fPFParticleModuleLabel);
    art::FindManyP<recob::Cluster> cluFMpfp(pfpListHandle, evt, fPFParticleModuleLabel);
    art::FindManyP<recob::PFParticle> pfpFMclu(cluListHandle, evt, fClusterModuleLabel);
    art::FindManyP<recob::Hit> hitFMclu(cluListHandle, evt, fClusterModuleLabel);
    art::FindManyP<recob::Hit> hitFMpfp(pfpListHandle, evt, fPFParticleModuleLabel);
    art::FindManyP<recob::Cluster> cluFMhit(hitListHandle, evt, fClusterModuleLabel);

    auto hitID = fMVAWriter.initOutputs<recob::Hit>(fHitModuleLabel, hitPtrList.size(), fPointIdAlg.outputLabels());

    auto pfpScores = std::make_unique<std::vector<PFPCNNScore>>();
    auto pfpAssns = std::make_unique<art::Assns<recob::PFParticle, PFPCNNScore>>();

    // loop over PFPs
    for (const art::Ptr<recob::PFParticle> &pfp : pfpList){
      int nClusters = 0; 
      std::vector<float> pfpTrackScore; 
      std::vector<float> pfpShowerScore; 
      std::vector<float> pfpNoiseScore;
      std::vector<float> pfpMichelScore;
      std::vector<float> pfpEndMichelScore;
      float pfpAvgTrackScore;
      float pfpAvgShowerScore;
      float pfpAvgNoiseScore;
      float pfpAvgMichelScore;
      float pfpAvgEndMichelScore;

      if (fSkipClearCosmics) {
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metas = metaFMpfp.at(pfp.key());
        auto const &properties = metas[0]->GetPropertiesMap();
        if (properties.count("IsClearCosmic")) {
           pfpAvgTrackScore = std::numeric_limits<float>::signaling_NaN();
           pfpAvgShowerScore = std::numeric_limits<float>::signaling_NaN();
           pfpAvgNoiseScore = std::numeric_limits<float>::signaling_NaN();
           pfpAvgMichelScore = std::numeric_limits<float>::signaling_NaN();
           pfpAvgEndMichelScore = std::numeric_limits<float>::signaling_NaN();
           pfpScores->emplace_back(pfpAvgTrackScore, pfpAvgShowerScore, pfpAvgNoiseScore, pfpAvgMichelScore, pfpAvgEndMichelScore, nClusters); //
           util::CreateAssn(*this, evt, *pfpScores, pfp, *pfpAssns);
           continue;
        }
      }

      auto const &clusters = cluFMpfp.at(pfp->Self());
      for (auto const &cluster : clusters){
        cluHitMap.clear();
        endMichelHitMap.clear();
        float cluTrackScore = 0.;
        float cluShowerScore = 0.;
        float cluNoiseScore = 0.;
        float cluMichelScore = 0.;
        float cluEndMichelScore = 0.;

        nClusters++;
        view = cluster->Plane().Plane;
        cryo = cluster->Plane().Cryostat;
        tpc = cluster->Plane().TPC;

        // Michel scores for hits around end of cluster
        if (fDoMichel) {
          std::vector<art::Ptr<recob::Hit>> michelPtrList;
          int clueEndWire = cluster->EndWire();
          auto clueEndTick = cluster->EndTick();

          // find hits around end of cluster
          for (auto const& hit : hitPtrList) {
            unsigned int hitView = hit->WireID().Plane;
            int hitWireID = hit->WireID().Wire;
            auto hitPeakTime = hit->PeakTime();
            if (hitView != view) continue;  // hit should be on the same plane
            if ((std::abs(hitWireID-clueEndWire) > fMichelRegionSize[0]) || (std::abs(hitPeakTime-clueEndTick) > fMichelRegionSize[1])) continue;
            michelPtrList.push_back(hit);
          }
          for (auto const& hit : michelPtrList){
            view = hit->WireID().Plane;
            cryo = hit->WireID().Cryostat;
            tpc = hit->WireID().TPC;
            endMichelHitMap[cryo][tpc][view].push_back(hit.key());
          }
          // run CNN over hits
          int nEndMichelHits = 0;
          for (auto const& pcryo : endMichelHitMap) {
            cryo = pcryo.first;
            for (auto const& ptpc : pcryo.second) {
              tpc = ptpc.first;
              for (auto const& pview : ptpc.second) {
                view = pview.first;
                fPointIdAlg.setWireDriftData(clockData, detProp, *wireHandle, view, tpc, cryo);
                for (size_t idx = 0; idx < pview.second.size(); idx += fBatchSize) {
                  std::vector<std::pair<unsigned int, float>> points;
                  std::vector<size_t> keys;
                  for (size_t k = 0; k < fBatchSize; ++k) {
                    if (idx + k >= pview.second.size()) { break; } // careful about the tail
                    size_t h = pview.second[idx + k]; // h is the Ptr< recob::Hit >::key()
                    const recob::Hit& hit = *(hitPtrList[h]);
                    points.emplace_back(hit.WireID().Wire, hit.PeakTime());
                    keys.push_back(h);
                  }
                  auto batch_out = fPointIdAlg.predictIdVectors(points);
                  if (points.size() != batch_out.size()) {
                    throw cet::exception("CNNID") << "hits processing failed" << std::endl;
                  }
                  for (size_t k = 0; k < points.size(); ++k) {
                    size_t h = keys[k];
                    fMVAWriter.setOutput(hitID, h, batch_out[k]);
                    nEndMichelHits++;
                    cluEndMichelScore += batch_out[k][3];
                  }
                }
              }
            }
          }  // run CNN end
          if (nEndMichelHits != 0) cluEndMichelScore = cluEndMichelScore/nEndMichelHits;
          else cluEndMichelScore = std::numeric_limits<float>::signaling_NaN();
        }  // cluster end Michel score end

        // pfp scores
        if (fDoPFP) {
          auto &cluHitPtrList = hitFMclu.at(cluster.key());
          std::vector<art::Ptr<recob::Hit>> cluHitPtrListApply; 
          int nCluHits = cluHitPtrList.size();
          // skip if no hits
          if (nCluHits == 0) {
            cluTrackScore = std::numeric_limits<float>::signaling_NaN();
            cluShowerScore = std::numeric_limits<float>::signaling_NaN();
            cluNoiseScore = std::numeric_limits<float>::signaling_NaN();
            cluMichelScore = std::numeric_limits<float>::signaling_NaN();
            continue;
          }
          // collect cluster hits
          if (nCluHits > fSparseLengthCut) {
            for (int i = 0; i < nCluHits; i++){
              if (i % fSparseRate != 0) continue;
              cluHitPtrListApply.push_back(cluHitPtrList[i]);
            }
          }
          else cluHitPtrListApply = cluHitPtrList;
          std::vector<std::pair<unsigned int, float>> cluApplyHits;
          std::vector<size_t> keys;
          for (auto const& hit : cluHitPtrListApply){
            cluApplyHits.emplace_back(hit->WireID().Wire, hit->PeakTime());
            keys.push_back(hit.key());
          }
          for (auto const& h : cluHitPtrListApply) {
            view = h->WireID().Plane;
            cryo = h->WireID().Cryostat;
            tpc = h->WireID().TPC;
            cluHitMap[cryo][tpc][view].push_back(h.key());
          }
          // run CNN over hits
          int nCluApplyHits = 0;
          for (auto const& pcryo : cluHitMap) {
            cryo = pcryo.first;
            for (auto const& ptpc : pcryo.second) {
              tpc = ptpc.first;
              for (auto const& pview : ptpc.second) {
                view = pview.first;
                fPointIdAlg.setWireDriftData(clockData, detProp, *wireHandle, view, tpc, cryo);
                for (size_t idx = 0; idx < pview.second.size(); idx += fBatchSize) {
                  std::vector<std::pair<unsigned int, float>> points;
                  std::vector<size_t> keys;
                  for (size_t k = 0; k < fBatchSize; ++k) {
                    if (idx + k >= pview.second.size()) { break; }
                    size_t h = pview.second[idx + k];
                    const recob::Hit& hit = *(hitPtrList[h]);
                    points.emplace_back(hit.WireID().Wire, hit.PeakTime());
                    keys.push_back(h);
                  }
                  auto batch_out = fPointIdAlg.predictIdVectors(points);
                  if (points.size() != batch_out.size()) {
                    throw cet::exception("CNNID") << "hits processing failed" << std::endl;
                  }
                  for (size_t k = 0; k < points.size(); ++k) {
                    size_t h = keys[k];
                    fMVAWriter.setOutput(hitID, h, batch_out[k]);
                    nCluApplyHits++;
                    cluTrackScore += batch_out[k][0];
                    cluShowerScore += batch_out[k][1];
                    cluNoiseScore += batch_out[k][2];
                    cluMichelScore += batch_out[k][3];
                  }
                }
              }
            }
          } // run CNN end
          pfpTrackScore.push_back(cluTrackScore);
          pfpShowerScore.push_back(cluShowerScore);
          pfpNoiseScore.push_back(cluNoiseScore);
          pfpMichelScore.push_back(cluMichelScore);
          pfpEndMichelScore.push_back(cluEndMichelScore);
        }  // cluster score end 
      }  // cluster loop end
      pfpAvgTrackScore = getAvgScore(pfpTrackScore);
      pfpAvgShowerScore = getAvgScore(pfpShowerScore);
      pfpAvgNoiseScore = getAvgScore(pfpNoiseScore);
      pfpAvgMichelScore = getAvgScore(pfpMichelScore);
      pfpAvgEndMichelScore = getAvgScore(pfpEndMichelScore);
      pfpScores->emplace_back(pfpAvgTrackScore, pfpAvgShowerScore, pfpAvgNoiseScore, pfpAvgMichelScore, pfpAvgEndMichelScore, nClusters); //
      util::CreateAssn(*this, evt, *pfpScores, pfp, *pfpAssns);
    }  // pfp score end
    fMVAWriter.saveOutputs(evt);
    evt.put(std::move(pfpScores));
    evt.put(std::move(pfpAssns));
  }

  // ------------------------------------------------------
  // get average score, drop outliers
  float CNNID::getAvgScore(std::vector<float> scores){
    // good score means a positive score (default score is set to 0)
    // if no good score, return NaN
    // if only 1 good score, use that score as "average"
    // if 2 good scores, use the average of the 2 scores
    // if 3 or more good scores, remove outliers if any and use the average of the remaining scores
    int nScores = scores.size();
    if (nScores == 0) return std::numeric_limits<float>::signaling_NaN();
    int nGoodScore = 0;
    std::vector<float> goodScores;
    for (int i = 0; i < nScores; i++) {
      if (scores[i] == 0.) continue;
      goodScores.push_back(scores[i]);
      nGoodScore++;
    }
    if (nGoodScore == 0) return std::numeric_limits<float>::signaling_NaN();
    else if (nGoodScore == 1) return goodScores[0];
    else if (nGoodScore == 2) return (goodScores[0]+goodScores[1])/2.;
    else {
      std::sort(goodScores.begin(), goodScores.end());
      float diffLow = goodScores[1] - goodScores[0];
      float diffHigh = goodScores[goodScores.size()-1] - goodScores[goodScores.size()-2];
      if (diffLow > 0.5) return (goodScores[1]+goodScores[2])/2.; // remove lowest score (outlier)
      else if (diffHigh > 0.5) return (goodScores[0]+goodScores[1])/2.; // remove highest score (outlier)
      else return (goodScores[0]+goodScores[1]+goodScores[2])/3.;
    }
  }
  // ------------------------------------------------------

  DEFINE_ART_MODULE(CNNID)
}
