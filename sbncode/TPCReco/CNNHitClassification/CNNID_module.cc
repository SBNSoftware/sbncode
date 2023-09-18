/////////////////////////////////////////////////////////////////////////////////////
// Class:       CNNID
// Module Type: producer
// File:        CNNID_module.cc
//
// apply CNN to hits associated with cluster to get track/shower/noise/endMichel score
// score of a PFP is the average score of all associated hits
// skip clear cosmic PFPs
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
    mf::LogVerbatim("CNNID") << "next event: " << evt.run() << " / " << evt.id().event();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

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

    auto hitID = fMVAWriter.initOutputs<recob::Hit>(
      fHitModuleLabel, hitPtrList.size(), fPointIdAlg.outputLabels());

    auto pfpScores = std::make_unique<std::vector<PFPCNNScore>>();
    auto pfpAssns = std::make_unique<art::Assns<recob::PFParticle, PFPCNNScore>>();

    // loop over PFPs
    for (const art::Ptr<recob::PFParticle> &pfp : pfpList){
      std::vector<float> pfpTrackScore; 
      std::vector<float> pfpShowerScore; 
      std::vector<float> pfpNoiseScore;
      std::vector<float> pfpEndMichelScore;
      int nClusters = 0; 

      if (fSkipClearCosmics) { // skip clear cosmic PFPs
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metas = metaFMpfp.at(pfp.key());
        auto const &properties = metas[0]->GetPropertiesMap();
        if (properties.count("IsClearCosmic")) continue;
      }

      // loop over clusters
      auto const &clusters = cluFMpfp.at(pfp->Self());
      for (auto const &cluster : clusters){
        nClusters++;
        unsigned int cluView = cluster->Plane().Plane;
        unsigned int cluCryo = cluster->Plane().Cryostat;
        unsigned int cluTPC = cluster->Plane().TPC;
        fPointIdAlg.setWireDriftData(clockData, detProp, *wireHandle, cluView, cluTPC, cluCryo);

        float cluTrackScore = 0.;
        float cluShowerScore = 0.;
        float cluNoiseScore = 0.;
        float cluEndMichelScore = 0.;

        if (fDoMichel) { // Michel scores for hits around end of cluster
          int clueEndWire = cluster->EndWire();
          auto clueEndTick = cluster->EndTick();
          std::vector<art::Ptr<recob::Hit>> michelPtrList;
          for (auto const& hit : hitPtrList) {
            unsigned int hitView = hit->WireID().Plane;
            int hitWireID = hit->WireID().Wire;
            auto hitPeakTime = hit->PeakTime();

            if (hitView != cluView) continue; // same plane
            if ((std::abs(hitWireID-clueEndWire) > fMichelRegionSize[0]) || (std::abs(hitPeakTime-clueEndTick) > fMichelRegionSize[1])) continue;
            bool isValidMichelHit = true;
            auto const &michelClusters = cluFMhit.at(hit.key());
            for (auto const &michelCluster : michelClusters){
              // exclude if hit is on the parent cluster 
              if (michelCluster.key() == cluster.key()) {
                isValidMichelHit = false;
                break;
              }
              // exclude if hit is on a cosmic cluster
              auto const &michelPFPs = pfpFMclu.at(michelCluster.key());
              for (auto const &michelPFP : michelPFPs){
                std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> michelMetas = metaFMpfp.at(michelPFP.key());
                auto const &michelProperties = michelMetas[0]->GetPropertiesMap();
                if (michelProperties.count("IsClearCosmic")) {
                  isValidMichelHit = false;
                  break;
                }
              }
            }
            if (isValidMichelHit) michelPtrList.push_back(hit);
          }
          std::vector<std::pair<unsigned int, float>> michelHits;
          std::vector<size_t> keys;
          for (auto const& hit : michelPtrList){
            michelHits.emplace_back(hit->WireID().Wire, hit->PeakTime());
            keys.push_back(hit.key());
          }
          auto batchOut = fPointIdAlg.predictIdVectors(michelHits);
          if (michelHits.size() != batchOut.size()) {
              throw cet::exception("CNNID") << "Size mismatch between input and output vectors";
          }
          size_t nMichelHits = michelHits.size();
          for (size_t h = 0; h < nMichelHits; h++) {
            fMVAWriter.setOutput(hitID, keys[h], batchOut[h]);
            cluEndMichelScore += batchOut[h][3];
          }
          if (nMichelHits != 0) cluEndMichelScore = cluEndMichelScore/nMichelHits;
        }
        else cluEndMichelScore = std::numeric_limits<float>::signaling_NaN();
        // Michel score end

        // non clear cosmic pfp track/shower/noise score
        if (fDoPFP) {
          auto &cluHitPtrList = hitFMclu.at(cluster.key());
          std::vector<art::Ptr<recob::Hit>> cluHitPtrListApply; 
          int nCluHits = cluHitPtrList.size();
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
          size_t nCluApplyHits = cluApplyHits.size();
          auto batchOut = fPointIdAlg.predictIdVectors(cluApplyHits);
          if (nCluApplyHits != batchOut.size()) {
            throw cet::exception("CNNID") << "Size mismatch between input and output vectors";
          }
          for (size_t h = 0; h < nCluApplyHits; h++) {
            fMVAWriter.setOutput(hitID, keys[h], batchOut[h]);
            cluTrackScore += batchOut[h][0];
            cluShowerScore += batchOut[h][1];
            cluNoiseScore += batchOut[h][2];
          }
          if (nCluApplyHits !=0) {
            cluTrackScore = cluTrackScore/nCluApplyHits;
            cluShowerScore = cluShowerScore/nCluApplyHits;
            cluNoiseScore = cluNoiseScore/nCluApplyHits;
          }
        }
        else {
          cluTrackScore = std::numeric_limits<float>::signaling_NaN();
          cluShowerScore = std::numeric_limits<float>::signaling_NaN();
          cluNoiseScore = std::numeric_limits<float>::signaling_NaN();
        }
        // cluster scores end

        pfpTrackScore.push_back(cluTrackScore);
        pfpShowerScore.push_back(cluShowerScore);
        pfpNoiseScore.push_back(cluNoiseScore);
        pfpEndMichelScore.push_back(cluEndMichelScore);
      }

      float pfpAvgTrackScore = getAvgScore(pfpTrackScore);
      float pfpAvgShowerScore = getAvgScore(pfpShowerScore);
      float pfpAvgNoiseScore = getAvgScore(pfpNoiseScore);
      float pfpAvgEndMichelScore = getAvgScore(pfpEndMichelScore);

    //   std::cout << "-----------" << " PFP " << pfp->Self() << "-----------" << std::endl;
    //   std::cout << " average track score: " << pfpAvgTrackScore << std::endl;
    //   std::cout << " average shower score: " << pfpAvgShowerScore << std::endl;
    //   std::cout << " average noise score: " << pfpAvgNoiseScore << std::endl;
    //   std::cout << " average michel score: " << pfpAvgEndMichelScore << std::endl;
    //   std::cout << "------------------------------------------------------" << std::endl;
      pfpScores->emplace_back(pfpAvgTrackScore, pfpAvgShowerScore, pfpAvgNoiseScore, pfpAvgEndMichelScore, nClusters); //
      util::CreateAssn(*this, evt, *pfpScores, pfp, *pfpAssns);
    }
    
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