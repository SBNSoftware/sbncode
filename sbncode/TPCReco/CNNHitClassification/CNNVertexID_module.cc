/////////////////////////////////////////////////////////////////////////////////////
// Class:       CNNVertexID
// Module Type: producer
// File:        CheckCNNScoreVertex_module.cc
//
// apply CNN to hits associated with cluster to get track/shower/noise/endVertex score
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
#include "sbnobj/Common/Reco/CNNVertexScore.h"

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sbn {

  class CNNVertexID : public art::EDProducer {
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
        Comment("tag of hits to be EM/track / Vertex tagged")};

      fhicl::Atom<art::InputTag> ClusterModuleLabel{
        Name("ClusterModuleLabel"),
        Comment("tag of clusters")};

      fhicl::Atom<art::InputTag> PFParticleModuleLabel{
        Name("PFParticleModuleLabel"),
        Comment("tag of PFParticles")};

      fhicl::Atom<bool> SkipClearCosmics{
        Name("SkipClearCosmics"),
        Comment("skip clear cosmic PFPs")};

      fhicl::Atom<bool> DoVertex{
        Name("DoVertex"),
        Comment("get Vertex Scores")};
    
      fhicl::Sequence<int> VertexRegionSize{
        Name("VertexRegionSize"),
        Comment("size of region around cluster end to get Vertex Scores")};

      fhicl::Sequence<int> Views{
        Name("Views"),
        Comment("tag clusters in selected views only, or in all views if empty list")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit CNNVertexID(Parameters const& p);

    CNNVertexID(CNNVertexID const&) = delete;
    CNNVertexID(CNNVertexID&&) = delete;
    CNNVertexID& operator=(CNNVertexID const&) = delete;
    CNNVertexID& operator=(CNNVertexID&&) = delete;

  private:
    void produce(art::Event& e) override;
    float getAvgScore(std::vector<float> scores);
    size_t fBatchSize;
    nnet::PointIdAlg fPointIdAlg;
    anab::MVAWriter<3> fMVAWriter; 
    art::InputTag fWireProducerLabel;
    art::InputTag fHitModuleLabel;
    art::InputTag fClusterModuleLabel;
    art::InputTag fPFParticleModuleLabel;
    bool fSkipClearCosmics;
    bool fDoVertex;
    std::vector<int> fVertexRegionSize;
    std::vector<int> fViews;

  };
  // ------------------------------------------------------

  CNNVertexID::CNNVertexID(CNNVertexID::Parameters const& config)
    : EDProducer{config}
    , fBatchSize(config().BatchSize())
    , fPointIdAlg(config().PointIdAlg())
    , fMVAWriter(producesCollector(), "cnnvertexid")
    , fWireProducerLabel(config().WireLabel())
    , fHitModuleLabel(config().HitModuleLabel())
    , fClusterModuleLabel(config().ClusterModuleLabel())
    , fPFParticleModuleLabel(config().PFParticleModuleLabel())
    , fSkipClearCosmics(config().SkipClearCosmics())
    , fDoVertex(config().DoVertex())
    , fVertexRegionSize(config().VertexRegionSize())
    , fViews(config().Views())
  {
    produces<std::vector<PFPCNNVertexScore>>();
    produces<art::Assns<recob::PFParticle, PFPCNNVertexScore>>();
    fMVAWriter.produces_using<recob::Hit>();
  }
  // ------------------------------------------------------

  void
  CNNVertexID::produce(art::Event& evt)
  {
    mf::LogVerbatim("CNNVertexID") << "next event: " << evt.run() << " / " << evt.id().event();

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

    auto pfpScores = std::make_unique<std::vector<PFPCNNVertexScore>>();
    auto pfpAssns = std::make_unique<art::Assns<recob::PFParticle, PFPCNNVertexScore>>();

    // loop over PFPs
    for (const art::Ptr<recob::PFParticle> &pfp : pfpList){
      std::vector<float> pfpEndVertexScore0;
      std::vector<float> pfpEndVertexScore1;
      std::vector<float> pfpEndVertexScore2;
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

        float cluEndVertexScore0 = 0.;
        float cluEndVertexScore1 = 0.;
        float cluEndVertexScore2 = 0.;

        if (fDoVertex) { // Vertex scores for hits around end of cluster
          int clueEndWire = cluster->EndWire();
          auto clueEndTick = cluster->EndTick();
          std::vector<art::Ptr<recob::Hit>> vertexPtrList;
          for (auto const& hit : hitPtrList) {
            unsigned int hitView = hit->WireID().Plane;
            int hitWireID = hit->WireID().Wire;
            auto hitPeakTime = hit->PeakTime();

            if (hitView != cluView) continue; // same plane
            if ((std::abs(hitWireID-clueEndWire) > fVertexRegionSize[0]) || (std::abs(hitPeakTime-clueEndTick) > fVertexRegionSize[1])) continue;
            bool isValidVertexHit = true;
            auto const &vertexClusters = cluFMhit.at(hit.key());
            for (auto const &vertexCluster : vertexClusters){
              // exclude if hit is on the parent cluster 
              if (vertexCluster.key() == cluster.key()) {
                isValidVertexHit = false;
                break;
              }
              // exclude if hit is on a cosmic cluster
              auto const &vertexPFPs = pfpFMclu.at(vertexCluster.key());
              for (auto const &vertexPFP : vertexPFPs){
                std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> vertexMetas = metaFMpfp.at(vertexPFP.key());
                auto const &vertexProperties = vertexMetas[0]->GetPropertiesMap();
                if (vertexProperties.count("IsClearCosmic")) {
                  isValidVertexHit = false;
                  break;
                }
              }
            }
            if (isValidVertexHit) vertexPtrList.push_back(hit);
          }
          std::vector<std::pair<unsigned int, float>> vertexHits;
          std::vector<size_t> keys;
          for (auto const& hit : vertexPtrList){
            vertexHits.emplace_back(hit->WireID().Wire, hit->PeakTime());
            keys.push_back(hit.key());
          }
          auto batchOut = fPointIdAlg.predictIdVectors(vertexHits);
          if (vertexHits.size() != batchOut.size()) {
              throw cet::exception("CNNVertexID") << "Size mismatch between input and output vectors";
          }
          size_t nVertexHits = vertexHits.size();
          for (size_t h = 0; h < nVertexHits; h++) {
            fMVAWriter.setOutput(hitID, keys[h], batchOut[h]);
            cluEndVertexScore0 += batchOut[h][0];
            cluEndVertexScore1 += batchOut[h][1];
            cluEndVertexScore2 += batchOut[h][2];
          }
          if (nVertexHits != 0) { 
            cluEndVertexScore0 = cluEndVertexScore0/nVertexHits;
            cluEndVertexScore1 = cluEndVertexScore1/nVertexHits;
            cluEndVertexScore2 = cluEndVertexScore2/nVertexHits;
          }
        }
        else {
          cluEndVertexScore0 = std::numeric_limits<float>::signaling_NaN();
          cluEndVertexScore1 = std::numeric_limits<float>::signaling_NaN();
          cluEndVertexScore2 = std::numeric_limits<float>::signaling_NaN();
        }
        // Vertex score end

        pfpEndVertexScore0.push_back(cluEndVertexScore0);
        pfpEndVertexScore1.push_back(cluEndVertexScore1);
        pfpEndVertexScore2.push_back(cluEndVertexScore2);
      }

      float pfpAvgEndVertexScore0 = getAvgScore(pfpEndVertexScore0);
      float pfpAvgEndVertexScore1 = getAvgScore(pfpEndVertexScore1);
      float pfpAvgEndVertexScore2 = getAvgScore(pfpEndVertexScore2);

      std::cout << "-----------" << " PFP " << pfp->Self() << "-----------" << std::endl;
      std::cout << " average end vertex score0: " << pfpAvgEndVertexScore0 << std::endl;
      std::cout << " average end vertex score1: " << pfpAvgEndVertexScore1 << std::endl;
      std::cout << " average end vertex score2: " << pfpAvgEndVertexScore2 << std::endl;
      std::cout << "------------------------------------------------------" << std::endl;
      pfpScores->emplace_back(pfpAvgEndVertexScore0, pfpAvgEndVertexScore1, pfpAvgEndVertexScore2, nClusters); //
      util::CreateAssn(*this, evt, *pfpScores, pfp, *pfpAssns);
    }
    
    fMVAWriter.saveOutputs(evt);
    evt.put(std::move(pfpScores));
    evt.put(std::move(pfpAssns));
  }
  // ------------------------------------------------------
  // get average score, drop outliers
  float CNNVertexID::getAvgScore(std::vector<float> scores){
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

  DEFINE_ART_MODULE(CNNVertexID)
}