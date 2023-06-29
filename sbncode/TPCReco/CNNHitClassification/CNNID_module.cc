/////////////////////////////////////////////////////////////////////////////////////
// Class:       CNNID
// Module Type: producer
// File:        CNNID_module.cc
//
// apply CNN to hits associated with cluster to get track/shower/noise score
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

    unsigned int cryo, tpc, view;

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

    art::FindManyP<recob::Cluster> cluFMpfp(pfpListHandle, evt, fClusterModuleLabel);

    art::FindManyP<recob::Hit> hitFMclu(cluListHandle, evt, fClusterModuleLabel);
    art::FindManyP<recob::Hit> hitFMpfp(pfpListHandle, evt, fPFParticleModuleLabel);

    auto hitID = fMVAWriter.initOutputs<recob::Hit>(
      fHitModuleLabel, hitPtrList.size(), fPointIdAlg.outputLabels());

    auto pfpScores = std::make_unique<std::vector<PFPCNNScore>>();
    auto pfpAssns = std::make_unique<art::Assns<recob::PFParticle, PFPCNNScore>>();

    // loop over PFPs
    for (const art::Ptr<recob::PFParticle> &thispfp : pfpList){
      std::vector<float> PFPTrackScore; 
      std::vector<float> PFPShowerScore; 
      std::vector<float> PFPNoiseScore;
      // std::vector<float> PFPMichelScore;
      int nClusters = 0; 

      // skip clear cosmic PFPs
      std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metas = metaFMpfp.at(thispfp.key());
      const art::Ptr<larpandoraobj::PFParticleMetadata> meta = metas[0];
      auto const &properties = meta->GetPropertiesMap();
      bool is_clear_cosmic = (properties.count("IsClearCosmic")); 
      if (is_clear_cosmic) continue;

      // loop over clusters
      auto const &clusters = cluFMpfp.at(thispfp->Self());
      for (auto const &cluster : clusters){
        float cluTrackScore = 0.;
        float cluShowerScore = 0.;
        float cluNoiseScore = 0.;
        // float cluMichelScore = 0.;

        auto &cluHitPtrList = hitFMclu.at(cluster.key());
        int nHitsClu = cluHitPtrList.size();
        if (nHitsClu < 5) continue; // at least 5 hits in cluster
        nClusters++;

        cryo = cluster->Plane().Cryostat;
        view = cluster->Plane().Plane;
        tpc = cluster->Plane().TPC;
        fPointIdAlg.setWireDriftData(clockData, detProp, *wireHandle, view, tpc, cryo);

        // hits to apply CNN to (sparsely if nHitsClu > fSparseLengthCut)
        std::vector<art::Ptr<recob::Hit>> cluHitPtrListApply; 
        if (nHitsClu > fSparseLengthCut) {
          for (int i = 0; i < nHitsClu; i++){
            if (i % fSparseRate != 0) continue;
            cluHitPtrListApply.push_back(cluHitPtrList[i]);
          }
        }
        else {
          cluHitPtrListApply = cluHitPtrList;
        }

        // apply CNN
        std::vector<std::pair<unsigned int, float>> points;
        std::vector<size_t> keys;
        for (auto const& hit : cluHitPtrListApply){
          points.emplace_back(hit->WireID().Wire, hit->PeakTime());
          keys.push_back(hit.key());
        }

        auto batch_out = fPointIdAlg.predictIdVectors(points);
        if (points.size() != batch_out.size()) {
          throw cet::exception("CNNID") << "Size mismatch between input and output vectors";
        }

        size_t nhits = points.size();
        for (size_t h = 0; h < nhits; h++) {
          fMVAWriter.setOutput(hitID, keys[h], batch_out[h]);
          cluTrackScore += batch_out[h][0];
          cluShowerScore += batch_out[h][1];
          cluNoiseScore += batch_out[h][2];
        }

        cluTrackScore = cluTrackScore/nhits;
        cluShowerScore = cluShowerScore/nhits;
        cluNoiseScore = cluNoiseScore/nhits;

        PFPTrackScore.push_back(cluTrackScore);
        PFPShowerScore.push_back(cluShowerScore);
        PFPNoiseScore.push_back(cluNoiseScore);
        // PFPMichelScore.push_back(cluMichelScore);

      }

      // get average of score in PFPTrackScore
      if (nClusters == 0) continue;
      float avgPFPTrackScore = getAvgScore(PFPTrackScore);
      float avgPFPShowerScore = getAvgScore(PFPShowerScore);
      float avgPFPNoiseScore = getAvgScore(PFPNoiseScore);
      // float avgPFPMichelScore = getAvgScore(PFPMichelScore);

      // std::cout << "--------------------------------------" << std::endl;
      // std::cout << "PFP " << thispfp->Self() << " has average track score: " << avgPFPTrackScore << std::endl;
      // std::cout << "PFP " << thispfp->Self() << " has average shower score: " << avgPFPShowerScore << std::endl;
      // std::cout << "PFP " << thispfp->Self() << " has average noise score: " << avgPFPNoiseScore << std::endl;
      // std::cout << "PFP " << thispfp->Self() << " has average michel score: " << avgPFPMichelScore << std::endl;
      // std::cout << "--------------------------------------" << std::endl;
      // pfpScores->emplace_back(avgPFPTrackScore, avgPFPShowerScore, avgPFPNoiseScore, avgPFPMichelScore, nClusters);
      pfpScores->emplace_back(avgPFPTrackScore, avgPFPShowerScore, avgPFPNoiseScore, nClusters);
      util::CreateAssn(*this, evt, *pfpScores, thispfp, *pfpAssns);
    }
    
    fMVAWriter.saveOutputs(evt);
    evt.put(std::move(pfpScores));
    evt.put(std::move(pfpAssns));
  }
  // ------------------------------------------------------
  // get average score, drop outliers
  float CNNID::getAvgScore(std::vector<float> scores){
    int nScores = scores.size();
    int nGoodScores = 0;
    float sumScore = 0.;
    float avgScore = 0.;
    int erase_index = -1;

    // if valid score on all 3 clusters, remove outliers
    if (nScores > 2) {
      std::sort(scores.begin(), scores.end());
      float diff_low = scores[1] - scores[0];
      float diff_high = scores[scores.size()-1] - scores[scores.size()-2];
      if (diff_low > 0.5) {
        erase_index = 0;
      }
      else if (diff_high > 0.5) {
        erase_index = scores.size()-1;
      }
    }

    for (int i = 0; i < nScores; i++) {
      if (i == erase_index) {
        continue;
      }
      sumScore += scores[i];
      nGoodScores++;
    }

    avgScore = sumScore/nGoodScores;
    return avgScore;
  }

  // ------------------------------------------------------

  DEFINE_ART_MODULE(CNNID)
}
