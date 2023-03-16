/////////////////////////////////////////////////////////////////////////////////////
// Class:       CNNID
// Module Type: producer
// File:        CNNID_module.cc
//
// Apply pre-trained CNN to hits associated with Cluster to get track/shower/noise/michel for each hit. 
// Store the average score for PFPs that are not clear cosmics
// Also calculate average Michel score around the endpoint of a PFP
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
    // these types to be replaced with use of feature proposed in redmine #12602
    typedef std::unordered_map<unsigned int, std::vector<size_t>> tpc_keymap;
    typedef std::unordered_map<unsigned int, double> wire_tick_map;
    typedef std::unordered_map<size_t, wire_tick_map> clukey_wire_tick_map;
    typedef std::unordered_map<size_t, std::vector<art::Ptr<recob::Hit>>> clukey_hitkey_map;
    typedef std::unordered_map<size_t, clukey_hitkey_map> pfpkey_clukey_hitkey_map;

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

      fhicl::Atom<bool> DoNeutrinoOnly{
        Name("DoNeutrinoOnly"),
        Comment("option to apply CNN to Pandora neutrino primary hits in the event")};

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
    bool fDoNeutrinoOnly;
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
    , fDoNeutrinoOnly(config().DoNeutrinoOnly())
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

    auto pfpScores = std::make_unique<std::vector<PFPCNNScore>>();
    auto pfpAssns = std::make_unique<art::Assns<recob::PFParticle, PFPCNNScore>>();

    auto wireHandle = evt.getValidHandle<std::vector<recob::Wire>>(fWireProducerLabel);
    unsigned int cryo, tpc, view;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

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

    // hitmap for clusterkey - endpointskey
    CNNID::clukey_hitkey_map cluPointsMap;
    CNNID::clukey_hitkey_map cluEndHits;
    CNNID::clukey_wire_tick_map cluEndPointsMap;
    CNNID::clukey_hitkey_map endPointsMap;
    std::vector<bool> UsedHits(hitPtrList.size(), false); // tag hits used in clusters
    std::vector<art::Ptr<recob::Hit>> cosmicHits;
                                                          
    auto hitID = fMVAWriter.initOutputs<recob::Hit>(
      fHitModuleLabel, hitPtrList.size(), fPointIdAlg.outputLabels());

    // find endpoints of clusters, calculate Michel score around endpoints
    for (const art::Ptr<recob::PFParticle> &thispfp : pfpList){
      std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metas = metaFMpfp.at(thispfp.key());
      const art::Ptr<larpandoraobj::PFParticleMetadata> meta = metas[0];
      auto const &properties = meta->GetPropertiesMap();
      bool is_clear_cosmic = (properties.count("IsClearCosmic")); 

      auto const &clusters = cluFMpfp.at(thispfp->Self());
      for (auto const &cluster : clusters){
        if (is_clear_cosmic) {
          auto &cluHitPtrList = hitFMclu.at(cluster.key());
          for (const art::Ptr<recob::Hit> &hit : cluHitPtrList) {
            cosmicHits.push_back(hit);
          }
          continue;
        }
        int end_wire = cluster->EndWire();
        auto end_tick = cluster->EndTick();
        cluEndPointsMap[cluster.key()][end_wire] = end_tick;
      }
    }

    auto nonCosmicHits = hitPtrList;
//    std::cout << "number of cosmic hits " << cosmicHits.size() << std::endl;
    nonCosmicHits.erase(std::remove_if(nonCosmicHits.begin(), nonCosmicHits.end(), [&cosmicHits](const art::Ptr<recob::Hit> &hit) 
          { return std::find(cosmicHits.begin(), cosmicHits.end(), hit) != cosmicHits.end(); }), nonCosmicHits.end());

    for (auto const& hit : nonCosmicHits){
      int wireid = hit->WireID().Wire;
      auto peakt = hit->PeakTime();
      auto view = hit->WireID().Plane;
      for (auto const& clukey: cluEndPointsMap) {
        auto clu = clukey.first;
        auto cluster = *(cluPtrList[clu]);
        auto cluview = cluster.Plane().Plane;
        for (auto const& wire: clukey.second) {
          int end_wire = wire.first;
          auto end_tick = wire.second;
          auto cluhits = hitFMclu.at(clu);
          if ((view==cluview)&&(std::abs(wireid-end_wire)<30)&&(std::abs(peakt-end_tick)<30)){
            auto it = std::find(cluhits.begin(), cluhits.end(), hit);
            if (it != cluhits.end()) continue;
            cluEndHits[clu].push_back(hit);
          }
        }
      }
    }

    // loop through PFPs, save to CNNScore
    for (const art::Ptr<recob::PFParticle> &thispfp : pfpList){
      std::vector<float> PFPTrackScore; 
      std::vector<float> PFPShowerScore; 
      std::vector<float> PFPNoiseScore;
      std::vector<float> PFPMichelScore;

      std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metas = metaFMpfp.at(thispfp.key());
      const art::Ptr<larpandoraobj::PFParticleMetadata> meta = metas[0];
      auto const &properties = meta->GetPropertiesMap();
      bool is_clear_cosmic = (properties.count("IsClearCosmic")); 
      if (is_clear_cosmic) continue;

      auto const &clusters = cluFMpfp.at(thispfp->Self());

      std::cout << "**********************************************************" << std::endl;
      int nClusters = 0; 
      for (auto const &cluster : clusters){
        auto &cluHitPtrList = hitFMclu.at(cluster.key());
        int nHitsClu = cluHitPtrList.size();
        if (nHitsClu == 0) continue;
        nClusters++;

        cryo = cluster->Plane().Cryostat;
        view = cluster->Plane().Plane;
        tpc = cluster->Plane().TPC;
//        std::cout << "cluster on cyro " << cryo << " view " << view << " tpc " << tpc << " has " << nHitsClu << " hits " <<std::endl;
        fPointIdAlg.setWireDriftData(clockData, detProp, *wireHandle, view, tpc, cryo);

        float cluTrackScore = 0.;
        float cluShowerScore = 0.;
        float cluNoiseScore = 0.;
        float cluMichelScore = 0.;

        std::vector<art::Ptr<recob::Hit>> cluSparseHits; // sparse vector of hits in this cluster
        for (int i = 0; i < nHitsClu; i++){
          if ((nHitsClu > fSparseLengthCut) && (i % fSparseRate != 0)) continue;
          cluSparseHits.push_back(cluHitPtrList[i]);
        }

        // apply CNN to cluster hits
        std::vector<std::pair<unsigned int, float>> points;
        std::vector<size_t> keys;
        for (auto const& hit : cluSparseHits){
          points.emplace_back(hit->WireID().Wire, hit->PeakTime());
          keys.push_back(hit.key());
        }

//        std::cout << "number of points to apply CNN " << points.size() << std::endl;
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

        if (isnan(cluTrackScore)) cluTrackScore = -1.;
        if (isnan(cluShowerScore)) cluShowerScore = -1.;
        if (isnan(cluNoiseScore)) cluNoiseScore = -1.;

        PFPTrackScore.push_back(cluTrackScore);
        PFPShowerScore.push_back(cluShowerScore);
        PFPNoiseScore.push_back(cluNoiseScore);

        // apply CNN to cluster end hits
        auto endHits = cluEndHits[cluster.key()];
        std::vector<std::pair<unsigned int, float>> michelPoints;
        std::vector<size_t> michelKeys;
        for (auto const& hit : endHits){
          michelPoints.emplace_back(hit->WireID().Wire, hit->PeakTime());
          michelKeys.push_back(hit.key());
        }

        size_t nMichelhits = michelPoints.size();
//        std::cout << "number of michel points to apply CNN " << nMichelhits << std::endl;
        auto batch_out_michel = fPointIdAlg.predictIdVectors(michelPoints);
        if (michelPoints.size() != batch_out_michel.size()) {
          throw cet::exception("CNNID") << "Size mismatch between input and output vectors";
        }
        for (size_t h = 0; h < nMichelhits; h++) {
          fMVAWriter.setOutput(hitID, michelKeys[h], batch_out_michel[h]);
          cluMichelScore += batch_out_michel[h][3];
        }
        cluMichelScore = cluMichelScore/nMichelhits;
        if (isnan(cluMichelScore)) cluMichelScore = -1.;
        PFPMichelScore.push_back(cluMichelScore);

      }

      std::cout << "PFP " << thispfp->Self() << " has " << nClusters << " clusters" << std::endl;

      // get average of score in PFPTrackScore
      float avgPFPTrackScore = getAvgScore(PFPTrackScore);
      float avgPFPShowerScore = getAvgScore(PFPShowerScore);
      float avgPFPNoiseScore = getAvgScore(PFPNoiseScore);
      float avgPFPMichelScore = getAvgScore(PFPMichelScore);

      std::cout << "--------------------------------------" << std::endl;
      std::cout << "PFP " << thispfp->Self() << " has average track score: " << avgPFPTrackScore << std::endl;
      std::cout << "PFP " << thispfp->Self() << " has average shower score: " << avgPFPShowerScore << std::endl;
      std::cout << "PFP " << thispfp->Self() << " has average noise score: " << avgPFPNoiseScore << std::endl;
      std::cout << "PFP " << thispfp->Self() << " has average michel score: " << avgPFPMichelScore << std::endl;
      std::cout << "--------------------------------------" << std::endl;
      pfpScores->emplace_back(avgPFPTrackScore, avgPFPShowerScore, avgPFPNoiseScore, avgPFPMichelScore, nClusters);
      util::CreateAssn(*this, evt, *pfpScores, thispfp, *pfpAssns);
    }
    
    fMVAWriter.saveOutputs(evt);
    evt.put(std::move(pfpScores));
    evt.put(std::move(pfpAssns));
  }
  // ------------------------------------------------------
  // get average score, drop outliers
  float CNNID::getAvgScore(std::vector<float> scores){
    float avgScore = 0.;
    float sumScore = 0.;
    int nScores = scores.size();
    int nGoodScores = 0;
    int erase_index = -1;

    // if has score on all 3 clusters, remove outliers
    //sort scores
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
        std::cout << "skipping outlier score " << scores[i] << std::endl;
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
