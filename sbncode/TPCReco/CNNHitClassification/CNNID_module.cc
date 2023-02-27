/////////////////////////////////////////////////////////////////////////////////////
// Class:       CNNID
// Module Type: producer
// File:        CNNID_module.cc
//
// Apply pre-trained CNN to hits associated with Cluster to get track/shower/noise/michel for each hit. 
// Store the average score for PFPs that are not clear cosmics
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
        Comment("tag of clusters to be used as a source of EM/track / Michel tagged new clusters "
                "(incl. single-hit clusters ) using accumulated results from hits")};

      fhicl::Atom<art::InputTag> PFParticleModuleLabel{
        Name("PFParticleModuleLabel"),
        Comment("tag of PFParticle to be used as a source of EM/track / Michel tagged new clusters "
                "(incl. single-hit clusters ) using accumulated results from hits")};

      fhicl::Atom<bool> DoNeutrinoOnly{
        Name("DoNeutrinoOnly"),
        Comment("option to apply CNN to Pandora neutrino primary hits in the event")};

      fhicl::Atom<bool> DoPFPs{
        Name("DoPFPs"),
        Comment("option to apply CNN to PFPs in the event")};

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

    bool isViewSelected(int view) const;

    size_t fBatchSize;
    nnet::PointIdAlg fPointIdAlg;
    anab::MVAWriter<4> fMVAWriter; 

    art::InputTag fWireProducerLabel;
    art::InputTag fHitModuleLabel;
    art::InputTag fClusterModuleLabel;
    art::InputTag fPFParticleModuleLabel;
    bool fDoNeutrinoOnly;
    bool fDoPFPs;

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
    , fDoPFPs(config().DoPFPs())
    , fViews(config().Views())
  {
    fMVAWriter.produces_using<recob::Hit>();

    if (fDoPFPs == true){
      produces<std::vector<PFPCNNScore>>();
      produces<art::Assns<recob::PFParticle, PFPCNNScore>>();
    }
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

    // create hitmap for non-clear-cosmic hits
    CNNID::cryo_tpc_view_keymap hitMap;
    for (const art::Ptr<recob::PFParticle> &thispfp : pfpList){
        auto const &clusters = cluFMpfp.at(thispfp->Self());
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metas = metaFMpfp.at(thispfp.key());
        const art::Ptr<larpandoraobj::PFParticleMetadata> meta = metas[0];
        auto const &properties = meta->GetPropertiesMap();
        bool is_clear_cosmic = (properties.count("IsClearCosmic")); 

        if (fDoNeutrinoOnly) {
            if (is_clear_cosmic) continue;
        }

        for (auto const &cluster : clusters){
            auto const &cluHitPtrList = hitFMclu.at(cluster.key());
            for (auto const& hit : cluHitPtrList) {
              view = hit->WireID().Plane;
              if (!isViewSelected(view)) continue;
              cryo = hit->WireID().Cryostat;
              tpc = hit->WireID().TPC;

              hitMap[cryo][tpc][view].push_back(hit.key());
            }
        }
    }

    // classify hits 
    auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);
    std::vector<art::Ptr<recob::Hit>> hitPtrList;
    art::fill_ptr_vector(hitPtrList, hitListHandle);

    auto hitID = fMVAWriter.initOutputs<recob::Hit>(
      fHitModuleLabel, hitPtrList.size(), fPointIdAlg.outputLabels());

    std::vector<char> hitInFA(hitPtrList.size(),0); 

    for (auto const& mcryo : hitMap) {
        cryo = mcryo.first;
        for (auto const& mtpc : mcryo.second) {
            tpc = mtpc.first;
            for (auto const& mview : mtpc.second) {
                view = mview.first;
                if (!isViewSelected(view)) continue; 
                fPointIdAlg.setWireDriftData(clockData, detProp, *wireHandle, view, tpc, cryo);
                // apply CNN
                for (size_t idx = 0; idx < mview.second.size(); idx += fBatchSize) {
                    std::vector<std::pair<unsigned int, float>> points;
                    std::vector<size_t> keys;
                    for (size_t k = 0; k < fBatchSize; ++k) {
                        if (idx + k >= mview.second.size()) { break; } 
                        size_t h = mview.second[idx + k]; 
                        const recob::Hit& hit = *(hitPtrList[h]);
                        points.emplace_back(hit.WireID().Wire, hit.PeakTime());
                        keys.push_back(h);
                    }
                    auto batch_out = fPointIdAlg.predictIdVectors(points);
                    if (points.size() != batch_out.size()) {
                        throw cet::exception("CNNID") << "hit processing with CNN failed" << std::endl;
                    }
                    for (size_t k = 0; k < points.size(); ++k) {
                        size_t h = keys[k];
                        fMVAWriter.setOutput(hitID, h, batch_out[k]);
                        if (fPointIdAlg.isInsideFiducialRegion(points[k].first, points[k].second)) {
                            hitInFA[h] = 1;
                        }
                    }
                } 
            }
        }
    }

    if (fDoPFPs) {
      std::vector<bool> hitUsed(hitPtrList.size(), false); // tag hits used in clusters
      for (const art::Ptr<recob::PFParticle> &pfpPtr : pfpList){
          auto const &clusters = cluFMpfp.at(pfpPtr->Self());
          int nClusters = clusters.size();
          float pfpTrackScore = 0.;
          float pfpShowerScore = 0.;
          float pfpNoiseScore = 0.;
          float pfpMichelScore = 0.;

          for (auto const &cluster : clusters){
              auto v = hitFMclu.at(cluster.key());
              if (v.empty()) continue;
              for (auto const& hit : v){
                  if (hitUsed[hit.key()]) {
                      mf::LogWarning("CNNID") << "hit already used in another cluster";
                  } // check if hit is already used
                  hitUsed[hit.key()] = true;
              } // hit

              auto vout = fMVAWriter.getOutput<recob::Hit>(
                      v, [&](art::Ptr<recob::Hit> const& ptr) { return (float)hitInFA[ptr.key()]; });

              if ((vout[0]>=0) && (vout[1]>=0) && (vout[2]>=0) && (vout[3]>=0)) {
                  pfpTrackScore = pfpTrackScore+vout[0];
                  pfpShowerScore = pfpShowerScore+vout[1];
                  pfpNoiseScore = pfpNoiseScore+vout[2];
                  pfpMichelScore = pfpMichelScore+vout[3];
              }
          }

          float avgTrackScore = pfpTrackScore/nClusters;
          float avgShowerScore = pfpShowerScore/nClusters;
          float avgNoiseScore = pfpNoiseScore/nClusters;
          float avgMichelScore = pfpMichelScore/nClusters;

          if (avgTrackScore == 0.) avgTrackScore = -999.;
          if (avgShowerScore == 0.) avgShowerScore = -999.;
          if (avgNoiseScore == 0.) avgNoiseScore = -999.;
          if (avgMichelScore == 0.) avgMichelScore = -999.;

          pfpScores->emplace_back(avgTrackScore, avgShowerScore, avgNoiseScore, avgMichelScore, nClusters);
          util::CreateAssn(*this, evt, *pfpScores, pfpPtr, *pfpAssns);
        } // pfp
          
      evt.put(std::move(pfpScores));
      evt.put(std::move(pfpAssns));

    } // fDoPFPs 

    fMVAWriter.saveOutputs(evt);
  }
  // ------------------------------------------------------

  bool
  CNNID::isViewSelected(int view) const
  {
    if (fViews.empty())
      return true;
    else {
      bool selected = false;
      for (auto k : fViews)
        if (k == view) {
          selected = true;
          break;
        }
      return selected;
    }
  }
  // ------------------------------------------------------

  DEFINE_ART_MODULE(CNNID)

}
