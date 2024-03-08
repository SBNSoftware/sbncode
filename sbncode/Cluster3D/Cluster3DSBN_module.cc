/**
 *  @file   Cluster3DSBN_module.cc
 *
 *          Class:       Cluster3DSBN
 *          Module Type: Producer
 *
 *  @brief  Producer module to create 3D clusters from input recob::Hit objects
 *
 *          This producer module will drive the 3D association of recob::Hit objects
 *          to form 3D clusters. This information will be output as:
 *          1) a PFParticle to anchor all the other objects (as associations)
 *          2) three recob::Cluster objects representing 2D hit clusterswith associations
 *             to the 2D hits comprising each of them
 *          3) One or more recob::PCAxis objects representing the Principal Components
 *             Analysis output of the space points associated to the 3D objects
 *          4) recob::SpacePoints representing the accepted 3D points for each PFParticle
 *          5) recob::Seed objects and associated seed hits representing candidate straight
 *             line segments in the space point collection.
 *
 *          The module has two main sections
 *          1) Find the 3D clusters of 3D hits
 *          2) Get the output objects for each of the 3D clusters
 *          See the code below for more detail on these steps
 *
 *          Note that the general 3D cluster suite of algorithms make extensive use of a set of data objects
 *          which contain volatile data members. At the end of the routine these are used to make the output
 *          LArSoft data products described above. See LarData/RecoObjects/Cluster3D.h
 *
 *          Configuration parameters:
 *          HitFinderModuleLabel:         the producer module responsible for making the recob:Hits to use
 *          EnableMonitoring:             if true then basic monitoring of the module performed
 *          ClusterAlg:                   Parameter block required by the 3D clustering algorithm
 *          PrincipalComponentsAlg:       Parameter block required by the Principal Components Analysis Algorithm
 *          SkeletonAlg:                  Parameter block required by the 3D skeletonization algorithm
 *          SeedFinderAlg:                Parameter block required by the Hough Seed Finder algorithm
 *          PCASeedFinderAlg:             Parameter block required by the PCA Seed Finder algorithm
 *          ParrallelHitsAlg:             Parameter block required by the parallel hits algorithm
 *
 *          The current producer module does not try to analyze or break apart PFParticles
 *          so, for example, all tracks emanating from a common vertex will be associated
 *          to a single PFParticle
 *
 *  @author usher@slac.stanford.edu
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Edge.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "larreco/RecoAlg/Cluster3DAlgs/HoughSeedFinderAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterModAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IClusterParamsBuilder.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IHit3DBuilder.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PCASeedFinderAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/ParallelHitsSeedFinderAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/PrincipalComponentsAlg.h"
#include "larreco/RecoAlg/Cluster3DAlgs/SkeletonAlg.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "larreco/RecoAlg/ClusterRecoUtil/OverriddenClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"

// ROOT includes
#include "TTree.h"
#include "TVector3.h"

// std includes
#include <iostream>
#include <memory>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d {
  using Hit3DToSPPtrMap =
    std::unordered_map<const reco::ClusterHit3D*, art::Ptr<recob::SpacePoint>>;
  using RecobHitVector = art::PtrVector<recob::Hit>;
  using RecobSpacePointVector = art::PtrVector<recob::SpacePoint>;

  //------------------------------------------------------------------------------------------------------------------------------------------
  // Definition of the producer module here

  /**
 *  @brief  Definition of the Cluster3D class
 */
  class Cluster3DSBN : public art::EDProducer {
  public:
    explicit Cluster3DSBN(fhicl::ParameterSet const& pset);

  private:
    void beginJob() override;
    void produce(art::Event& evt) override;

    class ArtOutputHandler {
    public:
      ArtOutputHandler(art::Event& evt,
                       std::string& pathName,
                       std::string& vertexName,
                       std::string& extremeName)
        : artSpacePointVector(std::make_unique<std::vector<recob::SpacePoint>>())
        , artSPHitAssociations(std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>())
        , artPCAxisVector(std::make_unique<std::vector<recob::PCAxis>>())
        , artPFParticleVector(std::make_unique<std::vector<recob::PFParticle>>())
        , artClusterVector(std::make_unique<std::vector<recob::Cluster>>())
        , artPathPointVector(std::make_unique<std::vector<recob::SpacePoint>>())
        , artVertexPointVector(std::make_unique<std::vector<recob::SpacePoint>>())
        , artExtremePointVector(std::make_unique<std::vector<recob::SpacePoint>>())
        , artSeedVector(std::make_unique<std::vector<recob::Seed>>())
        , artEdgeVector(std::make_unique<std::vector<recob::Edge>>())
        , artPathEdgeVector(std::make_unique<std::vector<recob::Edge>>())
        , artVertexEdgeVector(std::make_unique<std::vector<recob::Edge>>())
        , artClusterAssociations(std::make_unique<art::Assns<recob::Cluster, recob::Hit>>())
        , artPFPartAxisAssociations(std::make_unique<art::Assns<recob::PFParticle, recob::PCAxis>>())
        , artPFPartClusAssociations(std::make_unique<art::Assns<recob::PFParticle, recob::Cluster>>())
        , artPFPartSPAssociations(std::make_unique<art::Assns<recob::SpacePoint, recob::PFParticle>>())
        , artPFPartPPAssociations(std::make_unique<art::Assns<recob::SpacePoint, recob::PFParticle>>())
        , artPFPartSeedAssociations(std::make_unique<art::Assns<recob::PFParticle, recob::Seed>>())
        , artPFPartEdgeAssociations(std::make_unique<art::Assns<recob::Edge, recob::PFParticle>>())
        , artPFPartPathEdgeAssociations(std::make_unique<art::Assns<recob::Edge, recob::PFParticle>>())
        , artSeedHitAssociations(std::make_unique<art::Assns<recob::Seed, recob::Hit>>())
        , artPPHitAssociations(std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>())
        , artEdgeSPAssociations(std::make_unique<art::Assns<recob::SpacePoint, recob::Edge>>())
        , artEdgePPAssociations(std::make_unique<art::Assns<recob::SpacePoint, recob::Edge>>())
        , fEvt(evt)
        , fSPPtrMaker(evt)
        , fSPPtrMakerPath(evt, pathName)
        , fEdgePtrMaker(evt)
        , fEdgePtrMakerPath(evt, pathName)
        , fPathName(pathName)
        , fVertexName(vertexName)
        , fExtremeName(extremeName)
      {}

      void
      makeClusterHitAssns(RecobHitVector& recobHits)
      {
        util::CreateAssn(fEvt, *artClusterVector, recobHits, *artClusterAssociations);
      }

      void
      makeSpacePointHitAssns(std::vector<recob::SpacePoint>& spacePointVector,
                             RecobHitVector& recobHits,
                             art::Assns<recob::Hit, recob::SpacePoint>& spHitAssns,
                             const std::string& path = "")
      {
        for (auto& hit : recobHits)
          util::CreateAssn(fEvt, spacePointVector, hit, spHitAssns, path);
      }

      void
      makePFPartPCAAssns()
      {
        util::CreateAssn(fEvt,
                         *artPFParticleVector,
                         *artPCAxisVector,
                         *artPFPartAxisAssociations,
                         artPCAxisVector->size() - 2,
                         artPCAxisVector->size());
      }

      void
      makePFPartSeedAssns(size_t numSeedsStart)
      {
        util::CreateAssn(fEvt,
                         *artPFParticleVector,
                         *artSeedVector,
                         *artPFPartSeedAssociations,
                         numSeedsStart,
                         artSeedVector->size());
      }

      void
      makePFPartClusterAssns(size_t clusterStart)
      {
        util::CreateAssn(fEvt,
                         *artPFParticleVector,
                         *artClusterVector,
                         *artPFPartClusAssociations,
                         clusterStart,
                         artClusterVector->size());
      }

      void
      makePFPartSpacePointAssns(
        std::vector<recob::SpacePoint>& spacePointVector,
        art::Assns<recob::SpacePoint, recob::PFParticle>& pfPartSPAssociations,
        size_t spacePointStart,
        const std::string& instance = "")
      {
        for (size_t idx = spacePointStart; idx < spacePointVector.size(); idx++) {
          art::Ptr<recob::SpacePoint> spacePoint = makeSpacePointPtr(idx, instance);
          util::CreateAssn(fEvt, *artPFParticleVector, spacePoint, pfPartSPAssociations);
        }
      }

      void
      makePFPartEdgeAssns(std::vector<recob::Edge>& edgeVector,
                          art::Assns<recob::Edge, recob::PFParticle>& pfPartEdgeAssociations,
                          size_t edgeStart,
                          const std::string& instance = "")
      {
        for (size_t idx = edgeStart; idx < edgeVector.size(); idx++) {
          art::Ptr<recob::Edge> edge = makeEdgePtr(idx, instance);

          util::CreateAssn(fEvt, *artPFParticleVector, edge, pfPartEdgeAssociations);
        }
      }

      void
      makeEdgeSpacePointAssns(std::vector<recob::Edge>& edgeVector,
                              RecobSpacePointVector& spacePointVector,
                              art::Assns<recob::SpacePoint, recob::Edge>& edgeSPAssociations,
                              const std::string& path = "")
      {
        for (auto& spacePoint : spacePointVector)
          util::CreateAssn(fEvt, edgeVector, spacePoint, edgeSPAssociations, path);
      }

      void
      outputObjects(bool spacePointsOnly)
      {
        fEvt.put(std::move(artSpacePointVector));
        fEvt.put(std::move(artSPHitAssociations));

//        if (!spacePointsOnly)
        {
          fEvt.put(std::move(artPCAxisVector));
          fEvt.put(std::move(artPFParticleVector));
          fEvt.put(std::move(artClusterVector));
          fEvt.put(std::move(artSeedVector));
          fEvt.put(std::move(artEdgeVector));
          fEvt.put(std::move(artPFPartAxisAssociations));
          fEvt.put(std::move(artPFPartClusAssociations));
          fEvt.put(std::move(artClusterAssociations));
          fEvt.put(std::move(artPFPartSPAssociations));
          fEvt.put(std::move(artPFPartSeedAssociations));
          fEvt.put(std::move(artPFPartEdgeAssociations));
          fEvt.put(std::move(artSeedHitAssociations));
          fEvt.put(std::move(artEdgeSPAssociations));
          fEvt.put(std::move(artPathEdgeVector), fPathName);
          fEvt.put(std::move(artPathPointVector), fPathName);
          fEvt.put(std::move(artEdgePPAssociations), fPathName);
          fEvt.put(std::move(artPFPartPPAssociations), fPathName);
          fEvt.put(std::move(artPFPartPathEdgeAssociations), fPathName);
          fEvt.put(std::move(artPPHitAssociations), fPathName);
          fEvt.put(std::move(artVertexEdgeVector), fVertexName);
          fEvt.put(std::move(artVertexPointVector), fVertexName);
          fEvt.put(std::move(artExtremePointVector), fExtremeName);
        }
      }

      art::Ptr<recob::SpacePoint>
      makeSpacePointPtr(size_t index, const std::string& instance = "")
      {
        if (empty(instance)) { return fSPPtrMaker(index); }
        return fSPPtrMakerPath(index);
      };

      art::Ptr<recob::Edge>
      makeEdgePtr(size_t index, const std::string& instance = "")
      {
        if (empty(instance)) { return fEdgePtrMaker(index); }
        return fEdgePtrMakerPath(index);
      };

      std::unique_ptr<std::vector<recob::SpacePoint>> artSpacePointVector;
      std::unique_ptr<art::Assns<recob::Hit, recob::SpacePoint>> artSPHitAssociations;

      std::unique_ptr<std::vector<recob::PCAxis>> artPCAxisVector;
      std::unique_ptr<std::vector<recob::PFParticle>> artPFParticleVector;
      std::unique_ptr<std::vector<recob::Cluster>> artClusterVector;
      std::unique_ptr<std::vector<recob::SpacePoint>> artPathPointVector;
      std::unique_ptr<std::vector<recob::SpacePoint>> artVertexPointVector;
      std::unique_ptr<std::vector<recob::SpacePoint>> artExtremePointVector;
      std::unique_ptr<std::vector<recob::Seed>> artSeedVector;
      std::unique_ptr<std::vector<recob::Edge>> artEdgeVector;
      std::unique_ptr<std::vector<recob::Edge>> artPathEdgeVector;
      std::unique_ptr<std::vector<recob::Edge>> artVertexEdgeVector;

      std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>> artClusterAssociations;
      std::unique_ptr<art::Assns<recob::PFParticle, recob::PCAxis>> artPFPartAxisAssociations;
      std::unique_ptr<art::Assns<recob::PFParticle, recob::Cluster>> artPFPartClusAssociations;
      std::unique_ptr<art::Assns<recob::SpacePoint, recob::PFParticle>> artPFPartSPAssociations;
      std::unique_ptr<art::Assns<recob::SpacePoint, recob::PFParticle>> artPFPartPPAssociations;
      std::unique_ptr<art::Assns<recob::PFParticle, recob::Seed>> artPFPartSeedAssociations;
      std::unique_ptr<art::Assns<recob::Edge, recob::PFParticle>> artPFPartEdgeAssociations;
      std::unique_ptr<art::Assns<recob::Edge, recob::PFParticle>> artPFPartPathEdgeAssociations;
      std::unique_ptr<art::Assns<recob::Seed, recob::Hit>> artSeedHitAssociations;
      std::unique_ptr<art::Assns<recob::Hit, recob::SpacePoint>> artPPHitAssociations;
      std::unique_ptr<art::Assns<recob::SpacePoint, recob::Edge>> artEdgeSPAssociations;
      std::unique_ptr<art::Assns<recob::SpacePoint, recob::Edge>> artEdgePPAssociations;

    private:
      art::Event& fEvt;
      art::PtrMaker<recob::SpacePoint> fSPPtrMaker;
      art::PtrMaker<recob::SpacePoint> fSPPtrMakerPath;
      art::PtrMaker<recob::Edge> fEdgePtrMaker;
      art::PtrMaker<recob::Edge> fEdgePtrMakerPath;
      std::string& fPathName;
      std::string& fVertexName;
      std::string& fExtremeName;
    };

    /**
     *  @brief  Event Preparation
     *
     *  @param  evt  the ART event
     */
    void PrepareEvent(const art::Event& evt);

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();

    /**
     *  @brief An interface to the seed finding algorithm
     *
     *  @param evt          the ART event
     *  @param cluster      structure of information representing a single cluster
     *  @param hitToPtrMap  This maps our Cluster2D hits back to art Ptr's to reco Hits
     *  @param seedVec      the output vector of candidate seeds
     *  @param seedHitAssns the associations between the seeds and the 2D hits making them
     */
    void findTrackSeeds(art::Event& evt,
                        reco::ClusterParameters& cluster,
                        IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap,
                        std::vector<recob::Seed>& seedVec,
                        art::Assns<recob::Seed, recob::Hit>& seedHitAssns) const;

    /**
     *  @brief Attempt to split clusters by using a minimum spanning tree
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    void splitClustersWithMST(reco::ClusterParameters& clusterParameters,
                              reco::ClusterParametersList& clusterParametersList) const;

    /**
     *  @brief Attempt to split clusters using the output of the Hough Filter
     *
     *  @param clusterParameters     The given cluster parameters object to try to split
     *  @param clusterParametersList The list of clusters
     */
    void splitClustersWithHough(reco::ClusterParameters& clusterParameters,
                                reco::ClusterParametersList& clusterParametersList) const;

    /**
     *  @brief Special routine to handle creating and saving space points
     *
     *  @param output                the object containting the art output
     *  @param clusterParameters     Cluster info to output (in internal format)
     *  @param pfParticleParent      The parent ID reference for the output PFParticle
     *  @param hitToPtrMap           This maps our Cluster2D hits back to art Ptr's to reco Hits
     */
    void MakeAndSaveSpacePoints(ArtOutputHandler& output,
                                std::vector<recob::SpacePoint>& spacePointVec,
                                art::Assns<recob::Hit, recob::SpacePoint>& spHitAssns,
                                reco::HitPairListPtr& clusHitPairVector,
                                IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap,
                                Hit3DToSPPtrMap& hit3DToSPPtrMap,
                                const std::string& path = "") const;

    /**
     *  @brief Special routine to handle creating and saving space points
     *
     *  @param output                the object containting the art output
     *  @param clusHitPairVector     List of 3D hits to output as "extreme" space points
     */
    void MakeAndSaveKinkPoints(ArtOutputHandler& output,
                               reco::ConvexHullKinkTupleList& clusHitPairVector) const;

    /**
     *  @brief Special routine to handle creating and saving space points & edges associated to voronoi diagrams
     *
     *  @param output                the object containting the art output
     *  @param vertexList            list of vertices in the diagram
     *  @param HalfEdgeList          list of half edges in the diagram
     */
    void MakeAndSaveVertexPoints(ArtOutputHandler&,
                                 dcel2d::VertexList&,
                                 dcel2d::HalfEdgeList&) const;

    /**
     *  @brief Special routine to handle creating and saving space points & edges PCA points
     *
     *  @param output                the object containting the art output
     *  @param clusterParamsList     List of clusters to get PCA's from
     */
    using IdxToPCAMap = std::map<size_t, const reco::PrincipalComponents*>;

    void MakeAndSavePCAPoints(ArtOutputHandler&,
                              const reco::PrincipalComponents&,
                              IdxToPCAMap&) const;

    /**
     *  @brief This will produce art output for daughters noting that it needs to be done recursively
     *
     *  @param output                the object containting the art output
     *  @param clusterParameters     Cluster info to output (in internal format)
     *  @param pfParticleParent      The parent ID reference for the output PFParticle
     *  @param daughterList          List of PFParticle indices for stored daughters
     *  @param hitToPtrMap           This maps our Cluster2D hits back to art Ptr's to reco Hits
     */
    size_t FindAndStoreDaughters(util::GeometryUtilities const& gser,
                                 ArtOutputHandler& output,
                                 reco::ClusterParameters& clusterParameters,
                                 size_t pfParticleParent,
                                 IdxToPCAMap& idxToPCAMap,
                                 IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap,
                                 Hit3DToSPPtrMap& hit3DToSPPtrMap,
                                 Hit3DToSPPtrMap& best3DToSPPtrMap) const;

    /**
     *  @brief Top level output routine, allows checking cluster status
     *
     *  @param hitPairList           List of all 3D Hits in internal Cluster3D format
     *  @param clusterParametersList Data structure containing the cluster information to output
     *  @param  hitToPtrMap          This maps our Cluster2D hits back to art Ptr's to reco Hits
     */
    void ProduceArtClusters(util::GeometryUtilities const& gser,
                            ArtOutputHandler& output,
                            reco::HitPairList& hitPairList,
                            reco::ClusterParametersList& clusterParametersList,
                            IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap) const;

    /**
     *  @brief Produces the art output from all the work done in this producer module
     *
     *  @param output                the object containting the art output
     *  @param clusterParameters     Cluster info to output (in internal format)
     *  @param pfParticleParent      The parent ID reference for the output PFParticle
     *  @param hitToPtrMap           This maps our Cluster2D hits back to art Ptr's to reco Hits
     */
    size_t ConvertToArtOutput(util::GeometryUtilities const& gser,
                              ArtOutputHandler& output,
                              reco::ClusterParameters& clusterParameters,
                              size_t pfParticleParent,
                              IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap,
                              Hit3DToSPPtrMap& hit3DToSPPtrMap,
                              Hit3DToSPPtrMap& best3DToSPPtrMap) const;

    /**
     *  @brief There are several places we will want to know if a candidate cluster is a
     *         "parallel hits" type of cluster so encapsulate that here.
     *
     *  @param pca  The Principal Components Analysis parameters for the cluster
     */
    bool
    aParallelHitsCluster(const reco::PrincipalComponents& pca) const
    {
      return fabs(pca.getEigenVectors().row(2)(0)) > m_parallelHitsCosAng &&
             3. * sqrt(pca.getEigenValues()(1)) > m_parallelHitsTransWid;
    }

    /**
     *  @brief Count number of end of line daughters
     *
     *  @param clusterParams input cluster parameters to look at
     */
    size_t countUltimateDaughters(reco::ClusterParameters& clusterParameters) const;

    /**
     *   Algorithm parameters
     */
    bool m_onlyMakSpacePoints;    ///< If true we don't do the full cluster 3D processing
    bool m_enableMonitoring;      ///< Turn on monitoring of this algorithm
    float m_parallelHitsCosAng;   ///< Cut for PCA 3rd axis angle to X axis
    float m_parallelHitsTransWid; ///< Cut on transverse width of cluster (PCA 2nd eigenvalue)

    /**
     *   Tree variables for output
     */
    TTree* m_pRecoTree;            ///<
    int m_run;                     ///<
    int m_event;                   ///<
    int m_hits;                    ///< Keeps track of the number of hits seen
    int m_hits3D;                  ///< Keeps track of the number of 3D hits made
    float m_totalTime;             ///< Keeps track of total execution time
    float m_artHitsTime;           ///< Keeps track of time to recover hits
    float m_makeHitsTime;          ///< Keeps track of time to build 3D hits
    float m_buildNeighborhoodTime; ///< Keeps track of time to build epsilon neighborhood
    float m_dbscanTime;            ///< Keeps track of time to run DBScan
    float m_clusterMergeTime;      ///< Keeps track of the time to merge clusters
    float m_pathFindingTime;       ///< Keeps track of the path finding time
    float m_finishTime;            ///< Keeps track of time to run output module
    std::string m_pathInstance;    ///< Special instance for path points
    std::string m_vertexInstance;  ///< Special instance name for vertex points
    std::string m_extremeInstance; ///< Instance name for the extreme points

    // Algorithms
    std::unique_ptr<lar_cluster3d::IHit3DBuilder>
      m_hit3DBuilderAlg; ///<  Builds the 3D hits to operate on
    std::unique_ptr<lar_cluster3d::IClusterAlg>
      m_clusterAlg; ///<  Algorithm to do 3D space point clustering
    std::unique_ptr<lar_cluster3d::IClusterModAlg>
      m_clusterMergeAlg; ///<  Algorithm to do cluster merging
    std::unique_ptr<lar_cluster3d::IClusterModAlg>
      m_clusterPathAlg; ///<  Algorithm to do cluster path finding
    std::unique_ptr<lar_cluster3d::IClusterParametersBuilder>
      m_clusterBuilder;                          ///<  Common cluster builder tool
    PrincipalComponentsAlg m_pcaAlg;             ///<  Principal Components algorithm
    SkeletonAlg m_skeletonAlg;                   ///<  Skeleton point finder
    HoughSeedFinderAlg m_seedFinderAlg;          ///<  Seed finder
    PCASeedFinderAlg m_pcaSeedFinderAlg;         ///<  Use PCA axis to find seeds
    ParallelHitsSeedFinderAlg m_parallelHitsAlg; ///<  Deal with parallel hits clusters
  };

  DEFINE_ART_MODULE(Cluster3DSBN)

} // namespace lar_cluster3d

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

  Cluster3DSBN::Cluster3DSBN(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , m_pcaAlg(pset.get<fhicl::ParameterSet>("PrincipalComponentsAlg"))
    , m_skeletonAlg(pset.get<fhicl::ParameterSet>("SkeletonAlg"))
    , m_seedFinderAlg(pset.get<fhicl::ParameterSet>("SeedFinderAlg"))
    , m_pcaSeedFinderAlg(pset.get<fhicl::ParameterSet>("PCASeedFinderAlg"))
    , m_parallelHitsAlg(pset.get<fhicl::ParameterSet>("ParallelHitsAlg"))
  {
    m_onlyMakSpacePoints = pset.get<bool>("MakeSpacePointsOnly", false);
    m_enableMonitoring = pset.get<bool>("EnableMonitoring", false);
    m_parallelHitsCosAng = pset.get<float>("ParallelHitsCosAng", 0.999);
    m_parallelHitsTransWid = pset.get<float>("ParallelHitsTransWid", 25.0);
    m_pathInstance = pset.get<std::string>("PathPointsName", "Path");
    m_vertexInstance = pset.get<std::string>("VertexPointsName", "Vertex");
    m_extremeInstance = pset.get<std::string>("ExtremePointsName", "Extreme");

    m_hit3DBuilderAlg = art::make_tool<lar_cluster3d::IHit3DBuilder>(
      pset.get<fhicl::ParameterSet>("Hit3DBuilderAlg"));
    m_clusterAlg =
      art::make_tool<lar_cluster3d::IClusterAlg>(pset.get<fhicl::ParameterSet>("ClusterAlg"));
    m_clusterMergeAlg = art::make_tool<lar_cluster3d::IClusterModAlg>(
      pset.get<fhicl::ParameterSet>("ClusterMergeAlg"));
    m_clusterPathAlg = art::make_tool<lar_cluster3d::IClusterModAlg>(
      pset.get<fhicl::ParameterSet>("ClusterPathAlg"));
    m_clusterBuilder = art::make_tool<lar_cluster3d::IClusterParametersBuilder>(
      pset.get<fhicl::ParameterSet>("ClusterParamsBuilder"));

    // Handle special case of Space Point building outputting a new hit collection
    m_hit3DBuilderAlg->produces(producesCollector());

    // we always output space points
    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::Hit, recob::SpacePoint>>();

    // Do we output anything else?
//    if (!m_onlyMakSpacePoints)
    {
      produces<std::vector<recob::PCAxis>>();
      produces<std::vector<recob::PFParticle>>();
      produces<std::vector<recob::Cluster>>();
      produces<std::vector<recob::Seed>>();
      produces<std::vector<recob::Edge>>();

      produces<art::Assns<recob::PFParticle, recob::PCAxis>>();
      produces<art::Assns<recob::PFParticle, recob::Cluster>>();
      produces<art::Assns<recob::PFParticle, recob::Seed>>();
      produces<art::Assns<recob::Edge, recob::PFParticle>>();
      produces<art::Assns<recob::Seed, recob::Hit>>();
      produces<art::Assns<recob::Cluster, recob::Hit>>();

      produces<art::Assns<recob::SpacePoint, recob::PFParticle>>();
      produces<art::Assns<recob::SpacePoint, recob::Edge>>();

      produces<std::vector<recob::SpacePoint>>(m_pathInstance);
      produces<std::vector<recob::Edge>>(m_pathInstance);
      produces<art::Assns<recob::SpacePoint, recob::PFParticle>>(m_pathInstance);
      produces<art::Assns<recob::Edge, recob::PFParticle>>(m_pathInstance);
      produces<art::Assns<recob::Hit, recob::SpacePoint>>(m_pathInstance);
      produces<art::Assns<recob::SpacePoint, recob::Edge>>(m_pathInstance);

      produces<std::vector<recob::Edge>>(m_vertexInstance);
      produces<std::vector<recob::SpacePoint>>(m_vertexInstance);

      produces<std::vector<recob::SpacePoint>>(m_extremeInstance);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  Cluster3DSBN::beginJob()
  {
    /**
     *  @brief beginJob will be tasked with initializing monitoring, in necessary, but also to init the
     *         geometry and detector services (and this probably needs to go in a "beginEvent" method?)
     */
    if (m_enableMonitoring) this->InitializeMonitoring();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  Cluster3DSBN::produce(art::Event& evt)
  {
    mf::LogInfo("Cluster3DSBN") << " *** Cluster3D::produce(...)  [Run=" << evt.run()
                                << ", Event=" << evt.id().event() << "] Starting Now! *** "
                                << std::endl;

    // Set up for monitoring the timing... at some point this should be removed in favor of
    // external profilers
    cet::cpu_timer theClockTotal;
    cet::cpu_timer theClockFinish;

    if (m_enableMonitoring) theClockTotal.start();

    // This really only does anything if we are monitoring since it clears our tree variables
    this->PrepareEvent(evt);

    // Get instances of the primary data structures needed
    reco::ClusterParametersList clusterParametersList;
    IHit3DBuilder::RecobHitToPtrMap clusterHitToArtPtrMap;
    std::unique_ptr<reco::HitPairList> hitPairList(new reco::HitPairList); // Potentially lots of hits, use heap instead of stack

    // Call the algorithm that builds 3D hits and stores the hit collection
    m_hit3DBuilderAlg->Hit3DBuilder(evt, *hitPairList, clusterHitToArtPtrMap);

    // Only do the rest if we are not in the mode of only building space points (requested by ML folks)
    if (!m_onlyMakSpacePoints) 
    {
      // Call the main workhorse algorithm for building the local version of candidate 3D clusters
      m_clusterAlg->Cluster3DHits(*hitPairList, clusterParametersList);

      // Try merging clusters
      m_clusterMergeAlg->ModifyClusters(clusterParametersList);

      // Run the path finding
      m_clusterPathAlg->ModifyClusters(clusterParametersList);
    }

    if (m_enableMonitoring) theClockFinish.start();

    // Get the art ouput object
    ArtOutputHandler output(evt, m_pathInstance, m_vertexInstance, m_extremeInstance);

    // Call the module that does the end processing (of which there is quite a bit of work!)
    // This goes here to insure that something is always written to the data store
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    util::GeometryUtilities const gser{*lar::providerFrom<geo::Geometry>(), clockData, detProp};

    ProduceArtClusters(gser, output, *hitPairList, clusterParametersList, clusterHitToArtPtrMap);

    // Output to art
    output.outputObjects(m_onlyMakSpacePoints);

    // If monitoring then deal with the fallout
    if (m_enableMonitoring) {
      theClockFinish.stop();
      theClockTotal.stop();

      m_run = evt.run();
      m_event = evt.id().event();
      m_totalTime = theClockTotal.accumulated_real_time();
      m_artHitsTime = m_hit3DBuilderAlg->getTimeToExecute(IHit3DBuilder::COLLECTARTHITS);
      m_makeHitsTime = m_hit3DBuilderAlg->getTimeToExecute(IHit3DBuilder::BUILDTHREEDHITS);

      if (!m_onlyMakSpacePoints)
      {
        m_buildNeighborhoodTime = m_clusterAlg->getTimeToExecute(IClusterAlg::BUILDHITTOHITMAP);
        m_dbscanTime = m_clusterAlg->getTimeToExecute(IClusterAlg::RUNDBSCAN) +
                       m_clusterAlg->getTimeToExecute(IClusterAlg::BUILDCLUSTERINFO);
        m_clusterMergeTime = m_clusterMergeAlg->getTimeToExecute();
        m_pathFindingTime = m_clusterPathAlg->getTimeToExecute();
      }

      m_finishTime = theClockFinish.accumulated_real_time();
      m_hits = static_cast<int>(clusterHitToArtPtrMap.size());
      m_hits3D = static_cast<int>(hitPairList->size());

      m_pRecoTree->Fill();

      mf::LogDebug("Cluster3DSBN") << "*** Cluster3DSBN total time: " << m_totalTime
                                   << ", art: " << m_artHitsTime << ", make: " << m_makeHitsTime
                                   << ", build: " << m_buildNeighborhoodTime
                                   << ", clustering: " << m_dbscanTime
                                   << ", merge: " << m_clusterMergeTime
                                   << ", path: " << m_pathFindingTime << ", finish: " << m_finishTime
                                   << std::endl;
    }

    // Will we ever get here? ;-)
    return;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  Cluster3DSBN::InitializeMonitoring()
  {
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("hits", &m_hits, "hits/I");
    m_pRecoTree->Branch("hits3D", &m_hits3D, "hits3D/I");
    m_pRecoTree->Branch("totalTime", &m_totalTime, "time/F");
    m_pRecoTree->Branch("artHitsTime", &m_artHitsTime, "time/F");
    m_pRecoTree->Branch("makeHitsTime", &m_makeHitsTime, "time/F");
    m_pRecoTree->Branch("buildneigborhoodTime", &m_buildNeighborhoodTime, "time/F");
    m_pRecoTree->Branch("dbscanTime", &m_dbscanTime, "time/F");
    m_pRecoTree->Branch("clusterMergeTime", &m_clusterMergeTime, "time/F");
    m_pRecoTree->Branch("pathfindingtime", &m_pathFindingTime, "time/F");
    m_pRecoTree->Branch("finishTime", &m_finishTime, "time/F");

    m_clusterPathAlg->initializeHistograms(*tfs.get());
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  Cluster3DSBN::PrepareEvent(const art::Event& evt)
  {
    m_run = evt.run();
    m_event = evt.id().event();
    m_hits = 0;
    m_hits3D = 0;
    m_totalTime = 0.f;
    m_artHitsTime = 0.f;
    m_makeHitsTime = 0.f;
    m_buildNeighborhoodTime = 0.f;
    m_dbscanTime = 0.f;
    m_pathFindingTime = 0.f;
    m_finishTime = 0.f;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void
  Cluster3DSBN::findTrackSeeds(art::Event& evt,
                                  reco::ClusterParameters& cluster,
                                  IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap,
                                  std::vector<recob::Seed>& seedVec,
                                  art::Assns<recob::Seed, recob::Hit>& seedHitAssns) const
  {
    /**
     *  @brief This method provides an interface to various algorithms for finding candiate
     *         recob::Seed objects and, as well, their candidate related seed hits
     */

    // Make sure we are using the right pca
    reco::PrincipalComponents& fullPCA = cluster.getFullPCA();
    reco::PrincipalComponents& skeletonPCA = cluster.getSkeletonPCA();
    reco::HitPairListPtr& hitPairListPtr = cluster.getHitPairListPtr();
    reco::HitPairListPtr skeletonListPtr;

    // We want to work with the "skeleton" hits so first step is to call the algorithm to
    // recover only these hits from the entire input collection
    m_skeletonAlg.GetSkeletonHits(hitPairListPtr, skeletonListPtr);

    // Skeleton hits are nice but we can do better if we then make a pass through to "average"
    // the skeleton hits position in the Y-Z plane
    m_skeletonAlg.AverageSkeletonPositions(skeletonListPtr);

    SeedHitPairListPairVec seedHitPairVec;

    // Some combination of the elements below will be used to determine which seed finding algorithm
    // to pursue below
    float eigenVal0 = 3. * sqrt(skeletonPCA.getEigenValues()[0]);
    float eigenVal1 = 3. * sqrt(skeletonPCA.getEigenValues()[1]);
    float eigenVal2 = 3. * sqrt(skeletonPCA.getEigenValues()[2]);
    float transRMS = std::hypot(eigenVal0, eigenVal1);

    bool foundGoodSeed(false);

    // Choose a method for finding the seeds based on the PCA that was run...
    // Currently we have an ad hoc if-else block which I hope will be improved soon!
    if (aParallelHitsCluster(fullPCA)) {
      // In this case we have a track moving relatively parallel to the wire plane with lots of
      // ambiguous 3D hits. Your best bet here is to use the "parallel hits" algorithm to get the
      // best axis and seeds
      // This algorithm does not fail (foundGoodSeed will always return true)
      foundGoodSeed = m_parallelHitsAlg.findTrackSeeds(hitPairListPtr, skeletonPCA, seedHitPairVec);
    }
    else if (eigenVal2 > 40. && transRMS < 5.) {
      // If the input cluster is relatively "straight" then chances are it is a single straight track,
      // probably a CR muon, and we can simply use the PCA to determine the seed
      // This algorithm will check both "ends" of the input hits and if the angles become inconsistent
      // then it will "fail"
      foundGoodSeed =
        m_pcaSeedFinderAlg.findTrackSeeds(skeletonListPtr, skeletonPCA, seedHitPairVec);
    }

    // In the event the above two methods failed then we hit it with the real seed finder
    if (!foundGoodSeed) {
      // If here then we have a complicated 3D cluster and we'll use the hough transform algorithm to
      // return a list of candidate seeds and seed hits
      m_seedFinderAlg.findTrackSeeds(skeletonListPtr, skeletonPCA, seedHitPairVec);
    }

    // Go through the returned lists and build out the art friendly seeds and hits
    for (const auto& seedHitPair : seedHitPairVec) {
      seedVec.push_back(seedHitPair.first);

      // We use a set here because our 3D hits can share 2D hits
      // The set will make sure we get unique combinations of 2D hits
      std::set<art::Ptr<recob::Hit>> seedHitSet;
      for (const auto& hit3D : seedHitPair.second) {
        for (const auto& hit2D : hit3D->getHits()) {
          if (!hit2D) continue;

          const recob::Hit* recobHit = hit2D->getHit();
          seedHitSet.insert(hitToPtrMap[recobHit]);
        }
      }

      RecobHitVector seedHitVec;
      for (const auto& hit2D : seedHitSet)
        seedHitVec.push_back(hit2D);

      util::CreateAssn(evt, seedVec, seedHitVec, seedHitAssns);
    }
  }

  struct Hit3DDistanceOrder {
    bool
    operator()(const std::pair<float, const reco::ClusterHit3D*>& left,
               const std::pair<float, const reco::ClusterHit3D*>& right)
    {
      return left.first < right.first;
    }
  };

  void
  Cluster3DSBN::splitClustersWithMST(reco::ClusterParameters& clusterParameters,
                                        reco::ClusterParametersList& clusterParametersList) const
  {
    // This is being left in place for future development. Essentially, it was an attempt to implement
    // a Minimum Spanning Tree as a way to split a particular cluster topology, one where two straight
    // tracks cross closely enought to appear as one cluster. As of Feb 2, 2015 I think the idea is still
    // worth merit so am leaving this module in place for now.
    //
    // If this routine is called then we believe we have a cluster which needs splitting.
    // The way we will do this is to use a Minimum Spanning Tree algorithm to associate all
    // hits together by their distance apart. In theory, we should be able to split the cluster
    // by finding the largest distance and splitting at that point.
    //
    // Typedef some data structures that we will use.
    // Start with the adjacency map
    typedef std::pair<float, const reco::ClusterHit3D*> DistanceHit3DPair;
    typedef std::list<DistanceHit3DPair> DistanceHit3DPairList;
    typedef std::map<const reco::ClusterHit3D*, DistanceHit3DPairList> Hit3DToDistanceMap;

    // Now typedef the lists we'll keep
    typedef std::list<const reco::ClusterHit3D*> Hit3DList;
    typedef std::pair<Hit3DList::iterator, Hit3DList::iterator> Hit3DEdgePair;
    typedef std::pair<float, Hit3DEdgePair> DistanceEdgePair;
    typedef std::list<DistanceEdgePair> DistanceEdgePairList;

    struct DistanceEdgePairOrder {
      bool
      operator()(const DistanceEdgePair& left, const DistanceEdgePair& right) const
      {
        return left.first > right.first;
      }
    };

    // Recover the hits we'll work on.
    // Note that we use on the skeleton hits so will need to recover them
    reco::HitPairListPtr& hitPairListPtr = clusterParameters.getHitPairListPtr();
    reco::HitPairListPtr skeletonListPtr;

    // We want to work with the "skeleton" hits so first step is to call the algorithm to
    // recover only these hits from the entire input collection
    m_skeletonAlg.GetSkeletonHits(hitPairListPtr, skeletonListPtr);

    // Skeleton hits are nice but we can do better if we then make a pass through to "average"
    // the skeleton hits position in the Y-Z plane
    m_skeletonAlg.AverageSkeletonPositions(skeletonListPtr);

    // First task is to define and build the adjacency map
    Hit3DToDistanceMap hit3DToDistanceMap;

    for (reco::HitPairListPtr::const_iterator hit3DOuterItr = skeletonListPtr.begin();
         hit3DOuterItr != skeletonListPtr.end();) {
      const reco::ClusterHit3D* hit3DOuter = *hit3DOuterItr++;
      DistanceHit3DPairList& outerHitList = hit3DToDistanceMap[hit3DOuter];
      TVector3 outerPos(
        hit3DOuter->getPosition()[0], hit3DOuter->getPosition()[1], hit3DOuter->getPosition()[2]);

      for (reco::HitPairListPtr::const_iterator hit3DInnerItr = hit3DOuterItr;
           hit3DInnerItr != skeletonListPtr.end();
           hit3DInnerItr++) {
        const reco::ClusterHit3D* hit3DInner = *hit3DInnerItr;
        TVector3 innerPos(
          hit3DInner->getPosition()[0], hit3DInner->getPosition()[1], hit3DInner->getPosition()[2]);
        TVector3 deltaPos = innerPos - outerPos;
        float hitDistance(float(deltaPos.Mag()));

        if (hitDistance > 20.) continue;

        hit3DToDistanceMap[hit3DInner].emplace_back(hitDistance, hit3DOuter);
        outerHitList.emplace_back(hitDistance, hit3DInner);
      }

      // Make sure our membership bit is clear
      hit3DOuter->clearStatusBits(reco::ClusterHit3D::SELECTEDBYMST);
    }

    // Make pass through again to order each of the lists
    for (auto& mapPair : hit3DToDistanceMap) {
      mapPair.second.sort(Hit3DDistanceOrder());
    }

    // Get the containers for the MST to operate on/with
    Hit3DList hit3DList;
    DistanceEdgePairList distanceEdgePairList;

    // Initialize with first element
    hit3DList.emplace_back(skeletonListPtr.front());
    distanceEdgePairList.emplace_back(
      DistanceEdgePair(0., Hit3DEdgePair(hit3DList.begin(), hit3DList.begin())));

    skeletonListPtr.front()->setStatusBit(reco::ClusterHit3D::SELECTEDBYMST);

    float largestDistance(0.);

    // Now run the MST
    // Basically, we loop until the MST list is the same size as the input list
    while (hit3DList.size() < skeletonListPtr.size()) {
      Hit3DList::iterator bestHit3DIter = hit3DList.begin();
      float bestDist = 10000000.;

      // Loop through all hits currently in the list and look for closest hit not in the list
      for (Hit3DList::iterator hit3DIter = hit3DList.begin(); hit3DIter != hit3DList.end();
           hit3DIter++) {
        const reco::ClusterHit3D* hit3D = *hit3DIter;

        // For the given 3D hit, find the closest to it that is not already in the list
        DistanceHit3DPairList& nearestList = hit3DToDistanceMap[hit3D];

        while (!nearestList.empty()) {
          const reco::ClusterHit3D* hit3DToCheck = nearestList.front().second;

          if (!hit3DToCheck->bitsAreSet(reco::ClusterHit3D::SELECTEDBYMST)) {
            if (nearestList.front().first < bestDist) {
              bestHit3DIter = hit3DIter;
              bestDist = nearestList.front().first;
            }
            break;
          }

          nearestList.pop_front();
        }
      }

      if (bestDist > largestDistance) largestDistance = bestDist;

      // Now we add the best hit not in the list to our list, keep track of the distance
      // to the object it was closest to
      const reco::ClusterHit3D* bestHit3D = *bestHit3DIter; // "best" hit already in the list
      const reco::ClusterHit3D* nextHit3D =
        hit3DToDistanceMap[bestHit3D].front().second; // "next" hit we are adding to the list

      Hit3DList::iterator nextHit3DIter = hit3DList.insert(hit3DList.end(), nextHit3D);

      distanceEdgePairList.emplace_back(bestDist, Hit3DEdgePair(bestHit3DIter, nextHit3DIter));

      nextHit3D->setStatusBit(reco::ClusterHit3D::SELECTEDBYMST);
    }

    float thirdDist = 2. * sqrt(clusterParameters.getSkeletonPCA().getEigenValues()[2]);

    // Ok, find the largest distance in the iterator map
    distanceEdgePairList.sort(DistanceEdgePairOrder());

    DistanceEdgePairList::iterator largestDistIter = distanceEdgePairList.begin();

    for (DistanceEdgePairList::iterator edgeIter = distanceEdgePairList.begin();
         edgeIter != distanceEdgePairList.end();
         edgeIter++) {
      if (edgeIter->first < thirdDist) break;

      largestDistIter = edgeIter;
    }

    reco::HitPairListPtr::iterator breakIter = largestDistIter->second.second;
    reco::HitPairListPtr bestList;

    bestList.resize(std::distance(hit3DList.begin(), breakIter));

    std::copy(hit3DList.begin(), breakIter, bestList.begin());

    // Remove from the grand hit list and see what happens...
    // The pieces below are incomplete and were really for testing only.
    hitPairListPtr.sort();
    bestList.sort();

    reco::HitPairListPtr::iterator newListEnd = std::set_difference(hitPairListPtr.begin(),
                                                                    hitPairListPtr.end(),
                                                                    bestList.begin(),
                                                                    bestList.end(),
                                                                    hitPairListPtr.begin());

    hitPairListPtr.erase(newListEnd, hitPairListPtr.end());
  }

  class CopyIfInRange {
  public:
    CopyIfInRange(float maxRange) : m_maxRange(maxRange) {}

    bool
    operator()(const reco::ClusterHit3D* hit3D) const
    {
      return hit3D->getDocaToAxis() < m_maxRange;
    }

  private:
    float m_maxRange;
  };

  void
  Cluster3DSBN::splitClustersWithHough(reco::ClusterParameters& clusterParameters,
                                          reco::ClusterParametersList& clusterParametersList) const
  {
    // @brief A method for splitted "crossed tracks" clusters into separate clusters
    //
    // If this routine is called then we believe we have a cluster which needs splitting.
    // The specific topology we are looking for is two long straight tracks which cross at some
    // point in close proximity so their hits were joined into a single 3D cluster. The method
    // to split this topology is to let the hough transform algorithm find the two leading candidates
    // and then to see if we use those to build two clusters instead of one.

    // Recover the hits we'll work on.
    // Note that we use on the skeleton hits so will need to recover them
    reco::HitPairListPtr& hitPairListPtr = clusterParameters.getHitPairListPtr();
    reco::HitPairListPtr skeletonListPtr;

    // We want to work with the "skeleton" hits so first step is to call the algorithm to
    // recover only these hits from the entire input collection
    m_skeletonAlg.GetSkeletonHits(hitPairListPtr, skeletonListPtr);

    // Skeleton hits are nice but we can do better if we then make a pass through to "average"
    // the skeleton hits position in the Y-Z plane
    m_skeletonAlg.AverageSkeletonPositions(skeletonListPtr);

    // Define the container for our lists of hits
    reco::HitPairListPtrList hitPairListPtrList;

    // Now feed this to the Hough Transform to find candidate straight lines
    m_seedFinderAlg.findTrackHits(
      skeletonListPtr, clusterParameters.getSkeletonPCA(), hitPairListPtrList);

    // We need at least two lists or else there is nothing to do
    if (hitPairListPtrList.size() < 2) return;

    // The game plan will be the following:
    // 1) Take the first list of hits and run the PCA on this to get an axis
    //    - Then calculate the 3d doca for ALL hits in the cluster to this axis
    //    - Move all hits within "3 sigam" of the axis to a new list
    // 2) run the PCA on the second list of hits to get that axis
    //    - Then calculate the 3d doca for all hits in our first list
    //    - Copy hits in the first list which are within 3 sigma of the new axis
    //      back into our original cluster - these are shared hits
    reco::HitPairListPtrList::iterator hitPairListIter = hitPairListPtrList.begin();
    reco::HitPairListPtr& firstHitList = *hitPairListIter++;
    reco::PrincipalComponents firstHitListPCA;

    m_pcaAlg.PCAAnalysis_3D(firstHitList, firstHitListPCA);

    // Make sure we have a successful calculation.
    if (firstHitListPCA.getSvdOK()) {
      // The fill routines below will expect to see unused 2D hits so we need to clear the
      // status bits... and I am not sure of a better way...
      for (const auto& hit3D : hitPairListPtr) {
        for (const auto& hit2D : hit3D->getHits())
          if (hit2D) hit2D->clearStatusBits(0x1);
      }

      // Calculate the 3D doca's for the hits which were used to make this PCA
      m_pcaAlg.PCAAnalysis_calc3DDocas(firstHitList, firstHitListPCA);

      // Divine from the ether some maximum allowed range for transfering hits
      float allowedHitRange = 6. * firstHitListPCA.getAveHitDoca();

      // Now go through and calculate the 3D doca's for ALL the hits in the original cluster
      m_pcaAlg.PCAAnalysis_calc3DDocas(hitPairListPtr, firstHitListPCA);

      // Let's make a new cluster to hold the hits
      clusterParametersList.push_back(reco::ClusterParameters());

      // Can we get a reference to what we just created?
      reco::ClusterParameters& newClusterParams = clusterParametersList.back();

      reco::HitPairListPtr& newClusterHitList = newClusterParams.getHitPairListPtr();

      newClusterHitList.resize(hitPairListPtr.size());

      // Do the actual copy of the hits we want
      reco::HitPairListPtr::iterator newListEnd = std::copy_if(hitPairListPtr.begin(),
                                                               hitPairListPtr.end(),
                                                               newClusterHitList.begin(),
                                                               CopyIfInRange(allowedHitRange));

      // Shrink to fit
      newClusterHitList.resize(std::distance(newClusterHitList.begin(), newListEnd));

      // And now remove these hits from the original cluster
      hitPairListPtr.remove_if(CopyIfInRange(allowedHitRange));

      // Get an empty hit to cluster map...
      reco::Hit2DToClusterMap hit2DToClusterMap;

      // Now "fill" the cluster parameters but turn off the hit rejection
      m_clusterBuilder->FillClusterParams(newClusterParams, hit2DToClusterMap, 0., 1.);

      // Set the skeleton pca to the value calculated above on input
      clusterParameters.getSkeletonPCA() = firstHitListPCA;

      // We are done with splitting out one track. Because the two tracks cross in
      // close proximity, this is the one case where we might consider sharing 3D hits
      // So let's make a little detour here to try to copy some of those hits back into
      // the main hit list
      reco::HitPairListPtr& secondHitList = *hitPairListIter;
      reco::PrincipalComponents secondHitListPCA;

      m_pcaAlg.PCAAnalysis_3D(secondHitList, secondHitListPCA);

      // Make sure we have a successful calculation.
      if (secondHitListPCA.getSvdOK()) {
        // Calculate the 3D doca's for the hits which were used to make this PCA
        m_pcaAlg.PCAAnalysis_calc3DDocas(secondHitList, secondHitListPCA);

        // Since this is the "other" cluster, we'll be a bit more generous in adding back hits
        float newAllowedHitRange = 6. * secondHitListPCA.getAveHitDoca();

        // Go through and calculate the 3D doca's for the hits in our new candidate cluster
        m_pcaAlg.PCAAnalysis_calc3DDocas(newClusterHitList, secondHitListPCA);

        // Create a temporary list to fill with the hits we might want to save
        reco::HitPairListPtr tempHitList(newClusterHitList.size());

        // Do the actual copy of the hits we want...
        reco::HitPairListPtr::iterator tempListEnd =
          std::copy_if(newClusterHitList.begin(),
                       newClusterHitList.end(),
                       tempHitList.begin(),
                       CopyIfInRange(newAllowedHitRange));

        hitPairListPtr.insert(hitPairListPtr.end(), tempHitList.begin(), tempListEnd);
      }

      // Of course, now we need to modify the original cluster parameters
      reco::ClusterParameters originalParams(hitPairListPtr);

      // Now "fill" the cluster parameters but turn off the hit rejection
      m_clusterBuilder->FillClusterParams(originalParams, hit2DToClusterMap, 0., 1.);

      // Overwrite original cluster parameters with our new values
      clusterParameters.getClusterParams() = originalParams.getClusterParams();
      clusterParameters.getFullPCA() = originalParams.getFullPCA();
      clusterParameters.getSkeletonPCA() = secondHitListPCA;
    }
  }

  void
  Cluster3DSBN::ProduceArtClusters(util::GeometryUtilities const& gser,
                                      ArtOutputHandler& output,
                                      reco::HitPairList& hitPairVector,
                                      reco::ClusterParametersList& clusterParametersList,
                                      IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap) const
  {
    /**
     *  @brief The workhorse to take the candidate 3D clusters and produce all of the necessary art output
     */

    mf::LogDebug("Cluster3DSBN") << " *** Cluster3DSBN::ProduceArtClusters() *** " << std::endl;

    // Make sure there is something to do here!
    if (!clusterParametersList.empty()) {
      // This is the loop over candidate 3D clusters
      // Note that it might be that the list of candidate clusters is modified by splitting
      // So we use the following construct to make sure we get all of them
      for (auto& clusterParameters : clusterParametersList) {
        // It should be straightforward at this point to transfer information from our vector of clusters
        // to the larsoft objects... of course we still have some work to do first, in particular to
        // find the candidate seeds and their seed hits

        // The chances of getting here and this condition not being true are probably zero... but check anyway
        if (!clusterParameters.getFullPCA().getSvdOK()) {
          mf::LogDebug("Cluster3DSBN")
            << "--> no feature extraction done on this cluster!!" << std::endl;
          continue;
        }

        // Keep track of hit 3D to SP for when we do edges
        Hit3DToSPPtrMap hit3DToSPPtrMap;
        Hit3DToSPPtrMap best3DToSPPtrMap;

        // Do a special output of voronoi vertices here...
        dcel2d::VertexList& vertexList = clusterParameters.getVertexList();
        dcel2d::HalfEdgeList& halfEdgeList = clusterParameters.getHalfEdgeList();

        MakeAndSaveVertexPoints(output, vertexList, halfEdgeList);

        // Special case handling... if no daughters then call standard conversion routine to make sure space points
        // created, etc.
        if (clusterParameters.daughterList().empty()) {
          ConvertToArtOutput(gser,
                             output,
                             clusterParameters,
                             recob::PFParticle::kPFParticlePrimary,
                             hitToPtrMap,
                             hit3DToSPPtrMap,
                             best3DToSPPtrMap);

          // Get the extreme points
          MakeAndSaveKinkPoints(output,
                                clusterParameters.getConvexHull().getConvexHullKinkPoints());
        }
        // Otherwise, the cluster has daughters so we handle specially
        else {
          // Set up to keep track of parent/daughters
          IdxToPCAMap idxToPCAMap;
          size_t numTotalDaughters = countUltimateDaughters(clusterParameters);
          size_t pfParticleIdx(output.artPFParticleVector->size() + numTotalDaughters);

          FindAndStoreDaughters(gser,
                                output,
                                clusterParameters,
                                pfParticleIdx,
                                idxToPCAMap,
                                hitToPtrMap,
                                hit3DToSPPtrMap,
                                best3DToSPPtrMap);

          // Need to make a daughter vec from our map
          std::vector<size_t> daughterVec;

          for (auto& idxToPCA : idxToPCAMap)
            daughterVec.emplace_back(idxToPCA.first);

          // Now create/handle the parent PFParticle
          recob::PFParticle pfParticle(
            13, pfParticleIdx, recob::PFParticle::kPFParticlePrimary, daughterVec);
          output.artPFParticleVector->push_back(pfParticle);

          recob::PCAxis::EigenVectors eigenVecs;
          double eigenVals[] = {0., 0., 0.};
          double avePosition[] = {0., 0., 0.};

          eigenVecs.resize(3);

          reco::PrincipalComponents& skeletonPCA = clusterParameters.getSkeletonPCA();

          for (size_t outerIdx = 0; outerIdx < 3; outerIdx++) {
            avePosition[outerIdx] = skeletonPCA.getAvePosition()[outerIdx];
            eigenVals[outerIdx] = skeletonPCA.getEigenValues()[outerIdx];

            eigenVecs[outerIdx].resize(3);

            // Be careful here... eigen stores in column major order buy default
            for (size_t innerIdx = 0; innerIdx < 3; innerIdx++)
              eigenVecs[outerIdx][innerIdx] = skeletonPCA.getEigenVectors().row(outerIdx)[innerIdx];
          }

          recob::PCAxis skelPcAxis(skeletonPCA.getSvdOK(),
                                   skeletonPCA.getNumHitsUsed(),
                                   eigenVals,
                                   eigenVecs,
                                   avePosition,
                                   skeletonPCA.getAveHitDoca(),
                                   output.artPCAxisVector->size());

          output.artPCAxisVector->push_back(skelPcAxis);

          reco::PrincipalComponents& fullPCA = clusterParameters.getFullPCA();

          for (size_t outerIdx = 0; outerIdx < 3; outerIdx++) {
            avePosition[outerIdx] = fullPCA.getAvePosition()[outerIdx];
            eigenVals[outerIdx] = fullPCA.getEigenValues()[outerIdx];

            for (size_t innerIdx = 0; innerIdx < 3; innerIdx++)
              eigenVecs[outerIdx][innerIdx] = fullPCA.getEigenVectors().row(outerIdx)(innerIdx);
          }

          recob::PCAxis fullPcAxis(fullPCA.getSvdOK(),
                                   fullPCA.getNumHitsUsed(),
                                   eigenVals,
                                   eigenVecs,
                                   avePosition,
                                   fullPCA.getAveHitDoca(),
                                   output.artPCAxisVector->size());

          output.artPCAxisVector->push_back(fullPcAxis);

          // Create associations to the PFParticle
          output.makePFPartPCAAssns();

          // Make associations to all space points for this cluster
          MakeAndSaveSpacePoints(output,
                                 *output.artSpacePointVector.get(),
                                 *output.artSPHitAssociations.get(),
                                 clusterParameters.getHitPairListPtr(),
                                 hitToPtrMap,
                                 hit3DToSPPtrMap);

          // Get the extreme points
          MakeAndSaveKinkPoints(output,
                                clusterParameters.getConvexHull().getConvexHullKinkPoints());

          // Build the edges now
          size_t edgeStart(output.artEdgeVector->size());

          for (const auto& edge : clusterParameters.getConvexHull().getConvexHullEdgeList()) {
            RecobSpacePointVector spacePointVec;

            try {
              spacePointVec.push_back(hit3DToSPPtrMap.at(std::get<0>(edge)));
              spacePointVec.push_back(hit3DToSPPtrMap.at(std::get<1>(edge)));
            }
            catch (...) {
              std::cout << "Caught exception in looking up space point ptr... " << std::get<0>(edge)
                        << ", " << std::get<1>(edge) << std::endl;
              continue;
            }

            output.artEdgeVector->emplace_back(std::get<2>(edge),
                                               spacePointVec[0].key(),
                                               spacePointVec[1].key(),
                                               output.artEdgeVector->size());

            output.makeEdgeSpacePointAssns(
              *output.artEdgeVector.get(), spacePointVec, *output.artEdgeSPAssociations.get());
          }

          output.makePFPartEdgeAssns(
            *output.artEdgeVector.get(), *output.artPFPartEdgeAssociations.get(), edgeStart);
        }
      }
    }

    mf::LogDebug("Cluster3DSBN") << "--> Preparing to save the remaining space points, there are " << hitPairVector.size() << std::endl;

    // Right now error matrix is uniform...
    int nFreePoints(0);

    // Run through the HitPairVector and add any unused hit pairs to the list
    for (auto& hitPair : hitPairVector) 
    {
        if (hitPair.bitsAreSet(reco::ClusterHit3D::MADESPACEPOINT)) continue;

        double spacePointPos[] = {
          hitPair.getPosition()[0], hitPair.getPosition()[1], hitPair.getPosition()[2]};
        double spacePointErr[] = {1., 0., 0., 1., 0., 1.};
        double chisq = hitPair.getHitChiSquare();

        RecobHitVector recobHits;

        for (const auto hit : hitPair.getHits()) {
          if (!hit) continue;

          art::Ptr<recob::Hit> hitPtr = hitToPtrMap[hit->getHit()];
          recobHits.push_back(hitPtr);
        }

        nFreePoints++;

        spacePointErr[1] = hitPair.getTotalCharge();
        spacePointErr[3] = hitPair.getChargeAsymmetry();

        output.artSpacePointVector->push_back(
          recob::SpacePoint(spacePointPos, spacePointErr, chisq, output.artSpacePointVector->size()));

        if (!recobHits.empty())
          output.makeSpacePointHitAssns(
            *output.artSpacePointVector.get(), recobHits, *output.artSPHitAssociations.get());
    } 

    mf::LogDebug("Cluster3DSBN") << "++++>>>> total num hits: " << hitPairVector.size()
              << ", num free: " << nFreePoints;
  }

  size_t
  Cluster3DSBN::countUltimateDaughters(reco::ClusterParameters& clusterParameters) const
  {
    size_t localCount(0);

    if (!clusterParameters.daughterList().empty()) {
      for (auto& clusterParams : clusterParameters.daughterList())
        localCount += countUltimateDaughters(clusterParams);
    }
    else
      localCount++;

    return localCount;
  }

  size_t
  Cluster3DSBN::FindAndStoreDaughters(util::GeometryUtilities const& gser,
                                         ArtOutputHandler& output,
                                         reco::ClusterParameters& clusterParameters,
                                         size_t pfParticleParent,
                                         IdxToPCAMap& idxToPCAMap,
                                         IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap,
                                         Hit3DToSPPtrMap& hit3DToSPPtrMap,
                                         Hit3DToSPPtrMap& best3DToSPPtrMap) const
  {
    // This is a recursive routine, we keep calling ourself as long as the daughter list is non empty
    if (!clusterParameters.daughterList().empty()) {
      for (auto& clusterParams : clusterParameters.daughterList())
        FindAndStoreDaughters(gser,
                              output,
                              clusterParams,
                              pfParticleParent,
                              idxToPCAMap,
                              hitToPtrMap,
                              hit3DToSPPtrMap,
                              best3DToSPPtrMap);
    }
    // Otherwise we want to store the information
    else {
      size_t daughterIdx = ConvertToArtOutput(gser,
                                              output,
                                              clusterParameters,
                                              pfParticleParent,
                                              hitToPtrMap,
                                              hit3DToSPPtrMap,
                                              best3DToSPPtrMap);

      idxToPCAMap[daughterIdx] = &clusterParameters.getFullPCA();
    }

    return idxToPCAMap.size();
  }

  size_t
  Cluster3DSBN::ConvertToArtOutput(util::GeometryUtilities const& gser,
                                      ArtOutputHandler& output,
                                      reco::ClusterParameters& clusterParameters,
                                      size_t pfParticleParent,
                                      IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap,
                                      Hit3DToSPPtrMap& hit3DToSPPtrMap,
                                      Hit3DToSPPtrMap& best3DToSPPtrMap) const
  {
    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here, except that we override selected items
    // (so, thanks to metaprogramming, we finally have wrappers of wrappers);
    // configuration would happen here, but we are using the default
    // configuration for that algorithm
    using OverriddenClusterParamsAlg_t =
      cluster::OverriddenClusterParamsAlg<cluster::StandardClusterParamsAlg>;

    cluster::ClusterParamsImportWrapper<OverriddenClusterParamsAlg_t> ClusterParamAlgo;

    // It should be straightforward at this point to transfer information from our vector of clusters
    // to the larsoft objects... of course we still have some work to do first, in particular to
    // find the candidate seeds and their seed hits

    // We keep track of 2 PCA axes, the first is the "full" PCA run over all the 3D hits in the
    // candidate cluster. The second will be that derived from just using the "skeleton" hits.
    // Make a copy of the full PCA to keep that, then get a reference for the skeleton PCA
    reco::PrincipalComponents& fullPCA = clusterParameters.getFullPCA();
    reco::PrincipalComponents& skeletonPCA = clusterParameters.getSkeletonPCA();

    size_t clusterStart = output.artClusterVector->size();

    // Start loop over views to build out the hit lists and the 2D cluster objects
    for (const auto& clusParametersPair : clusterParameters.getClusterParams()) {
      const reco::RecobClusterParameters& clusParams = clusParametersPair.second;

      // Protect against a missing view
      if (clusParams.m_view == geo::kUnknown) continue;

      // Get a vector of hit pointers to seed the cluster algorithm
      std::vector<const recob::Hit*> recobHitVec(clusParams.m_hitVector.size());

      std::transform(clusParams.m_hitVector.begin(),
                     clusParams.m_hitVector.end(),
                     recobHitVec.begin(),
                     [](const auto& hit2D) { return hit2D->getHit(); });

      // Get the tdc/wire slope... from the unit vector...
      double startWire(clusParams.m_startWire);
      double endWire(clusParams.m_endWire);
      double startTime(clusParams.m_startTime);
      double endTime(clusParams.m_endTime);

      // feed the algorithm with all the cluster hits
      ClusterParamAlgo.ImportHits(gser, recobHitVec);

      // create the recob::Cluster directly in the vector
      cluster::ClusterCreator artCluster(gser,
                                         ClusterParamAlgo,                // algo
                                         startWire,                       // start_wire
                                         0.,                              // sigma_start_wire
                                         startTime,                       // start_tick
                                         clusParams.m_sigmaStartTime,     // sigma_start_tick
                                         endWire,                         // end_wire
                                         0.,                              // sigma_end_wire,
                                         endTime,                         // end_tick
                                         clusParams.m_sigmaEndTime,       // sigma_end_tick
                                         output.artClusterVector->size(), // ID
                                         clusParams.m_view,               // view
                                         clusParams.m_plane,              // plane
                                         recob::Cluster::Sentry           // sentry
      );

      output.artClusterVector->emplace_back(artCluster.move());

      // We love looping. In this case, our list of hits is comprised of "ClusterHits" and we need to get a RecobHitVector instead...
      RecobHitVector recobHits;

      for (const auto& hit2D : clusParams.m_hitVector) {
        if (hit2D == nullptr || hit2D->getHit() == nullptr) continue;

        art::Ptr<recob::Hit> hitPtr = hitToPtrMap[hit2D->getHit()];
        recobHits.push_back(hitPtr);
      }

      output.makeClusterHitAssns(recobHits);
    }

    // Last, let's try to get seeds for tracking..
    // Keep track of how many we have so far
    size_t numSeedsStart = output.artSeedVector->size();

    // Empty daughter vector for now
    std::vector<size_t> nullVector;
    size_t pfParticleIdx(output.artPFParticleVector->size());

    recob::PFParticle pfParticle(13, pfParticleIdx, pfParticleParent, nullVector);
    output.artPFParticleVector->push_back(pfParticle);

    // Assume that if we are a daughter particle then there is a set of "best" 3D points and the convex hull will have been calculated from them
    // If there is no best list then assume the convex hull is from all of the points...
    std::string instance("");
    reco::HitPairListPtr* hit3DListPtr = &clusterParameters.getHitPairListPtr();
    Hit3DToSPPtrMap* hit3DToPtrMap = &hit3DToSPPtrMap;

    auto spVector = output.artSpacePointVector.get();
    auto edgeVector = output.artEdgeVector.get();
    auto pfPartEdgeAssns = output.artPFPartEdgeAssociations.get();
    auto spEdgeAssns = output.artEdgeSPAssociations.get();
    auto spHitAssns = output.artSPHitAssociations.get();
    auto pfPartSPAssns = output.artPFPartSPAssociations.get();

    // Make associations to all space points for this cluster
    int spacePointStart = spVector->size();

    MakeAndSaveSpacePoints(
      output, *spVector, *spHitAssns, *hit3DListPtr, hitToPtrMap, *hit3DToPtrMap);

    // And make sure to have associations to the PFParticle
    output.makePFPartSpacePointAssns(*spVector, *pfPartSPAssns, spacePointStart);

    // If they exist, make space points for the Path points
    if (!clusterParameters.getBestHitPairListPtr().empty()) {
      instance = m_pathInstance;
      hit3DListPtr = &clusterParameters.getBestHitPairListPtr();
      hit3DToPtrMap = &best3DToSPPtrMap;

      spVector = output.artPathPointVector.get();
      edgeVector = output.artPathEdgeVector.get();
      pfPartEdgeAssns = output.artPFPartPathEdgeAssociations.get();
      spEdgeAssns = output.artEdgePPAssociations.get();
      spHitAssns = output.artPPHitAssociations.get();
      pfPartSPAssns = output.artPFPartPPAssociations.get();

      spacePointStart = spVector->size();

      MakeAndSaveSpacePoints(
        output, *spVector, *spHitAssns, *hit3DListPtr, hitToPtrMap, *hit3DToPtrMap, instance);

      // And make sure to have associations to the PFParticle
      output.makePFPartSpacePointAssns(*spVector, *pfPartSPAssns, spacePointStart, instance);
    }

    // Now handle the edges according to whether associated with regular or "best" space points
    if (!clusterParameters.getBestEdgeList().empty()) {
      size_t edgeStart = edgeVector->size();

      for (const auto& edge : clusterParameters.getBestEdgeList()) {
        RecobSpacePointVector spacePointVec;

        try {
          spacePointVec.push_back(hit3DToPtrMap->at(std::get<0>(edge)));
          spacePointVec.push_back(hit3DToPtrMap->at(std::get<1>(edge)));
        }
        catch (...) {
          std::cout << "Caught exception in looking up space point ptr... " << std::get<0>(edge)
                    << ", " << std::get<1>(edge) << std::endl;
          continue;
        }

        edgeVector->emplace_back(
          std::get<2>(edge), spacePointVec[0].key(), spacePointVec[1].key(), edgeVector->size());

        output.makeEdgeSpacePointAssns(*edgeVector, spacePointVec, *spEdgeAssns, instance);
      }

      output.makePFPartEdgeAssns(*edgeVector, *pfPartEdgeAssns, edgeStart, instance);
    }

    // Look at making the PCAxis and associations - for both the skeleton (the first) and the full
    // First need some float to double conversion containers
    recob::PCAxis::EigenVectors eigenVecs;
    double eigenVals[] = {0., 0., 0.};
    double avePosition[] = {0., 0., 0.};

    eigenVecs.resize(3);

    for (size_t outerIdx = 0; outerIdx < 3; outerIdx++) {
      avePosition[outerIdx] = skeletonPCA.getAvePosition()[outerIdx];
      eigenVals[outerIdx] = skeletonPCA.getEigenValues()[outerIdx];

      eigenVecs[outerIdx].resize(3);

      for (size_t innerIdx = 0; innerIdx < 3; innerIdx++)
        eigenVecs[outerIdx][innerIdx] = skeletonPCA.getEigenVectors().row(outerIdx)(innerIdx);
    }

    recob::PCAxis skelPcAxis(skeletonPCA.getSvdOK(),
                             skeletonPCA.getNumHitsUsed(),
                             eigenVals,
                             eigenVecs,
                             avePosition,
                             skeletonPCA.getAveHitDoca(),
                             output.artPCAxisVector->size());

    output.artPCAxisVector->push_back(skelPcAxis);

    for (size_t outerIdx = 0; outerIdx < 3; outerIdx++) {
      avePosition[outerIdx] = fullPCA.getAvePosition()[outerIdx];
      eigenVals[outerIdx] = fullPCA.getEigenValues()[outerIdx];

      for (size_t innerIdx = 0; innerIdx < 3; innerIdx++)
        eigenVecs[outerIdx][innerIdx] = fullPCA.getEigenVectors().row(outerIdx)(innerIdx);
    }

    recob::PCAxis fullPcAxis(fullPCA.getSvdOK(),
                             fullPCA.getNumHitsUsed(),
                             eigenVals,   //fullPCA.getEigenValues(),
                             eigenVecs,   //fullPCA.getEigenVectors(),
                             avePosition, //fullPCA.getAvePosition(),
                             fullPCA.getAveHitDoca(),
                             output.artPCAxisVector->size());

    output.artPCAxisVector->push_back(fullPcAxis);

    // Create associations to the PFParticle
    output.makePFPartPCAAssns();
    output.makePFPartSeedAssns(numSeedsStart);
    output.makePFPartClusterAssns(clusterStart);

    return pfParticleIdx;
  }

  void
  Cluster3DSBN::MakeAndSaveSpacePoints(ArtOutputHandler& output,
                                          std::vector<recob::SpacePoint>& spacePointVec,
                                          art::Assns<recob::Hit, recob::SpacePoint>& spHitAssns,
                                          reco::HitPairListPtr& clusHitPairVector,
                                          IHit3DBuilder::RecobHitToPtrMap& hitToPtrMap,
                                          Hit3DToSPPtrMap& hit3DToSPPtrMap,
                                          const std::string& instance) const
  {
    // Reserve space...
    spacePointVec.reserve(spacePointVec.size() + clusHitPairVector.size());

    // Right now error matrix is uniform...
    double spError[] = {0., 0., 0., 0., 0., 0.};

    // Copy these hits to the vector to be stored with the event
    for (auto& hitPair : clusHitPairVector) {
      // Skip those space points that have already been created
      if (hit3DToSPPtrMap.find(hitPair) != hit3DToSPPtrMap.end()) continue;

      // Don't make space point if this hit was "rejected"
      if (hitPair->bitsAreSet(reco::ClusterHit3D::REJECTEDHIT)) continue;

      double chisq = hitPair->getHitChiSquare(); // secret handshake...

      // Mark this hit pair as in use
      hitPair->setStatusBit(reco::ClusterHit3D::MADESPACEPOINT);

      // Create and store the space point
      size_t spacePointID = spacePointVec.size();
      double spacePointPos[] = {
        hitPair->getPosition()[0], hitPair->getPosition()[1], hitPair->getPosition()[2]};

      spError[1] = hitPair->getTotalCharge();
      spError[3] = hitPair->getChargeAsymmetry();

      spacePointVec.emplace_back(spacePointPos, spError, chisq, spacePointID);

      // Update mapping
      hit3DToSPPtrMap[hitPair] = output.makeSpacePointPtr(spacePointID, instance);

      // space point hits associations
      RecobHitVector recobHits;

      for (const auto& hit : hitPair->getHits()) {
        if (!hit) continue;

        art::Ptr<recob::Hit> hitPtr = hitToPtrMap[hit->getHit()];
        recobHits.push_back(hitPtr);
      }

      if (!recobHits.empty())
        output.makeSpacePointHitAssns(spacePointVec, recobHits, spHitAssns, instance);
    }

    return;
  }

  void
  Cluster3DSBN::MakeAndSaveKinkPoints(ArtOutputHandler& output,
                                         reco::ConvexHullKinkTupleList& kinkTupleVec) const
  {
    // Right now error matrix is uniform...
    double spError[] = {1., 0., 1., 0., 0., 1.};

    // Copy these hits to the vector to be stored with the event
    for (auto& kinkTuple : kinkTupleVec) {
      const reco::ClusterHit3D* hit = std::get<2>(std::get<0>(kinkTuple));

      double chisq = hit->getHitChiSquare(); // secret handshake...

      // Create and store the space point
      double spacePointPos[] = {
        hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]};

      spError[1] = hit->getTotalCharge();
      spError[3] = hit->getChargeAsymmetry();

      output.artExtremePointVector->push_back(
        recob::SpacePoint(spacePointPos, spError, chisq, output.artExtremePointVector->size()));
    }

    return;
  }

  void
  Cluster3DSBN::MakeAndSaveVertexPoints(ArtOutputHandler& output,
                                           dcel2d::VertexList& vertexList,
                                           dcel2d::HalfEdgeList& halfEdgeList) const
  {
    // We actually do two things here:
    // 1) Create space points to represent the vertex locations of the voronoi diagram
    // 2) Create the edges that link the space points together

    // Set up the space point creation
    // Right now error matrix is uniform...
    double spError[] = {1., 0., 1., 0., 0., 1.};
    double chisq = 1.;

    // Keep track of the vertex to space point association
    std::map<const dcel2d::Vertex*, size_t> vertexToSpacePointMap;

    // Copy these hits to the vector to be stored with the event
    for (auto& vertex : vertexList) {
      const dcel2d::Coords& coords = vertex.getCoords();

      // Create and store the space point
      double spacePointPos[] = {coords[0], coords[1], coords[2]};

      vertexToSpacePointMap[&vertex] = output.artVertexPointVector->size();

      output.artVertexPointVector->emplace_back(
        spacePointPos, spError, chisq, output.artVertexPointVector->size());
    }

    // Try to avoid double counting
    std::set<const dcel2d::HalfEdge*> halfEdgeSet;

    // Build the edges now
    for (const auto& halfEdge : halfEdgeList) {
      // Recover twin
      const dcel2d::HalfEdge* twin = halfEdge.getTwinHalfEdge();

      // It can happen that we have no twin... and also check that we've not been here before
      if (twin && halfEdgeSet.find(twin) == halfEdgeSet.end()) {
        // Recover the vertices
        const dcel2d::Vertex* fromVertex = twin->getTargetVertex();
        const dcel2d::Vertex* toVertex = halfEdge.getTargetVertex();

        // It can happen for the open edges that there is no target vertex
        if (!toVertex || !fromVertex) continue;

        if (vertexToSpacePointMap.find(fromVertex) == vertexToSpacePointMap.end() ||
            vertexToSpacePointMap.find(toVertex) == vertexToSpacePointMap.end())
          continue;

        // Need the distance between vertices
        Eigen::Vector3f distVec = toVertex->getCoords() - fromVertex->getCoords();

        output.artVertexEdgeVector->emplace_back(distVec.norm(),
                                                 vertexToSpacePointMap[fromVertex],
                                                 vertexToSpacePointMap[toVertex],
                                                 output.artEdgeVector->size());

        halfEdgeSet.insert(&halfEdge);
      }
    }

    return;
  }

  void
  Cluster3DSBN::MakeAndSavePCAPoints(ArtOutputHandler& output,
                                        const reco::PrincipalComponents& fullPCA,
                                        IdxToPCAMap& idxToPCAMap) const
  {
    // We actually do two things here:
    // 1) Create space points from the centroids of the PCA for each cluster
    // 2) Create the edges that link the space points together

    // The first task is to put the list of PCA's into some semblance of order... they may be
    // preordered by likely they are piecewise ordered so fix that here

    // We'll need the current PCA axis to determine doca and arclen
    Eigen::Vector3f avePosition(
      fullPCA.getAvePosition()[0], fullPCA.getAvePosition()[1], fullPCA.getAvePosition()[2]);
    Eigen::Vector3f axisDirVec(fullPCA.getEigenVectors().row(2));

    using DocaToPCAPair = std::pair<float, const reco::PrincipalComponents*>;
    using DocaToPCAVec = std::vector<DocaToPCAPair>;

    DocaToPCAVec docaToPCAVec;

    // Outer loop over views
    for (const auto& idxToPCA : idxToPCAMap) {
      const reco::PrincipalComponents* pca = idxToPCA.second;

      // Now we need to calculate the doca and poca...
      // Start by getting this hits position
      Eigen::Vector3f pcaPos(
        pca->getAvePosition()[0], pca->getAvePosition()[1], pca->getAvePosition()[2]);

      // Form a TVector from this to the cluster average position
      Eigen::Vector3f pcaToAveVec = pcaPos - avePosition;

      // With this we can get the arclength to the doca point
      float arclenToPoca = pcaToAveVec.dot(axisDirVec);

      docaToPCAVec.emplace_back(arclenToPoca, pca);
    }

    std::sort(docaToPCAVec.begin(), docaToPCAVec.end(), [](const auto& left, const auto& right) {
      return left.first < right.first;
    });

    // Set up the space point creation
    // Right now error matrix is uniform...
    double spError[] = {1., 0., 1., 0., 0., 1.};
    double chisq = 1.;

    const reco::PrincipalComponents* lastPCA(nullptr);

    // Set up to loop through the clusters
    for (const auto& docaToPCAPair : docaToPCAVec) {
      // Recover the PCA for this cluster
      const reco::PrincipalComponents* curPCA = docaToPCAPair.second;

      if (lastPCA) {
        double lastPointPos[] = {
          lastPCA->getAvePosition()[0], lastPCA->getAvePosition()[1], lastPCA->getAvePosition()[2]};
        size_t lastPointBin = output.artVertexPointVector->size();
        double curPointPos[] = {
          curPCA->getAvePosition()[0], curPCA->getAvePosition()[1], curPCA->getAvePosition()[2]};
        size_t curPointBin = lastPointBin + 1;

        output.artVertexPointVector->emplace_back(lastPointPos, spError, chisq, lastPointBin);
        output.artVertexPointVector->emplace_back(curPointPos, spError, chisq, curPointBin);

        Eigen::Vector3f distVec(curPointPos[0] - lastPointPos[0],
                                curPointPos[1] - lastPointPos[1],
                                curPointPos[2] - lastPointPos[2]);

        output.artVertexEdgeVector->emplace_back(
          distVec.norm(), lastPointBin, curPointBin, output.artEdgeVector->size());
      }

      lastPCA = curPCA;
    }
  }

} // namespace lar_Cluster3DSBN
