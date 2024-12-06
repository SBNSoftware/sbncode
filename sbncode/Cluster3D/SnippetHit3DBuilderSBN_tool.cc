/**
 *  @file   SnippetHit3DBuilderSBN_tool.cc
 *
 *  @brief  This tool provides "standard" 3D hits built (by this tool) from 2D hits
 *
 */

// Framework Includes
#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IHit3DBuilder.h"

// Eigen
#include <Eigen/Core>

// std includes
#include <string>
#include <iostream>
#include <memory>
#include <numeric> // std::accumulate

// Ack!
#include "TTree.h"

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

/**
 *   @brief What follows are several highly useful typedefs which we
 *          want to expose to the outside world
 */

// forward declaration to define an ordering function for our hit set
struct Hit2DSetCompare
{
    bool operator() (const reco::ClusterHit2D*, const reco::ClusterHit2D*) const;
};

using HitVector                    = std::vector<const reco::ClusterHit2D*>;
using HitStartEndPair              = std::pair<raw::TDCtick_t,raw::TDCtick_t>;
using SnippetHitMap                = std::map<HitStartEndPair,HitVector>;
using PlaneToSnippetHitMap         = std::map<geo::PlaneID, SnippetHitMap>;
using TPCToPlaneToSnippetHitMap    = std::map<geo::TPCID, PlaneToSnippetHitMap>;
using Hit2DList                    = std::list<reco::ClusterHit2D>;
using Hit2DSet                     = std::set<const reco::ClusterHit2D*, Hit2DSetCompare>;
using WireToHitSetMap              = std::map<unsigned int, Hit2DSet>;
using PlaneToWireToHitSetMap       = std::map<geo::PlaneID, WireToHitSetMap>;
using TPCToPlaneToWireToHitSetMap  = std::map<geo::TPCID, PlaneToWireToHitSetMap>;
using HitVectorMap                 = std::map<size_t, HitVector>;
using SnippetHitMapItrPair         = std::pair<SnippetHitMap::iterator,SnippetHitMap::iterator>;
using PlaneSnippetHitMapItrPairVec = std::vector<SnippetHitMapItrPair>;

/**
 *  @brief  SnippetHit3DBuilderSBN class definiton
 */
class SnippetHit3DBuilderSBN : virtual public IHit3DBuilder
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit SnippetHit3DBuilderSBN(fhicl::ParameterSet const &pset);

    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::ProducesCollector&) override;

    void configure(const fhicl::ParameterSet&);

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param hitPairList           The input list of 3D hits to run clustering on
     *  @param clusterParametersList A list of cluster objects (parameters from associated hits)
     */
    virtual void Hit3DBuilder(art::Event&, reco::HitPairList&, RecobHitToPtrMap&) override;

    /**
     *  @brief If monitoring, recover the time to execute a particular function
     */
    virtual float getTimeToExecute(IHit3DBuilder::TimeValues index) const override {return m_timeVector[index];}

private:

    /**
     *  @brief  Extract the ART hits and the ART hit-particle relationships
     *
     *  @param  evt - the ART event
     */
    void CollectArtHits(const art::Event& evt) const;

    /**
     *  @brief Given the ClusterHit2D objects, build the HitPairMap
     */
    void BuildHit3D(reco::HitPairList& hitPairList) const;

    /**
     *  @brief Create a new 2D hit collection from hits associated to 3D space points
     */
    void CreateNewRecobHitCollection(art::Event&, reco::HitPairList&, std::vector<recob::Hit>&, RecobHitToPtrMap&);

    /**
     *  @brief Create recob::Wire to recob::Hit associations
     */
    void makeWireAssns(const art::Event&, art::Assns<recob::Wire, recob::Hit>&, RecobHitToPtrMap&) const;

    /**
     *  @brief Create raw::RawDigit to recob::Hit associations
     */
    void makeRawDigitAssns(const art::Event&, art::Assns<raw::RawDigit, recob::Hit>&, RecobHitToPtrMap&) const;

    /**
     *  @brief Given the ClusterHit2D objects, build the HitPairMap
     */
    size_t BuildHitPairMap(PlaneToSnippetHitMap& planeToHitVectorMap, reco::HitPairList& hitPairList) const;

    /**
     * @brief This defines a structure allowing us to see which hits have already been matched. Generally, if two
     *        hits have been matched then there is no point considering them anymore as it simply means duplicate
     *        space points will be created
     */
    using MatchedHitMap = std::unordered_map<const reco::ClusterHit2D*,std::set<const reco::ClusterHit2D*>>;

    /**
     *  @brief Given the ClusterHit2D objects, build the HitPairMap
     */
    size_t BuildHitPairMapByTPC(PlaneSnippetHitMapItrPairVec& planeSnippetHitMapItrPairVec, reco::HitPairList& hitPairList) const;

    /**
     *  @brief This builds a list of candidate hit pairs from lists of hits on two planes
     */
    using HitMatchTriplet       = std::tuple<const reco::ClusterHit2D*,const reco::ClusterHit2D*,const reco::ClusterHit3D>;
    using HitMatchTripletVec    = std::vector<HitMatchTriplet>;
    using HitMatchTripletVecMap = std::map<geo::WireID,HitMatchTripletVec>;

    int findGoodHitPairs(SnippetHitMap::iterator&, SnippetHitMap::iterator&, SnippetHitMap::iterator&, HitMatchTripletVecMap&) const;

    /**
     *  @brief This algorithm takes lists of hit pairs and finds good triplets
     */
    void findGoodTriplets(HitMatchTripletVecMap&, HitMatchTripletVecMap&, MatchedHitMap&, reco::HitPairList&, bool = false) const;

    /**
     * @brief This will look at storing pair "orphans" where the 2D hits are otherwise unused
     */

    int saveOrphanPairs(HitMatchTripletVecMap&, MatchedHitMap&, reco::HitPairList&) const;

    /**
     *  @brief Make a HitPair object by checking two hits
     */
    bool makeHitPair(reco::ClusterHit3D&       pairOut,
                     const reco::ClusterHit2D* hit1,
                     const reco::ClusterHit2D* hit2,
                     float                     hitWidthSclFctr = 1.,
                     size_t                    hitPairCntr = 0) const;

    /**
     *  @brief Make a 3D HitPair object by checking two hits
     */
    bool makeHitTriplet(reco::ClusterHit3D&       pairOut,
                        MatchedHitMap&            matchedHitMap,
                        const reco::ClusterHit3D& pairIn,
                        const reco::ClusterHit2D* hit2) const;

    /**
     *  @brief Make a 3D HitPair object from a valid pair and a dead channel in the missing plane
     */
    bool makeDeadChannelPair(reco::ClusterHit3D& pairOut, const reco::ClusterHit3D& pair, MatchedHitMap&, size_t maxStatus = 4, size_t minStatus = 0, float minOverlap=0.2) const;

    /**
     * @brief function to detemine if two wires "intersect" (in the 2D sense)
     */

    bool WireIDsIntersect(const geo::WireID&, const geo::WireID&, geo::WireIDIntersection&) const;

    /**
     * @brief function to compute the distance of closest approach and the arc length to the points of closest approach
     */
    float closestApproach(const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&, float&, float&) const;

    /**
     *  @brief A utility routine for finding a 2D hit closest in time to the given pair
     */
    const reco::ClusterHit2D* FindBestMatchingHit(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float pairDeltaTimeLimits) const;

    /**
     *  @brief A utility routine for returning the number of 2D hits from the list in a given range
     */
    int FindNumberInRange(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float range) const;

    /**
     *  @brief Jacket the calls to finding the nearest wire in order to intercept the exceptions if out of range
     */
    geo::WireID NearestWireID(const Eigen::Vector3f& position, const geo::WireID& wireID) const;

    /**
     *  @brief Jacket the calls to finding the nearest wire in order to intercept the exceptions if out of range
     */
    float DistanceFromPointToHitWire(const Eigen::Vector3f& position, const geo::WireID& wireID) const;

    /**
     *  @brief Create the internal channel status vector (assume will eventually be event-by-event)
     */
    void BuildChannelStatusVec(PlaneToWireToHitSetMap& planeToWiretoHitSetMap) const;

    /**
     * @brief Perform charge integration between limits
     */
    float chargeIntegral(float,float,float,float,int,int) const;

    /**
     *  @brief define data structure for keeping track of channel status
     */
    using ChannelStatusVec        = std::vector<size_t>;
    using ChannelStatusByPlaneVec = std::vector<ChannelStatusVec>;

    /**
     *  @brief clear the tuple vectors before processing next event
     */
    void clear();

    /**
     *  @brief Data members to follow
     */
    using TickCorrectionArray = std::vector<std::vector<std::vector<float>>>;

    std::vector<art::InputTag>              m_hitFinderTagVec;
    float                                   m_hitWidthSclFctr;
    float                                   m_deltaPeakTimeSig;
    float                                   m_rangeNumSig;
    float                                   m_LongHitStretchFctr;
    float                                   m_pulseHeightFrac;
    float                                   m_PHLowSelection;
    float                                   m_minPHFor2HitPoints;    ///< Set a minimum pulse height for 2 hit space point candidates
    std::vector<int>                        m_invalidTPCVec;
    float                                   m_wirePitchScaleFactor;  ///< Scaling factor to determine max distance allowed between candidate pairs
    float                                   m_maxHit3DChiSquare;     ///< Provide ability to select hits based on "chi square"
    bool                                    m_saveMythicalPoints;    ///< Should we save valid 2 hit space points? 
    float                                   m_maxMythicalChiSquare;  ///< Selection cut on mythical points
    bool                                    m_useT0Offsets;          ///< If true then we will use the LArSoft interplane offsets
    bool                                    m_outputHistograms;      ///< Take the time to create and fill some histograms for diagnostics
    bool                                    m_makeAssociations;      ///< Do we make wire/rawdigit associations to space points?
   
    bool                                    m_enableMonitoring;      ///<
    float                                   m_wirePitch[3];
    mutable std::vector<float>              m_timeVector;            ///<
   
    float                                   m_zPosOffset;
   
    using PlaneToT0OffsetMap = std::map<geo::PlaneID,float>;

    PlaneToT0OffsetMap                      m_PlaneToT0OffsetMap;
   
    // Define some basic histograms   
    TTree*                                  m_tupleTree;             ///< output analysis tree

    mutable std::vector<float>              m_deltaPeakTimePlane0Vec;
    mutable std::vector<float>              m_deltaPeakSigmaPlane0Vec;
    mutable std::vector<float>              m_deltaPeakTimePlane1Vec;
    mutable std::vector<float>              m_deltaPeakSigmaPlane1Vec;
    mutable std::vector<float>              m_deltaPeakTimePlane2Vec;
    mutable std::vector<float>              m_deltaPeakSigmaPlane2Vec;

    mutable std::vector<float>              m_deltaTimeVec;
    mutable std::vector<float>              m_deltaTime0Vec;
    mutable std::vector<float>              m_deltaSigma0Vec;
    mutable std::vector<float>              m_deltaTime1Vec;
    mutable std::vector<float>              m_deltaSigma1Vec;
    mutable std::vector<float>              m_deltaTime2Vec;
    mutable std::vector<float>              m_deltaSigma2Vec;
    mutable std::vector<float>              m_chiSquare3DVec;
    mutable std::vector<float>              m_maxPullVec;
    mutable std::vector<float>              m_overlapFractionVec;
    mutable std::vector<float>              m_overlapRangeVec;
    mutable std::vector<float>              m_maxDeltaPeakVec;
    mutable std::vector<float>              m_maxSideVecVec;
    mutable std::vector<float>              m_pairWireDistVec;
    mutable std::vector<float>              m_smallChargeDiffVec;
    mutable std::vector<int>                m_smallIndexVec;
    mutable std::vector<float>              m_qualityMetricVec;
    mutable std::vector<float>              m_spacePointChargeVec;
    mutable std::vector<float>              m_hitAsymmetryVec;
    mutable std::vector<float>              m_2hit1stPHVec;
    mutable std::vector<float>              m_2hit2ndPHVec;
    mutable std::vector<float>              m_2hitDeltaPHVec;
    mutable std::vector<float>              m_2hitSumPHVec;
   
    // Get instances of the primary data structures needed
    mutable Hit2DList                       m_clusterHit2DMasterList;
    mutable PlaneToSnippetHitMap            m_planeToSnippetHitMap;
    mutable PlaneToWireToHitSetMap          m_planeToWireToHitSetMap;


    mutable ChannelStatusByPlaneVec         m_channelStatus;
    mutable size_t                          m_numBadChannels;

    mutable bool                            m_weHaveAllBeenHereBefore = false;

    const geo::GeometryCore*                m_geometry;
    const geo::WireReadoutGeom*             m_wireReadout;
    const lariov::ChannelStatusProvider*    m_channelFilter;
};

SnippetHit3DBuilderSBN::SnippetHit3DBuilderSBN(fhicl::ParameterSet const &pset) :
    m_channelFilter(&art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider())
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SnippetHit3DBuilderSBN::produces(art::ProducesCollector& collector)
{
    collector.produces< std::vector<recob::Hit>>();
    collector.produces< art::Assns<recob::Wire,   recob::Hit>>();
    collector.produces< art::Assns<raw::RawDigit, recob::Hit>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SnippetHit3DBuilderSBN::configure(fhicl::ParameterSet const &pset)
{
    m_hitFinderTagVec      = pset.get<std::vector<art::InputTag>>("HitFinderTagVec",        {"gaushit"});
    m_enableMonitoring     = pset.get<bool                      >("EnableMonitoring",       true);
    m_hitWidthSclFctr      = pset.get<float                     >("HitWidthScaleFactor",    6.  );
    m_rangeNumSig          = pset.get<float                     >("RangeNumSigma",          3.  );
    m_LongHitStretchFctr   = pset.get<float                     >("LongHitsStretchFactor",  1.5 );
    m_pulseHeightFrac      = pset.get<float                     >("PulseHeightFraction",    0.5 );
    m_PHLowSelection       = pset.get<float                     >("PHLowSelection",         20. );
    m_minPHFor2HitPoints   = pset.get<float                     >("MinPHFor2HitPoints",     15. );
    m_deltaPeakTimeSig     = pset.get<float                     >("DeltaPeakTimeSig",       1.7 );
    m_zPosOffset           = pset.get<float                     >("ZPosOffset",             0.0 );
    m_invalidTPCVec        = pset.get<std::vector<int>          >("InvalidTPCVec",          std::vector<int>());
    m_wirePitchScaleFactor = pset.get<float                     >("WirePitchScaleFactor",   1.9 );
    m_maxHit3DChiSquare    = pset.get<float                     >("MaxHitChiSquare",        6.0 );
    m_saveMythicalPoints   = pset.get<bool                      >("SaveMythicalPoints",     true);
    m_maxMythicalChiSquare = pset.get<float                     >("MaxMythicalChiSquare",    10.);
    m_useT0Offsets         = pset.get<bool                      >("UseT0Offsets",           true);
    m_outputHistograms     = pset.get<bool                      >("OutputHistograms",      false);
    m_makeAssociations     = pset.get<bool                      >("MakeAssociations",      false);

    m_wireReadout = &art::ServiceHandle<geo::WireReadout const>{}->Get();

    // Returns the wire pitch per plane assuming they will be the same for all TPCs
    constexpr geo::TPCID tpcid{0, 0};
    m_wirePitch[0] = m_wireReadout->Plane({tpcid, 0}).WirePitch();
    m_wirePitch[1] = m_wireReadout->Plane({tpcid, 1}).WirePitch();
    m_wirePitch[2] = m_wireReadout->Plane({tpcid, 2}).WirePitch();

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    if (m_outputHistograms)
    {
        art::ServiceHandle<art::TFileService> tfs;

        m_tupleTree = tfs->make<TTree>("Hit3DBuilderTree", "Tree by SnippetHit3DBuilderSBN");

        clear();

        m_tupleTree->Branch("DeltaPeakTimePair0",  "std::vector<float>", &m_deltaPeakTimePlane0Vec);
        m_tupleTree->Branch("DeltaPeakSigmaPair0", "std::vector<float>", &m_deltaPeakSigmaPlane0Vec);
        m_tupleTree->Branch("DeltaPeakTimePair1",  "std::vector<float>", &m_deltaPeakTimePlane1Vec);
        m_tupleTree->Branch("DeltaPeakSigmaPair1", "std::vector<float>", &m_deltaPeakSigmaPlane1Vec);
        m_tupleTree->Branch("DeltaPeakTimePair2",  "std::vector<float>", &m_deltaPeakTimePlane2Vec);
        m_tupleTree->Branch("DeltaPeakSigmaPair2", "std::vector<float>", &m_deltaPeakSigmaPlane2Vec);

        m_tupleTree->Branch("DeltaTime2D",     "std::vector<float>", &m_deltaTimeVec);
        m_tupleTree->Branch("DeltaTime2D0",    "std::vector<float>", &m_deltaTime0Vec);
        m_tupleTree->Branch("DeltaSigma2D0",   "std::vector<float>", &m_deltaSigma0Vec);
        m_tupleTree->Branch("DeltaTime2D1",    "std::vector<float>", &m_deltaTime1Vec);
        m_tupleTree->Branch("DeltaSigma2D1",   "std::vector<float>", &m_deltaSigma1Vec);
        m_tupleTree->Branch("DeltaTime2D2",    "std::vector<float>", &m_deltaTime2Vec);
        m_tupleTree->Branch("DeltaSigma2D2",   "std::vector<float>", &m_deltaSigma2Vec);
        m_tupleTree->Branch("ChiSquare3D",     "std::vector<float>", &m_chiSquare3DVec);
        m_tupleTree->Branch("MaxPullValue",    "std::vector<float>", &m_maxPullVec);
        m_tupleTree->Branch("OverlapFraction", "std::vector<float>", &m_overlapFractionVec);
        m_tupleTree->Branch("OverlapRange",    "std::vector<float>", &m_overlapRangeVec);
        m_tupleTree->Branch("MaxDeltaPeak",    "std::vector<float>", &m_maxDeltaPeakVec);
        m_tupleTree->Branch("MaxSideVec",      "std::vector<float>", &m_maxSideVecVec);
        m_tupleTree->Branch("PairWireDistVec", "std::vector<float>", &m_pairWireDistVec);
        m_tupleTree->Branch("SmallChargeDiff", "std::vector<float>", &m_smallChargeDiffVec);
        m_tupleTree->Branch("SmallChargeIdx",  "std::vector<int>",   &m_smallIndexVec);
        m_tupleTree->Branch("QualityMetric",   "std::vector<float>", &m_qualityMetricVec);
        m_tupleTree->Branch("SPCharge",        "std::vector<float>", &m_spacePointChargeVec);
        m_tupleTree->Branch("HitAsymmetry",    "std::vector<float>", &m_hitAsymmetryVec);

        m_tupleTree->Branch("2hit1stPH",       "std::vector<float>", &m_2hit1stPHVec);
        m_tupleTree->Branch("2hit2ndPH",       "std::vector<float>", &m_2hit2ndPHVec);
        m_tupleTree->Branch("2hitDeltaPH",     "std::vector<float>", &m_2hitDeltaPHVec);
        m_tupleTree->Branch("2hitSumPH",       "std::vector<float>", &m_2hitSumPHVec);

    }

    return;
}

void SnippetHit3DBuilderSBN::clear()
{
    m_deltaPeakTimePlane0Vec.clear();
    m_deltaPeakSigmaPlane0Vec.clear();
    m_deltaPeakTimePlane1Vec.clear();
    m_deltaPeakSigmaPlane1Vec.clear();
    m_deltaPeakTimePlane1Vec.clear();
    m_deltaPeakSigmaPlane1Vec.clear();

    m_deltaTimeVec.clear();
    m_deltaTime0Vec.clear();
    m_deltaSigma0Vec.clear();
    m_deltaTime1Vec.clear();
    m_deltaSigma1Vec.clear();
    m_deltaTime2Vec.clear();
    m_deltaSigma2Vec.clear();
    m_chiSquare3DVec.clear();
    m_maxPullVec.clear();
    m_overlapFractionVec.clear();
    m_overlapRangeVec.clear();
    m_maxDeltaPeakVec.clear();
    m_maxSideVecVec.clear();
    m_pairWireDistVec.clear();
    m_smallChargeDiffVec.clear();
    m_smallIndexVec.clear();
    m_qualityMetricVec.clear();
    m_spacePointChargeVec.clear();
    m_hitAsymmetryVec.clear();

    m_2hit1stPHVec.clear();
    m_2hit2ndPHVec.clear();
    m_2hitDeltaPHVec.clear();
    m_2hitSumPHVec.clear();

    return;
}

void SnippetHit3DBuilderSBN::BuildChannelStatusVec(PlaneToWireToHitSetMap& planeToWireToHitSetMap) const
{
    // This is called each event, clear out the previous version and start over
    m_channelStatus.clear();

    m_numBadChannels = 0;
    m_channelStatus.resize(m_wireReadout->Nplanes());

    // Loop through views/planes to set the wire length vectors
    constexpr geo::TPCID tpcid{0, 0};
    for(size_t idx = 0; idx < m_channelStatus.size(); idx++)
    {
        m_channelStatus[idx] = ChannelStatusVec(m_wireReadout->Nwires(geo::PlaneID(tpcid, idx)), 5);
    }

    // Loop through the channels and mark those that are "bad"
    for(size_t channel = 0; channel < m_wireReadout->Nchannels(); channel++)
    {
        try
        {  
            if( m_channelFilter->IsPresent(channel) && !m_channelFilter->IsGood(channel))
            {
                std::vector<geo::WireID>                wireIDVec = m_wireReadout->ChannelToWire(channel);
                geo::WireID                             wireID    = wireIDVec[0];
                lariov::ChannelStatusProvider::Status_t chanStat  = m_channelFilter->Status(channel);

                m_channelStatus[wireID.Plane][wireID.Wire] = chanStat;
                m_numBadChannels++;
            }
        }
        catch(...)
        {
            mf::LogDebug("SnippetHit3D") << "--> Channel: " << channel << " threw exception so we will skip" << std::endl;
        }
    }

    std::cout << "SnippetHit3D finds " << m_numBadChannels << " bad channels" << std::endl;

    return;
}


bool SetPeakHitPairIteratorOrder(const reco::HitPairList::iterator& left, const reco::HitPairList::iterator& right)
{
    return (*left).getAvePeakTime() < (*right).getAvePeakTime();
}

struct HitPairClusterOrder
{
    bool operator()(const reco::HitPairClusterMap::iterator& left, const reco::HitPairClusterMap::iterator& right)
    {
        // Watch out for the case where two clusters can have the same number of hits!
        if (left->second.size() == right->second.size())
            return left->first < right->first;

        return left->second.size() > right->second.size();
    }
};

void SnippetHit3DBuilderSBN::Hit3DBuilder(art::Event& evt, reco::HitPairList& hitPairList, RecobHitToPtrMap& clusterHitToArtPtrMap)
{
    // Clear the internal data structures
    m_clusterHit2DMasterList.clear();
    m_planeToSnippetHitMap.clear();
    m_planeToWireToHitSetMap.clear();

    // Do the one time initialization of the tick offsets. 
    if (m_PlaneToT0OffsetMap.empty())
    {
        // Need the detector properties which needs the clocks 
        auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
        auto const det_prop   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);

        // Initialize the plane to hit vector map
        for(auto const& planeID : m_wireReadout->Iterate<geo::PlaneID>())
            {
          if (m_useT0Offsets) m_PlaneToT0OffsetMap[planeID] = det_prop.GetXTicksOffset(planeID) - det_prop.GetXTicksOffset(geo::PlaneID(planeID.parentID(),0));
                    else                m_PlaneToT0OffsetMap[planeID] = 0.;
                }
            }

    m_timeVector.resize(NUMTIMEVALUES, 0.);

    // Get a hit refiner
    std::unique_ptr<std::vector<recob::Hit>> outputHitPtrVec(new std::vector<recob::Hit>);

    // Recover the 2D hits and then organize them into data structures which will be used in the
    // DBscan algorithm for building the 3D clusters
    this->CollectArtHits(evt);

    // If there are no hits in our view/wire data structure then do not proceed with the full analysis
    if (!m_planeToWireToHitSetMap.empty())
    {
        // Call the algorithm that builds 3D hits
        this->BuildHit3D(hitPairList);

        // If we built 3D points then attempt to output a new hit list as well
        if (!hitPairList.empty())
            CreateNewRecobHitCollection(evt, hitPairList, *outputHitPtrVec, clusterHitToArtPtrMap);
    }

    // Set up to make the associations (if desired)
    /// Associations with wires.
    std::unique_ptr<art::Assns<recob::Wire, recob::Hit>> wireAssns(new art::Assns<recob::Wire, recob::Hit>);

    if (m_makeAssociations) makeWireAssns(evt, *wireAssns, clusterHitToArtPtrMap);

    /// Associations with raw digits.
    std::unique_ptr<art::Assns<raw::RawDigit, recob::Hit>> rawDigitAssns(new art::Assns<raw::RawDigit, recob::Hit>);

    if (m_makeAssociations) makeRawDigitAssns(evt, *rawDigitAssns, clusterHitToArtPtrMap);

    // Move everything into the event
    evt.put(std::move(outputHitPtrVec));
    evt.put(std::move(wireAssns));
    evt.put(std::move(rawDigitAssns));

    // Handle tree output too
    if (m_outputHistograms)
    {
        m_tupleTree->Fill();

        clear();
    }

    return;
}

void SnippetHit3DBuilderSBN::BuildHit3D(reco::HitPairList& hitPairList) const
{
    /**
     *  @brief Driver for processing input 2D hits, transforming to 3D hits and building lists
     *         of associated 3D hits (candidate 3D clusters)
     */
    cet::cpu_timer theClockMakeHits;

    if (m_enableMonitoring) theClockMakeHits.start();

    // The first task is to take the lists of input 2D hits (a map of view to sorted lists of 2D hits)
    // and then to build a list of 3D hits to be used in downstream processing
    BuildChannelStatusVec(m_planeToWireToHitSetMap);

    size_t numHitPairs = BuildHitPairMap(m_planeToSnippetHitMap, hitPairList);

    if (m_enableMonitoring)
    {
        theClockMakeHits.stop();

        m_timeVector[BUILDTHREEDHITS] = theClockMakeHits.accumulated_real_time();
    }

    mf::LogDebug("SnippetHit3D") << ">>>>> 3D hit building done, found " << numHitPairs << " 3D Hits" << std::endl;

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------
class SetStartTimeOrder
{
public:
    SetStartTimeOrder() {}

    bool operator()(const SnippetHitMapItrPair& left, const SnippetHitMapItrPair& right) const
    {
        // Special case handling, there is nothing to compare for the left or right
        if (left.first  == left.second)  return false;
        if (right.first == right.second) return true;

        // de-referencing a bunch of pairs here...
        return left.first->first.first < right.first->first.first;
    }

private:
};

bool SetPairStartTimeOrder(const reco::ClusterHit3D& left, const reco::ClusterHit3D& right)
{
    // Sort by "modified start time" of pulse
    return left.getAvePeakTime() - left.getSigmaPeakTime() < right.getAvePeakTime() - right.getSigmaPeakTime();
}


//------------------------------------------------------------------------------------------------------------------------------------------

size_t SnippetHit3DBuilderSBN::BuildHitPairMap(PlaneToSnippetHitMap& planeToSnippetHitMap, reco::HitPairList& hitPairList) const
{
    /**
     *  @brief Given input 2D hits, build out the lists of possible 3D hits
     *
     *         The current strategy: ideally all 3D hits would be comprised of a triplet of 2D hits, one from each view
     *         However, we have concern that, in particular, the v-plane may have some inefficiency which we have to be
     *         be prepared to deal with. The idea, then, is to first make the association of hits in the U and W planes
     *         and then look for the match in the V plane. In the event we don't find the match in the V plane then we
     *         will evaluate the situation and in some instances keep the U-W pairs in order to keep efficiency high.
     */
    size_t totalNumHits(0);
    size_t hitPairCntr(0);

    size_t nTriplets(0);

    // Set up to loop over cryostats and tpcs...
    for(auto const& tpcID : m_geometry->Iterate<geo::TPCID>())
        {
            //************************************
            // Kludge
//            if (!(cryoIdx == 1 && tpcIdx == 0)) continue;

      PlaneToSnippetHitMap::iterator mapItr0 = planeToSnippetHitMap.find(geo::PlaneID(tpcID,0));
      PlaneToSnippetHitMap::iterator mapItr1 = planeToSnippetHitMap.find(geo::PlaneID(tpcID,1));
      PlaneToSnippetHitMap::iterator mapItr2 = planeToSnippetHitMap.find(geo::PlaneID(tpcID,2));

            size_t nPlanesWithHits = (mapItr0 != planeToSnippetHitMap.end() && !mapItr0->second.empty() ? 1 : 0)
                                   + (mapItr1 != planeToSnippetHitMap.end() && !mapItr1->second.empty() ? 1 : 0)
                                   + (mapItr2 != planeToSnippetHitMap.end() && !mapItr2->second.empty() ? 1 : 0);

            if (nPlanesWithHits < 2) continue;

            SnippetHitMap& snippetHitMap0 = mapItr0->second;
            SnippetHitMap& snippetHitMap1 = mapItr1->second;
            SnippetHitMap& snippetHitMap2 = mapItr2->second;

            PlaneSnippetHitMapItrPairVec hitItrVec = {SnippetHitMapItrPair(snippetHitMap0.begin(),snippetHitMap0.end()),
                                                      SnippetHitMapItrPair(snippetHitMap1.begin(),snippetHitMap1.end()),
                                                      SnippetHitMapItrPair(snippetHitMap2.begin(),snippetHitMap2.end())};

            totalNumHits += BuildHitPairMapByTPC(hitItrVec, hitPairList);
    }

    // Return the hit pair list but sorted by z and y positions (faster traversal in next steps)
    hitPairList.sort(SetPairStartTimeOrder);

    // Where are we?
    mf::LogDebug("SnippetHit3D") << "Total number hits: " << totalNumHits << std::endl;
    mf::LogDebug("SnippetHit3D") << "Created a total of " << hitPairList.size() << " hit pairs, counted: " << hitPairCntr << std::endl;
    mf::LogDebug("SnippetHit3D") << "-- Triplets: " << nTriplets << std::endl;

    return hitPairList.size();
}

size_t SnippetHit3DBuilderSBN::BuildHitPairMapByTPC(PlaneSnippetHitMapItrPairVec& snippetHitMapItrVec, reco::HitPairList& hitPairList) const
{
    /**
     *  @brief Given input 2D hits, build out the lists of possible 3D hits
     *
     *         The current strategy: ideally all 3D hits would be comprised of a triplet of 2D hits, one from each view
     *         However, we have concern that, in particular, the v-plane may have some inefficiency which we have to be
     *         be prepared to deal with. The idea, then, is to first make the association of hits in the U and W planes
     *         and then look for the match in the V plane. In the event we don't find the match in the V plane then we
     *         will evaluate the situation and in some instances keep the U-W pairs in order to keep efficiency high.
     */

    // Define functions to set start/end iterators in the loop below
    auto SetStartIterator = [](SnippetHitMap::iterator startItr, SnippetHitMap::iterator endItr, float startTime)
    {
        while(startItr != endItr)
        {
            if (startItr->first.second < startTime) startItr++;
            else break;
        }
        return startItr;
    };

    auto SetEndIterator = [](SnippetHitMap::iterator lastItr, SnippetHitMap::iterator endItr, float endTime)
    {
        while(lastItr != endItr)
        {
            if (lastItr->first.first < endTime) lastItr++;
            else break;
        }
        return lastItr;
    };

    size_t nTriplets(0);
    size_t nOrphanPairs(0);

    // Structure to keep track of hit associations
    MatchedHitMap matchedHitMap;

    //*********************************************************************************
    // Basically, we try to loop until done...
    while(1)
    {
        // Sort so that the earliest hit time will be the first element, etc.
        std::sort(snippetHitMapItrVec.begin(),snippetHitMapItrVec.end(),SetStartTimeOrder());

        // Make sure there are still hits on at least
        int nPlanesWithHits(0);

        for(auto& pair : snippetHitMapItrVec)
            if (pair.first != pair.second) nPlanesWithHits++;

        if (nPlanesWithHits < 2) break;

        // End condition: no more hit snippets
//        if (snippetHitMapItrVec.front().first == snippetHitMapItrVec.front().second) break;

        // This loop iteration's snippet iterator
        SnippetHitMap::iterator firstSnippetItr = snippetHitMapItrVec.front().first;

        // Set iterators to insure we'll be in the overlap ranges
        SnippetHitMap::iterator snippetHitMapItr1Start = SetStartIterator(snippetHitMapItrVec[1].first, snippetHitMapItrVec[1].second, firstSnippetItr->first.first);
        SnippetHitMap::iterator snippetHitMapItr1End   = SetEndIterator( snippetHitMapItr1Start,        snippetHitMapItrVec[1].second, firstSnippetItr->first.second);
        SnippetHitMap::iterator snippetHitMapItr2Start = SetStartIterator(snippetHitMapItrVec[2].first, snippetHitMapItrVec[2].second, firstSnippetItr->first.first);
        SnippetHitMap::iterator snippetHitMapItr2End   = SetEndIterator( snippetHitMapItr2Start,        snippetHitMapItrVec[2].second, firstSnippetItr->first.second);

        // Since we'll use these many times in the internal loops, pre make the pairs for the second set of hits
        size_t                curHitListSize(hitPairList.size());
        HitMatchTripletVecMap pair12Map;
        HitMatchTripletVecMap pair13Map;

        size_t n12Pairs = findGoodHitPairs(firstSnippetItr, snippetHitMapItr1Start, snippetHitMapItr1End, pair12Map);
        size_t n13Pairs = findGoodHitPairs(firstSnippetItr, snippetHitMapItr2Start, snippetHitMapItr2End, pair13Map);

        curHitListSize  = hitPairList.size();

        if (n12Pairs > n13Pairs) findGoodTriplets(pair12Map, pair13Map, matchedHitMap, hitPairList);
        else                     findGoodTriplets(pair13Map, pair12Map, matchedHitMap, hitPairList);

        if (m_saveMythicalPoints)
        {
            nOrphanPairs += saveOrphanPairs(pair12Map, matchedHitMap, hitPairList);
            nOrphanPairs += saveOrphanPairs(pair13Map, matchedHitMap, hitPairList);
        }

        nTriplets += hitPairList.size() - curHitListSize;

        snippetHitMapItrVec.front().first++;
    }

    mf::LogDebug("SnippetHit3D") << "--> Created " << nTriplets << " triplets of which " << nOrphanPairs << " are orphans" << std::endl;

    return hitPairList.size();
}

int SnippetHit3DBuilderSBN::findGoodHitPairs(SnippetHitMap::iterator& firstSnippetItr,
                                                SnippetHitMap::iterator& startItr,
                                                SnippetHitMap::iterator& endItr,
                                                HitMatchTripletVecMap&   hitMatchMap) const
{
    int numPairs(0);

    HitVector::iterator firstMaxItr = std::max_element(firstSnippetItr->second.begin(),firstSnippetItr->second.end(),[](const auto& left, const auto& right){return left->getHit()->PeakAmplitude() < right->getHit()->PeakAmplitude();});
    float               firstPHCut  = firstMaxItr != firstSnippetItr->second.end() ? m_pulseHeightFrac * (*firstMaxItr)->getHit()->PeakAmplitude() : 4096.;

    // Loop through the hits on the first snippet
    for(const auto& hit1 : firstSnippetItr->second)
    {
        // Let's focus on the largest hit in the chain
        if (hit1->getHit()->DegreesOfFreedom() > 1 && hit1->getHit()->PeakAmplitude() < firstPHCut && hit1->getHit()->PeakAmplitude() < m_PHLowSelection) continue;

        // Inside loop iterator
        SnippetHitMap::iterator secondHitItr = startItr;

        // Loop through the input secon hits and make pairs
        while(secondHitItr != endItr)
        {
            HitVector::iterator secondMaxItr = std::max_element(secondHitItr->second.begin(),secondHitItr->second.end(),[](const auto& left, const auto& right){return left->getHit()->PeakAmplitude() < right->getHit()->PeakAmplitude();});
            float               secondPHCut  = secondMaxItr != secondHitItr->second.end() ? m_pulseHeightFrac * (*secondMaxItr)->getHit()->PeakAmplitude() : 0.;

            for(const auto& hit2 : secondHitItr->second)
            {
                // Again, focus on the large hits
                if (hit2->getHit()->DegreesOfFreedom() > 1 && hit2->getHit()->PeakAmplitude() < secondPHCut && hit2->getHit()->PeakAmplitude() < m_PHLowSelection) continue;

                reco::ClusterHit3D  pair;

                // pair returned with a negative ave time is signal of failure
                if (!makeHitPair(pair, hit1, hit2, m_hitWidthSclFctr)) continue;

                std::vector<const recob::Hit*> recobHitVec = {nullptr,nullptr,nullptr};

                recobHitVec[hit1->WireID().Plane] = hit1->getHit();
                recobHitVec[hit2->WireID().Plane] = hit2->getHit();

                geo::WireID wireID = hit2->WireID();

                hitMatchMap[wireID].emplace_back(hit1,hit2,pair);

                numPairs++;
            }

            secondHitItr++;
        }
    }

    return numPairs;
}

void SnippetHit3DBuilderSBN::findGoodTriplets(HitMatchTripletVecMap& pair12Map, 
                                                 HitMatchTripletVecMap& pair13Map, 
                                                 MatchedHitMap&         matchedHitMap,
                                                 reco::HitPairList&     hitPairList, 
                                                 bool                   tagged) const
{
    // Build triplets from the two lists of hit pairs
    if (!pair12Map.empty())
    {
        // temporary container for dead channel hits
        std::vector<reco::ClusterHit3D> tempDeadChanVec;
        reco::ClusterHit3D              deadChanPair;

        // Keep track of which third plane hits have been used
        std::map<const reco::ClusterHit3D*,bool> usedPairMap;

        // Initial population of this map with the pair13Map hits
        for(const auto& pair13 : pair13Map)
        {
            for(const auto& hit2Dhit3DPair : pair13.second) usedPairMap[&std::get<2>(hit2Dhit3DPair)] = false;
        }

        // The outer loop is over all hit pairs made from the first two plane combinations
        for(const auto& pair12 : pair12Map)
        {
            if (pair12.second.empty()) continue;

            // This loop is over hit pairs that share the same first two plane wires but may have different
            // hit times on those wires
            for(const auto& hit2Dhit3DPair12 : pair12.second)
            {
                const reco::ClusterHit3D& pair1  = std::get<2>(hit2Dhit3DPair12);

                // populate the map with initial value
                usedPairMap[&pair1] = false;

                // The simplest approach here is to loop over all possibilities and let the triplet builder weed out the weak candidates
                for(const auto& pair13 : pair13Map)
                {
                    if (pair13.second.empty()) continue;

                    for(const auto& hit2Dhit3DPair13 : pair13.second)
                    {
                        // Protect against double counting
                        if (std::get<0>(hit2Dhit3DPair12) != std::get<0>(hit2Dhit3DPair13)) continue;

                        const reco::ClusterHit2D* hit2  = std::get<1>(hit2Dhit3DPair13);
                        const reco::ClusterHit3D& pair2 = std::get<2>(hit2Dhit3DPair13);

                        // If success try for the triplet
                        reco::ClusterHit3D triplet;

                        if (makeHitTriplet(triplet, matchedHitMap, pair1, hit2))
                        {
                            triplet.setID(hitPairList.size());
                            hitPairList.emplace_back(triplet);
                            usedPairMap[&pair1] = true;
                            usedPairMap[&pair2] = true;
                        }
                    }
                }
            }
        }

        // One more loop through the other pairs to check for sick channels
        if (m_numBadChannels > 0)
        {
            for(const auto& pairMapPair : usedPairMap)
            {
                if (pairMapPair.second) continue;

                const reco::ClusterHit3D* pair = pairMapPair.first;

                // Here we look to see if we failed to make a triplet because the partner wire was dead/noisy/sick
                if (makeDeadChannelPair(deadChanPair, *pair, matchedHitMap, 4, 0, 0.)) tempDeadChanVec.emplace_back(deadChanPair);
            }

            // Handle the dead wire triplets
            if(!tempDeadChanVec.empty())
            {
                // If we have many then see if we can trim down a bit by keeping those with time significance
                if (tempDeadChanVec.size() > 1)
                {
                    // Sort by "significance" of agreement
                    std::sort(tempDeadChanVec.begin(),tempDeadChanVec.end(),[](const auto& left, const auto& right){return left.getDeltaPeakTime()/left.getSigmaPeakTime() < right.getDeltaPeakTime()/right.getSigmaPeakTime();});

                    // What is the range of "significance" from first to last?
                    float firstSig = tempDeadChanVec.front().getDeltaPeakTime() / tempDeadChanVec.front().getSigmaPeakTime();
                    float lastSig  = tempDeadChanVec.back().getDeltaPeakTime()  / tempDeadChanVec.back().getSigmaPeakTime();
                    float sigRange = lastSig - firstSig;

                    if (lastSig > 0.5 * m_deltaPeakTimeSig && sigRange > 0.5)
                    {
                        // Declare a maximum of 1.5 * the average of the first and last pairs...
                        float maxSignificance = std::max(0.75 * (firstSig + lastSig),1.0);

                        std::vector<reco::ClusterHit3D>::iterator firstBadElem = std::find_if(tempDeadChanVec.begin(),tempDeadChanVec.end(),[&maxSignificance](const auto& pair){return pair.getDeltaPeakTime()/pair.getSigmaPeakTime() > maxSignificance;});

                        // But only keep the best 10?
                        if (std::distance(tempDeadChanVec.begin(),firstBadElem) > 20) firstBadElem = tempDeadChanVec.begin() + 20;
                        // Keep at least one hit...
                        else if (firstBadElem == tempDeadChanVec.begin()) firstBadElem++;

                        tempDeadChanVec.resize(std::distance(tempDeadChanVec.begin(),firstBadElem));
                    }
                }

                for(auto& pair : tempDeadChanVec)
                {
                    pair.setID(hitPairList.size());
                    hitPairList.emplace_back(pair);
                }
            }
        }
    }

    return;
}

int SnippetHit3DBuilderSBN::saveOrphanPairs(HitMatchTripletVecMap& pairMap, MatchedHitMap& matchedHitMap, reco::HitPairList& hitPairList) const
{
    int curTripletCount = hitPairList.size();

    // Build triplets from the two lists of hit pairs
    if (!pairMap.empty())
    {
        // Initial population of this map with the pair13Map hits
        for(const auto& pair : pairMap)
        {
            if (pair.second.empty()) continue;

            // This loop is over hit pairs that share the same first two plane wires but may have different
            // hit times on those wires
            for(const auto& hit2Dhit3DPair : pair.second)
            {
                const reco::ClusterHit3D& hit3D = std::get<2>(hit2Dhit3DPair);

                // No point considering a 3D hit that has been used to make a space point already
                if (hit3D.getStatusBits() & reco::ClusterHit3D::MADESPACEPOINT) continue;

                const reco::ClusterHit2D* hit1 = std::get<0>(hit2Dhit3DPair);
                const reco::ClusterHit2D* hit2 = std::get<1>(hit2Dhit3DPair);

                // Have these already been used in some fashion?
                if (matchedHitMap[hit1].find(hit2) != matchedHitMap[hit1].end()) continue;

                if (m_outputHistograms)
                {
                    m_2hit1stPHVec.emplace_back(hit1->getHit()->PeakAmplitude());
                    m_2hit2ndPHVec.emplace_back(hit2->getHit()->PeakAmplitude());
                    m_2hitDeltaPHVec.emplace_back(hit2->getHit()->PeakAmplitude() - hit1->getHit()->PeakAmplitude());
                    m_2hitSumPHVec.emplace_back(hit2->getHit()->PeakAmplitude() + hit1->getHit()->PeakAmplitude());
                }

                // If both hits already appear in a triplet then there is no gain here so reject
                if (hit1->getHit()->PeakAmplitude() < m_minPHFor2HitPoints || hit2->getHit()->PeakAmplitude() < m_minPHFor2HitPoints) continue;

                // Require that one of the hits is on the collection plane
                if (hit1->WireID().Plane == 2 || hit2->WireID().Plane == 2)
                {
                    // Allow cut on the quality of the space point
                    if (hit3D.getHitChiSquare() < m_maxMythicalChiSquare)
                    {
                        // Add to the list
                        hitPairList.emplace_back(hit3D);
                        hitPairList.back().setID(hitPairList.size()-1);
                        matchedHitMap[hit1].insert(hit2);
                        matchedHitMap[hit2].insert(hit1);
                    }
                }
            }
        }
    }

    return hitPairList.size() - curTripletCount;
}

bool SnippetHit3DBuilderSBN::makeHitPair(reco::ClusterHit3D&       hitPair,
                                            const reco::ClusterHit2D* hit1,
                                            const reco::ClusterHit2D* hit2,
                                            float                     hitWidthSclFctr,
                                            size_t                    hitPairCntr) const
{
    // Assume failure
    bool result(false);

    // Start by checking time consistency since this is fastest
    // Wires intersect so now we can check the timing
    float hit1Peak  = hit1->getTimeTicks();
    float hit1Sigma = hit1->getHit()->RMS();

    float hit2Peak  = hit2->getTimeTicks();
    float hit2Sigma = hit2->getHit()->RMS();

    // "Long hits" are an issue... so we deal with these differently
    int   hit1NDF   = hit1->getHit()->DegreesOfFreedom();
    int   hit2NDF   = hit2->getHit()->DegreesOfFreedom();

    // Basically, allow the range to extend to the nearest end of the snippet
    if (hit1NDF < 2) hit1Sigma *= m_LongHitStretchFctr; // This sets the range to the width of the pulse
    if (hit2NDF < 2) hit2Sigma *= m_LongHitStretchFctr;

    // The "hit sigma" is the gaussian fit sigma of the hit, we need to expand a bit to allow hit overlap efficiency
    float hit1Width = hitWidthSclFctr * hit1Sigma;
    float hit2Width = hitWidthSclFctr * hit2Sigma;

    // Coarse check hit times are "in range"
    if (fabs(hit1Peak - hit2Peak) <= (hit1Width + hit2Width))
    {
        // Check to see that hit peak times are consistent with each other
        float hit1SigSq     = std::pow(hit1Sigma,2);
        float hit2SigSq     = std::pow(hit2Sigma,2);
        float deltaPeakTime = std::fabs(hit1Peak - hit2Peak);
        float sigmaPeakTime = std::sqrt(hit1SigSq + hit2SigSq);

        if (m_outputHistograms)
        {
            // brute force... sigh... 
            int plane1   = hit1->WireID().Plane;
            int plane2   = hit2->WireID().Plane;
            int planeIdx = (plane1 + plane2) - 1;     // should be 0 for 0-1, 1 for 0-2 and 2 for 1-2

            float deltaTicks = hit2Peak - hit1Peak;

            if (plane1 > plane2) deltaTicks = -deltaTicks;

            if (planeIdx == 0)
            {
                m_deltaPeakTimePlane0Vec.emplace_back(deltaTicks);
                m_deltaPeakSigmaPlane0Vec.emplace_back(sigmaPeakTime);
            }
            else if (planeIdx == 1)
            {
                m_deltaPeakTimePlane1Vec.emplace_back(deltaTicks);
                m_deltaPeakSigmaPlane1Vec.emplace_back(sigmaPeakTime);
            }
            else
            {
                m_deltaPeakTimePlane2Vec.emplace_back(deltaTicks);
                m_deltaPeakSigmaPlane2Vec.emplace_back(sigmaPeakTime);
            }
        }

        // delta peak time consistency check here
        if (deltaPeakTime < m_deltaPeakTimeSig * sigmaPeakTime)    // 2 sigma consistency? (do this way to avoid divide)
        {
            // We assume in this routine that we are looking at hits in different views
            // The first mission is to check that the wires intersect
            const geo::WireID& hit1WireID = hit1->WireID();
            const geo::WireID& hit2WireID = hit2->WireID();

            geo::WireIDIntersection widIntersect;

            if (WireIDsIntersect(hit1WireID, hit2WireID, widIntersect))
            {
                float oneOverWghts  = hit1SigSq * hit2SigSq / (hit1SigSq + hit2SigSq);
                float avePeakTime   = (hit1Peak / hit1SigSq + hit2Peak / hit2SigSq) * oneOverWghts;
                float averageCharge = 0.5 * (hit1->getHit()->Integral() + hit2->getHit()->Integral());
                float hitChiSquare  = std::pow((hit1Peak - avePeakTime),2) / hit1SigSq
                                    + std::pow((hit2Peak - avePeakTime),2) / hit2SigSq;

                float xPositionHit1(hit1->getXPosition());
                float xPositionHit2(hit2->getXPosition());
                float xPosition = (xPositionHit1 / hit1SigSq + xPositionHit2 / hit2SigSq) * hit1SigSq * hit2SigSq / (hit1SigSq + hit2SigSq);

                Eigen::Vector3f position(xPosition, float(widIntersect.y), float(widIntersect.z)-m_zPosOffset);

                // If to here then we need to sort out the hit pair code telling what views are used
                unsigned statusBits = 1 << hit1->WireID().Plane | 1 << hit2->WireID().Plane;

                // handle status bits for the 2D hits
                if (hit1->getStatusBits() & reco::ClusterHit2D::USEDINPAIR) hit1->setStatusBit(reco::ClusterHit2D::SHAREDINPAIR);
                if (hit2->getStatusBits() & reco::ClusterHit2D::USEDINPAIR) hit2->setStatusBit(reco::ClusterHit2D::SHAREDINPAIR);

                hit1->setStatusBit(reco::ClusterHit2D::USEDINPAIR);
                hit2->setStatusBit(reco::ClusterHit2D::USEDINPAIR);

                reco::ClusterHit2DVec hitVector;

                hitVector.resize(3, NULL);

                hitVector[hit1->WireID().Plane] = hit1;
                hitVector[hit2->WireID().Plane] = hit2;

                unsigned int cryostatIdx = hit1->WireID().Cryostat;
                unsigned int tpcIdx      = hit1->WireID().TPC;

                // Initialize the wireIdVec
                std::vector<geo::WireID> wireIDVec = {geo::WireID(cryostatIdx,tpcIdx,0,0),
                                                      geo::WireID(cryostatIdx,tpcIdx,1,0),
                                                      geo::WireID(cryostatIdx,tpcIdx,2,0)};

                wireIDVec[hit1->WireID().Plane] = hit1->WireID();
                wireIDVec[hit2->WireID().Plane] = hit2->WireID();

                // For compiling at the moment
                std::vector<float> hitDelTSigVec = {0.,0.,0.};

                hitDelTSigVec[hit1->WireID().Plane] = deltaPeakTime / sigmaPeakTime;
                hitDelTSigVec[hit2->WireID().Plane] = deltaPeakTime / sigmaPeakTime;

                // Create the 3D cluster hit
                hitPair.initialize(hitPairCntr,
                                   statusBits,
                                   position,
                                   averageCharge,
                                   avePeakTime,
                                   deltaPeakTime,
                                   sigmaPeakTime,
                                   hitChiSquare,
                                   0.,
                                   0.,
                                   0.,
                                   0.,
                                   hitVector,
                                   hitDelTSigVec,
                                   wireIDVec);

                result = true;
            }
        }
    }

    // Send it back
    return result;
}


bool SnippetHit3DBuilderSBN::makeHitTriplet(reco::ClusterHit3D&       hitTriplet,
                                               MatchedHitMap&            matchedHitMap,
                                               const reco::ClusterHit3D& pair,
                                               const reco::ClusterHit2D* hit) const
{
    // Assume failure
    bool result(false);

    // We are going to force the wire pitch here, some time in the future we need to fix
    static const double wirePitch = 0.5 * m_wirePitchScaleFactor * *std::max_element(m_wirePitch,m_wirePitch+3);

    // Recover hit info
    float hitTimeTicks = hit->getTimeTicks();
    float hitSigma     = hit->getHit()->RMS();

    // Special case check
    if (hitSigma > 2. * hit->getHit()->PeakAmplitude()) hitSigma = 2. * hit->getHit()->PeakAmplitude();

    // Let's do a quick consistency check on the input hits to make sure we are in range...
    // Require the W hit to be "in range" with the UV Pair
    if (fabs(hitTimeTicks - pair.getAvePeakTime()) < m_hitWidthSclFctr * (pair.getSigmaPeakTime() + hitSigma))
    {
        // Check the distance from the point to the wire the hit is on
        float hitWireDist = DistanceFromPointToHitWire(pair.getPosition(), hit->WireID());

        if (m_outputHistograms) m_maxSideVecVec.push_back(hitWireDist);

        // Reject hits that are not within range
        if (hitWireDist < wirePitch)
        {
            if (m_outputHistograms) m_pairWireDistVec.push_back(hitWireDist);

            // Let's mark that this pair had a valid 3 wire combination
            pair.setStatusBit(reco::ClusterHit3D::MADESPACEPOINT);

            // Use the existing code to see the U and W hits are willing to pair with the V hit
            reco::ClusterHit3D pair0h;
            reco::ClusterHit3D pair1h;

            // Recover all the hits involved
            const reco::ClusterHit2DVec& pairHitVec = pair.getHits();
            const reco::ClusterHit2D*    hit0       = pairHitVec[0];
            const reco::ClusterHit2D*    hit1       = pairHitVec[1];

            if      (!hit0) hit0 = pairHitVec[2];
            else if (!hit1) hit1 = pairHitVec[2];

            bool notMatched =  matchedHitMap[hit].find(hit0)  == matchedHitMap[hit].end() 
                            && matchedHitMap[hit].find(hit1)  == matchedHitMap[hit].end()
                            && matchedHitMap[hit0].find(hit1) == matchedHitMap[hit0].end();

            // If good pairs made here then we can try to make a triplet
            if (notMatched && makeHitPair(pair0h, hit0, hit, m_hitWidthSclFctr) && makeHitPair(pair1h, hit1, hit, m_hitWidthSclFctr))
            {
                // Get a copy of the input hit vector (note the order is by plane - by definition)
                reco::ClusterHit2DVec hitVector = pair.getHits();

                // include the new hit
                hitVector[hit->WireID().Plane] = hit;

                // Set up to get average peak time, hitChiSquare, etc.
                unsigned int statusBits(0x7);
                float        avePeakTime(0.);
                float        weightSum(0.);
                float        xPosition(0.);

                // And get the wire IDs
                std::vector<geo::WireID> wireIDVec = {geo::WireID(), geo::WireID(), geo::WireID()};

                // First loop through the hits to get WireIDs and calculate the averages
                for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
                {
                    const reco::ClusterHit2D* hit2D = hitVector[planeIdx];

                    wireIDVec[planeIdx] = hit2D->WireID();

                    float hitRMS   = hit2D->getHit()->RMS();
                    float peakTime = hit2D->getTimeTicks();

                    // Basically, allow the range to extend to the nearest end of the snippet
                    if (hit2D->getHit()->DegreesOfFreedom() < 2) hitRMS *= m_LongHitStretchFctr;
                        //hitRMS = std::min(hit2D->getTimeTicks() - float(hit2D->getHit()->StartTick()),float(hit2D->getHit()->EndTick())-hit2D->getTimeTicks());

                    float weight = 1. / (hitRMS * hitRMS);

                    avePeakTime += peakTime * weight;
                    xPosition   += hit2D->getXPosition() * weight;
                    weightSum   += weight;
                }

                avePeakTime /= weightSum;
                xPosition   /= weightSum;

                Eigen::Vector2f pair0hYZVec(pair0h.getPosition()[1],pair0h.getPosition()[2]);
                Eigen::Vector2f pair1hYZVec(pair1h.getPosition()[1],pair1h.getPosition()[2]);
                Eigen::Vector2f pairYZVec(pair.getPosition()[1],pair.getPosition()[2]);
                Eigen::Vector3f position(xPosition,
                                         float((pairYZVec[0] + pair0hYZVec[0] + pair1hYZVec[0]) / 3.),
                                         float((pairYZVec[1] + pair0hYZVec[1] + pair1hYZVec[1]) / 3.));

                // Armed with the average peak time, now get hitChiSquare and the sig vec
                float              hitChiSquare(0.);
                float              sigmaPeakTime(std::sqrt(1./weightSum));
                std::vector<float> hitDelTSigVec;

                for(const auto& hit2D : hitVector)
                {
                    float hitRMS = hit2D->getHit()->RMS();

                    // Basically, allow the range to extend to the nearest end of the snippet
                    if (hit2D->getHit()->DegreesOfFreedom() < 2) hitRMS *= m_LongHitStretchFctr;
                        //hitRMS = std::min(hit2D->getTimeTicks() - float(hit2D->getHit()->StartTick()),float(hit2D->getHit()->EndTick())-hit2D->getTimeTicks());

                    float combRMS   = std::sqrt(hitRMS*hitRMS - sigmaPeakTime*sigmaPeakTime);
                    float peakTime  = hit2D->getTimeTicks();
                    float deltaTime = peakTime - avePeakTime;
                    float hitSig    = deltaTime / combRMS;

                    hitChiSquare += hitSig * hitSig;

                    hitDelTSigVec.emplace_back(std::fabs(hitSig));
                }

                if (m_outputHistograms) m_chiSquare3DVec.push_back(hitChiSquare);

                int lowMinIndex(std::numeric_limits<int>::max());
                int lowMaxIndex(std::numeric_limits<int>::min());
                int hiMinIndex(std::numeric_limits<int>::max());
                int hiMaxIndex(std::numeric_limits<int>::min());

                // First task is to get the min/max values for the common overlap region
                for(const auto& hit2D : hitVector)
                {
                    float range = m_rangeNumSig * hit2D->getHit()->RMS();

                    // Basically, allow the range to extend to the nearest end of the snippet
                    if (hit2D->getHit()->DegreesOfFreedom() < 2) range *= m_LongHitStretchFctr;
                        //range = std::min(hit2D->getTimeTicks() - float(hit2D->getHit()->StartTick()),float(hit2D->getHit()->EndTick())-hit2D->getTimeTicks());

                    int hitStart = hit2D->getHit()->PeakTime() - range - 0.5;
                    int hitStop  = hit2D->getHit()->PeakTime() + range + 0.5;

                    lowMinIndex = std::min(hitStart,    lowMinIndex);
                    lowMaxIndex = std::max(hitStart,    lowMaxIndex);
                    hiMinIndex  = std::min(hitStop + 1, hiMinIndex);
                    hiMaxIndex  = std::max(hitStop + 1, hiMaxIndex);
                }

                // Keep only "good" hits...
                if (hitChiSquare < m_maxHit3DChiSquare && hiMinIndex > lowMaxIndex)
                {
                    // One more pass through hits to get charge
                    std::vector<float> chargeVec;

                    for(const auto& hit2D : hitVector)
                        chargeVec.push_back(chargeIntegral(hit2D->getHit()->PeakTime(),hit2D->getHit()->PeakAmplitude(),hit2D->getHit()->RMS(),1.,lowMaxIndex,hiMinIndex));

                    float totalCharge     = std::accumulate(chargeVec.begin(),chargeVec.end(),0.) / float(chargeVec.size());
                    float overlapRange    = float(hiMinIndex - lowMaxIndex);
                    float overlapFraction = overlapRange / float(hiMaxIndex - lowMinIndex);

                    // Set up to compute the charge asymmetry
                    std::vector<float> smallestChargeDiffVec;
                    std::vector<float> chargeAveVec;
                    float              smallestDiff(std::numeric_limits<float>::max());
                    float              maxDeltaPeak(0.);
                    size_t             chargeIndex(0);

                    for(size_t idx = 0; idx < 3; idx++)
                    {
                        size_t leftIdx  = (idx + 2) % 3;
                        size_t rightIdx = (idx + 1) % 3;

                        smallestChargeDiffVec.push_back(std::abs(chargeVec[leftIdx] - chargeVec[rightIdx]));
                        chargeAveVec.push_back(float(0.5 * (chargeVec[leftIdx] + chargeVec[rightIdx])));

                        if (smallestChargeDiffVec.back() < smallestDiff)
                        {
                            smallestDiff = smallestChargeDiffVec.back();
                            chargeIndex  = idx;
                        }

                        // Take opportunity to look at peak time diff
                        float deltaPeakTime = hitVector[rightIdx]->getTimeTicks() - hitVector[leftIdx]->getTimeTicks();

                        if (std::abs(deltaPeakTime) > maxDeltaPeak) maxDeltaPeak = std::abs(deltaPeakTime);

                        if (m_outputHistograms) 
                        {
                            int   deltaTimeIdx = (hitVector[leftIdx]->WireID().Plane + hitVector[rightIdx]->WireID().Plane) - 1;
                            float combRMS      = std::sqrt(std::pow(hitVector[leftIdx]->getHit()->RMS(),2) + std::pow(hitVector[rightIdx]->getHit()->RMS(),2));

                            m_deltaTimeVec.push_back(deltaPeakTime);

                            // Want to get the sign of the difference correct...
                            if (deltaTimeIdx == 0)  // This is planes 1 and 0
                            {
                                m_deltaTime0Vec.emplace_back(float(hitVector[1]->getTimeTicks() - hitVector[0]->getTimeTicks()));
                                m_deltaSigma0Vec.emplace_back(combRMS);
                            }
                            else if (deltaTimeIdx == 1)  // This is planes 0 and 2
                            {
                                m_deltaTime1Vec.emplace_back(float(hitVector[2]->getTimeTicks() - hitVector[0]->getTimeTicks()));
                                m_deltaSigma1Vec.emplace_back(combRMS);
                            }
                            else // must be planes 1 and 2
                            {
                                m_deltaTime2Vec.emplace_back(float(hitVector[2]->getTimeTicks() - hitVector[1]->getTimeTicks()));
                                m_deltaSigma2Vec.emplace_back(combRMS);
                            }
                        }
                    }

                    float chargeAsymmetry = (chargeAveVec[chargeIndex] - chargeVec[chargeIndex]) / (chargeAveVec[chargeIndex] + chargeVec[chargeIndex]);

                    // If this is true there has to be a negative charge that snuck in somehow
                    if (chargeAsymmetry < -1. || chargeAsymmetry > 1.)
                    {
                        const geo::WireID& hitWireID = hitVector[chargeIndex]->WireID();

                        std::cout << "============> Charge asymmetry out of range: " << chargeAsymmetry << " <============" << std::endl;
                        std::cout << "     hit C: " << hitWireID.Cryostat << ", TPC: " << hitWireID.TPC << ", Plane: " << hitWireID.Plane << ", Wire: " << hitWireID.Wire << std::endl;
                        std::cout << "     charge: " << chargeVec[0] << ", " << chargeVec[1] << ", " << chargeVec[2] << std::endl;
                        std::cout << "     index: " << chargeIndex << ", smallest diff: " << smallestDiff << std::endl;
                        return result;
                    }

                    // Usurping "deltaPeakTime" to be the maximum pull
                    float deltaPeakTime = *std::max_element(hitDelTSigVec.begin(),hitDelTSigVec.end());

                    if (m_outputHistograms)
                    {
                        m_smallChargeDiffVec.push_back(smallestDiff);
                        m_smallIndexVec.push_back(chargeIndex);
                        m_maxPullVec.push_back(deltaPeakTime);
                        m_qualityMetricVec.push_back(hitChiSquare);
                        m_spacePointChargeVec.push_back(totalCharge);
                        m_overlapFractionVec.push_back(overlapFraction);
                        m_overlapRangeVec.push_back(overlapRange);
                        m_maxDeltaPeakVec.push_back(maxDeltaPeak);
                        m_hitAsymmetryVec.push_back(chargeAsymmetry);
                    }

                    // Try to weed out cases where overlap doesn't match peak separation
                    if (maxDeltaPeak > 3. * overlapRange) return result;

                    // Create the 3D cluster hit
                    hitTriplet.initialize(0,
                                          statusBits,
                                          position,
                                          totalCharge,
                                          avePeakTime,
                                          deltaPeakTime,
                                          sigmaPeakTime,
                                          hitChiSquare,
                                          overlapFraction,
                                          chargeAsymmetry,
                                          0.,
                                          0.,
                                          hitVector,
                                          hitDelTSigVec,
                                          wireIDVec);
                    
                    // Since we are keeping the triplet, mark the hits as used
                    for(size_t hit2DIdx = 0; hit2DIdx < hitVector.size(); hit2DIdx++)
                    {
                        const reco::ClusterHit2D* hit2D = hitVector[hit2DIdx];

                        matchedHitMap[hit2D].insert(hitVector[(hit2DIdx+1)%hitVector.size()]);
                        matchedHitMap[hit2D].insert(hitVector[(hit2DIdx+2)%hitVector.size()]);

                        if (hit2D->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit2D->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);

                        hit2D->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);
                    }

                    // Mark the input pair
                    pair.setStatusBit(reco::ClusterHit3D::MADESPACEPOINT);

                    result = true;
                }
//                else std::cout << "-Rejecting triple with chiSquare: " << hitChiSquare << " and hiMinIndex: " << hiMinIndex << ", loMaxIndex: " << lowMaxIndex << std::endl;
            }
        }
    }
//    else std::cout << "-MakeTriplet hit cut, delta: " << hitTimeTicks - pair.getAvePeakTime() << ", min scale fctr: " <<m_hitWidthSclFctr << ", pair sig: " << pair.getSigmaPeakTime() << ", hitSigma: " << hitSigma << std::endl;

    // return success/fail
    return result;
}

bool SnippetHit3DBuilderSBN::WireIDsIntersect(const geo::WireID& wireID0, const geo::WireID& wireID1, geo::WireIDIntersection& widIntersection) const
{
    bool success(false);

    // Do quick check that things are in the same logical TPC
    if (wireID0.Cryostat != wireID1.Cryostat || wireID0.TPC != wireID1.TPC || wireID0.Plane == wireID1.Plane) return success;
        
    // Recover wire geometry information for each wire
    const geo::WireGeo& wireGeo0 = m_wireReadout->Wire(wireID0);
    const geo::WireGeo& wireGeo1 = m_wireReadout->Wire(wireID1);

    // Get wire position and direction for first wire
    auto wirePosArr = wireGeo0.GetCenter();

    Eigen::Vector3f wirePos0(wirePosArr.X(),wirePosArr.Y(),wirePosArr.Z());
    Eigen::Vector3f wireDir0(wireGeo0.Direction().X(),wireGeo0.Direction().Y(),wireGeo0.Direction().Z());

    //*********************************
    // Kludge
//    if (wireID0.Plane > 0) wireDir0[2] = -wireDir0[2];

    // And now the second one
    wirePosArr = wireGeo1.GetCenter();

    Eigen::Vector3f wirePos1(wirePosArr.X(),wirePosArr.Y(),wirePosArr.Z());
    Eigen::Vector3f wireDir1(wireGeo1.Direction().X(),wireGeo1.Direction().Y(),wireGeo1.Direction().Z());

    //**********************************
    // Kludge
//    if (wireID1.Plane > 0) wireDir1[2] = -wireDir1[2];

    // Get the distance of closest approach
    float arcLen0;
    float arcLen1;

    if (closestApproach(wirePos0, wireDir0, wirePos1, wireDir1, arcLen0, arcLen1))
    {
        // Now check that arc lengths are within range
        if (std::abs(arcLen0) < wireGeo0.HalfL() && std::abs(arcLen1) < wireGeo1.HalfL())
        {
            Eigen::Vector3f poca0 = wirePos0 + arcLen0 * wireDir0;

            widIntersection.y = poca0[1];
            widIntersection.z = poca0[2];

            success = true;
        }
    }

    return success;
}

float SnippetHit3DBuilderSBN::closestApproach(const Eigen::Vector3f& P0,
                                                 const Eigen::Vector3f& u0,
                                                 const Eigen::Vector3f& P1,
                                                 const Eigen::Vector3f& u1,
                                                 float&                 arcLen0,
                                                 float&                 arcLen1) const
{
    // Technique is to compute the arclength to each point of closest approach
    Eigen::Vector3f w0 = P0 - P1;
    float a(1.);
    float b(u0.dot(u1));
    float c(1.);
    float d(u0.dot(w0));
    float e(u1.dot(w0));
    float den(a * c - b * b);

    arcLen0 = (b * e - c * d) / den;
    arcLen1 = (a * e - b * d) / den;

    Eigen::Vector3f poca0 = P0 + arcLen0 * u0;
    Eigen::Vector3f poca1 = P1 + arcLen1 * u1;

    return (poca0 - poca1).norm();
}

float SnippetHit3DBuilderSBN::chargeIntegral(float peakMean,
                                           float peakAmp,
                                           float peakSigma,
                                           float areaNorm,
                                           int   low,
                                           int   hi) const
{
    float integral(0);

    for(int sigPos = low; sigPos < hi; sigPos++)
    {
        float arg = (float(sigPos) - peakMean + 0.5) / peakSigma;
        integral += peakAmp * std::exp(-0.5 * arg * arg);
    }

    return integral;
}

bool SnippetHit3DBuilderSBN::makeDeadChannelPair(reco::ClusterHit3D&       pairOut,
                                                    const reco::ClusterHit3D& pair,
                                                    MatchedHitMap&            matchedHitMap,
                                                    size_t                    maxChanStatus,
                                                    size_t                    minChanStatus,
                                                    float                     minOverlap) const
{           
    // Assume failure (most common result)
    bool result(false);

    const reco::ClusterHit2D* hit0 = pair.getHits()[0];
    const reco::ClusterHit2D* hit1 = pair.getHits()[1];

    size_t missPlane(2);

    // u plane hit is missing
    if (!hit0)
    {
        hit0      = pair.getHits()[2];
        missPlane = 0;
    }
    // v plane hit is missing
    else if (!hit1)
    {
        hit1      = pair.getHits()[2];
        missPlane = 1;
    }

    // It can be the case that these hits have already been associated to a space point
    if (matchedHitMap[hit0].find(hit1) != matchedHitMap[hit0].end())
    {
        // Which plane is missing?
        geo::WireID wireID0 = hit0->WireID();
        geo::WireID wireID1 = hit1->WireID();

        // Ok, recover the wireID expected in the third plane...
        geo::WireID wireIn(wireID0.Cryostat,wireID0.TPC,missPlane,0);
        geo::WireID wireID = NearestWireID(pair.getPosition(), wireIn);

        // There can be a round off issue so check the next wire as well
        bool wireStatus    = m_channelStatus[wireID.Plane][wireID.Wire]   < maxChanStatus && m_channelStatus[wireID.Plane][wireID.Wire]   >= minChanStatus;
        bool wireOneStatus = m_channelStatus[wireID.Plane][wireID.Wire+1] < maxChanStatus && m_channelStatus[wireID.Plane][wireID.Wire+1] >= minChanStatus;

        // Make sure they are of at least the minimum status
        if(wireStatus || wireOneStatus)
        {
            // Sort out which is the wire we're dealing with
            if (!wireStatus) wireID.Wire += 1;

            // Want to refine position since we "know" the missing wire
            if (auto widIntersect0 = m_wireReadout->WireIDsIntersect(wireID0, wireID))
            {
                if (auto widIntersect1 = m_wireReadout->WireIDsIntersect(wireID1, wireID))
                {
                    Eigen::Vector3f newPosition(pair.getPosition()[0],pair.getPosition()[1],pair.getPosition()[2]);

                    newPosition[1] = (newPosition[1] + widIntersect0->y + widIntersect1->y) / 3.;
                    newPosition[2] = (newPosition[2] + widIntersect0->z + widIntersect1->z - 2. * m_zPosOffset) / 3.;

                    pairOut = pair;
                    pairOut.setWireID(wireID);
                    pairOut.setPosition(newPosition);

                    if (hit0->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit0->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);
                    if (hit1->getStatusBits() & reco::ClusterHit2D::USEDINTRIPLET) hit1->setStatusBit(reco::ClusterHit2D::SHAREDINTRIPLET);

                    hit0->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);
                    hit1->setStatusBit(reco::ClusterHit2D::USEDINTRIPLET);

                    matchedHitMap[hit0].insert(hit1);
                    matchedHitMap[hit1].insert(hit0);

                    result  = true;
                }
            }
        }
    }

    return result;
}

const reco::ClusterHit2D* SnippetHit3DBuilderSBN::FindBestMatchingHit(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float pairDeltaTimeLimits) const
{
    static const float minCharge(0.);

    const reco::ClusterHit2D* bestVHit(0);

    float pairAvePeakTime(pair.getAvePeakTime());

    // Idea is to loop through the input set of hits and look for the best combination
    for (const auto& hit2D : hit2DSet)
    {
        if (hit2D->getHit()->Integral() < minCharge) continue;

        float hitVPeakTime(hit2D->getTimeTicks());
        float deltaPeakTime(pairAvePeakTime-hitVPeakTime);

        if (deltaPeakTime >  pairDeltaTimeLimits) continue;

        if (deltaPeakTime < -pairDeltaTimeLimits) break;

        pairDeltaTimeLimits = fabs(deltaPeakTime);
        bestVHit            = hit2D;
    }

    return bestVHit;
}

int SnippetHit3DBuilderSBN::FindNumberInRange(const Hit2DSet& hit2DSet, const reco::ClusterHit3D& pair, float range) const
{
    static const float minCharge(0.);

    int    numberInRange(0);
    float pairAvePeakTime(pair.getAvePeakTime());

    // Idea is to loop through the input set of hits and look for the best combination
    for (const auto& hit2D : hit2DSet)
    {
        if (hit2D->getHit()->Integral() < minCharge) continue;

        float hitVPeakTime(hit2D->getTimeTicks());
        float deltaPeakTime(pairAvePeakTime-hitVPeakTime);

        if (deltaPeakTime >  range) continue;

        if (deltaPeakTime < -range) break;

        numberInRange++;
    }

    return numberInRange;
}

geo::WireID SnippetHit3DBuilderSBN::NearestWireID(const Eigen::Vector3f& position, const geo::WireID& wireIDIn) const
{
    int wire(0);

    // Embed the call to the geometry's services nearest wire id method in a try-catch block
    try
    {
        // Not sure the thinking above but is wrong... switch back to NearestWireID...
        wire = m_wireReadout->NearestWireID(geo::vect::toPoint(position.data()),wireIDIn).Wire;
    }
    catch(std::exception& exc) 
    {
        // This can happen, almost always because the coordinates are **just** out of range
        mf::LogWarning("Cluster3D") << "Exception caught finding nearest wire, position - " << exc.what() << std::endl;

        // Assume extremum for wire number depending on z coordinate
        if (position[2] < m_geometry->TPC().HalfLength()) wire = 0;
        else                                              wire = m_wireReadout->Nwires(wireIDIn.asPlaneID()) - 1;
    }

    geo::WireID wireID(wireIDIn.Cryostat,wireIDIn.TPC,wireIDIn.Plane,wire);

    return wireID;
}

float SnippetHit3DBuilderSBN::DistanceFromPointToHitWire(const Eigen::Vector3f& position, const geo::WireID& wireIDIn) const
{
    float distance = std::numeric_limits<float>::max();

    // Embed the call to the geometry's services nearest wire id method in a try-catch block
    try
    {
        // Recover wire geometry information for each wire
        const geo::WireGeo& wireGeo = m_wireReadout->Wire(wireIDIn);

        // Get wire position and direction for first wire
        auto const wirePosArr = wireGeo.GetCenter();

        Eigen::Vector3f wirePos(wirePosArr.X(),wirePosArr.Y(),wirePosArr.Z());
        Eigen::Vector3f wireDir(wireGeo.Direction().X(),wireGeo.Direction().Y(),wireGeo.Direction().Z());

        //*********************************
        // Kludge
//        if (wireIDIn.Plane > 0) wireDir[2] = -wireDir[2];


        // Want the hit position to have same x value as wire coordinates
        Eigen::Vector3f hitPosition(wirePos[0],position[1],position[2]);

        // Get arc length to doca
        double arcLen = (hitPosition - wirePos).dot(wireDir);

        // Make sure arclen is in range
        if (abs(arcLen) < wireGeo.HalfL())
        {
            Eigen::Vector3f docaVec = hitPosition - (wirePos + arcLen * wireDir);

            distance = docaVec.norm();
        }
    }
    catch(std::exception& exc)
    {
        // This can happen, almost always because the coordinates are **just** out of range
        mf::LogWarning("Cluster3D") << "Exception caught finding nearest wire, position - " << exc.what() << std::endl;

        // Assume extremum for wire number depending on z coordinate
        distance = 0.;
    }

    return distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool SetHitTimeOrder(const reco::ClusterHit2D* left, const reco::ClusterHit2D* right)
{
    // Sort by "modified start time" of pulse
    return left->getHit()->PeakTime() < right->getHit()->PeakTime();
}

bool Hit2DSetCompare::operator() (const reco::ClusterHit2D* left, const reco::ClusterHit2D* right) const
{
    return left->getHit()->PeakTime() < right->getHit()->PeakTime();
}

//------------------------------------------------------------------------------------------------------------------------------------------
void SnippetHit3DBuilderSBN::CollectArtHits(const art::Event& evt) const
{
    /**
     *  @brief Recover the 2D hits from art and fill out the local data structures for the 3D clustering
     */

    // Start by getting a vector of valid, non empty hit collections to make sure we really have something to do here...
    // Here is a container for the hits...
    std::vector<const recob::Hit*> recobHitVec;

    // Loop through the list of input sources
    for(const auto& inputTag : m_hitFinderTagVec)
    {
        art::Handle< std::vector<recob::Hit> > recobHitHandle;
        evt.getByLabel(inputTag, recobHitHandle);

        if (!recobHitHandle.isValid() || recobHitHandle->size() == 0) continue;

        recobHitVec.reserve(recobHitVec.size() + recobHitHandle->size());

        for(const auto& hit : *recobHitHandle) recobHitVec.push_back(&hit);
    }

    // If the vector is empty there is nothing to do
    if (recobHitVec.empty()) return;

    cet::cpu_timer theClockMakeHits;

    if (m_enableMonitoring) theClockMakeHits.start();

    // Need the detector properties which needs the clocks
    auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const det_prop   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);

    // Try to output a formatted string
    std::string debugMessage("");

    // Keep track of x position limits
    std::map<geo::PlaneID,double> planeIDToPositionMap;

    // Initialize the plane to hit vector map
    for(auto const& tpcID : m_geometry->Iterate<geo::TPCID>())
    {
      m_planeToSnippetHitMap[geo::PlaneID(tpcID,0)] = SnippetHitMap();
      m_planeToSnippetHitMap[geo::PlaneID(tpcID,1)] = SnippetHitMap();
      m_planeToSnippetHitMap[geo::PlaneID(tpcID,2)] = SnippetHitMap();

            // Should we provide output?
            if (!m_weHaveAllBeenHereBefore)
            {
                std::ostringstream outputString;

        outputString << "***> plane 0 offset: " << m_PlaneToT0OffsetMap.find(geo::PlaneID(tpcID,0))->second
                     << ", plane 1: " << m_PlaneToT0OffsetMap.find(geo::PlaneID(tpcID,1))->second
                     << ", plane    2: " << m_PlaneToT0OffsetMap.find(geo::PlaneID(tpcID,2))->second << "\n";
        outputString << "     Det prop plane 0: " << det_prop.GetXTicksOffset(geo::PlaneID(tpcID,0)) << ", plane 1: "  << det_prop.GetXTicksOffset(geo::PlaneID(tpcID,1)) << ", plane 2: " << det_prop.GetXTicksOffset(geo::PlaneID(tpcID,2)) << ", Trig: " << trigger_offset(clock_data) << "\n";
                debugMessage += outputString.str() + "\n";
            }

      geo::PlaneID const planeID2{tpcID, 2};
      double xPosition(det_prop.ConvertTicksToX(0., planeID2));
      planeIDToPositionMap[planeID2] = xPosition;
    }

    if (!m_weHaveAllBeenHereBefore)
    {
        for(const auto& planeToPositionPair : planeIDToPositionMap)
        {
            std::ostringstream outputString;

            outputString << "***> Plane " << planeToPositionPair.first << " has time=0 position: " << planeToPositionPair.second << "\n";

            debugMessage += outputString.str();
        }

        mf::LogDebug("SnippetHit3D") << debugMessage << std::endl;

        m_weHaveAllBeenHereBefore = true;
    }

    // Cycle through the recob hits to build ClusterHit2D objects and insert
    // them into the map
    for (const auto& recobHit : recobHitVec)
    {
        // Reject hits with negative charge, these are misreconstructed
        if (recobHit->Integral() < 0.) continue;

        // For some detectors we can have multiple wire ID's associated to a given channel.
        // So we recover the list of these wire IDs
        const std::vector<geo::WireID>& wireIDs = m_wireReadout->ChannelToWire(recobHit->Channel());

        // Start/End ticks to identify the snippet
        HitStartEndPair hitStartEndPair(recobHit->StartTick(),recobHit->EndTick());

        // Can this really happen?
        if (hitStartEndPair.second <= hitStartEndPair.first)
        {
            mf::LogInfo("SnippetHit3D") << "Yes, found a hit with end time less than start time: " << hitStartEndPair.first << "/" << hitStartEndPair.second << ", mult: " << recobHit->Multiplicity();
            continue;
        }

        // And then loop over all possible to build out our maps
        //for(const auto& wireID : wireIDs)
        for(auto wireID : wireIDs)
        {
            // Check if this is an invalid TPC
            // (for example, in protoDUNE there are logical TPC's which see no signal)
            if (std::find(m_invalidTPCVec.begin(),m_invalidTPCVec.end(),wireID.TPC) != m_invalidTPCVec.end()) continue;

            // Note that a plane ID will define cryostat, TPC and plane
            const geo::PlaneID& planeID = wireID.planeID();

            double hitPeakTime(recobHit->PeakTime() - m_PlaneToT0OffsetMap.find(planeID)->second);
            double xPosition(det_prop.ConvertTicksToX(recobHit->PeakTime(), planeID.Plane, planeID.TPC, planeID.Cryostat));

            m_clusterHit2DMasterList.emplace_back(0, 0., 0., xPosition, hitPeakTime, wireID, recobHit);

            m_planeToSnippetHitMap[planeID][hitStartEndPair].emplace_back(&m_clusterHit2DMasterList.back());
            m_planeToWireToHitSetMap[planeID][wireID.Wire].insert(&m_clusterHit2DMasterList.back());
        }
    }

    // Make a loop through to sort the recover hits in time order
//    for(auto& hitVectorMap : m_planeToSnippetHitMap)
//        std::sort(hitVectorMap.second.begin(), hitVectorMap.second.end(), SetHitTimeOrder);

    if (m_enableMonitoring)
    {
        theClockMakeHits.stop();

        m_timeVector[COLLECTARTHITS] = theClockMakeHits.accumulated_real_time();
    }

    mf::LogDebug("SnippetHit3D") << ">>>>> Number of ART hits: " << m_clusterHit2DMasterList.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SnippetHit3DBuilderSBN::CreateNewRecobHitCollection(art::Event&              event,
                                                            reco::HitPairList&       hitPairList,
                                                            std::vector<recob::Hit>& hitPtrVec,
                                                            RecobHitToPtrMap&        recobHitToPtrMap)
{
    // Set up the timing
    cet::cpu_timer theClockBuildNewHits;

    if (m_enableMonitoring) theClockBuildNewHits.start();

    // We want to build a unique list of hits from the input 3D points which, we know, reuse 2D points frequently
    // At the same time we need to create a new 2D hit with the "correct" WireID and replace the old 2D hit in the
    // ClusterHit2D object with this new one... while keeping track of the use of the old ones. My head is spinning...
    // Declare a set which will allow us to keep track of those CusterHit2D objects we have seen already
    std::set<const reco::ClusterHit2D*> visitedHit2DSet;

    // Use this handy art utility to make art::Ptr objects to the new recob::Hits for use in the output phase
    art::PtrMaker<recob::Hit> ptrMaker(event);

    // Reserve enough memory to replace every recob::Hit we have considered (this is upper limit)
    hitPtrVec.reserve(m_clusterHit2DMasterList.size());

    // Scheme is to loop through all 3D hits, then through each associated ClusterHit2D object
    for(reco::ClusterHit3D& hit3D : hitPairList)
    {
        reco::ClusterHit2DVec& hit2DVec = hit3D.getHits();

        // The loop is over the index so we can recover the correct WireID to associate to the new hit when made
        for(size_t idx = 0; idx < hit3D.getHits().size(); idx++)
        {
            const reco::ClusterHit2D* hit2D = hit2DVec[idx];

            if (!hit2D) continue;

            // Have we seen this 2D hit already?
            if (visitedHit2DSet.find(hit2D) == visitedHit2DSet.end())
            {
                visitedHit2DSet.insert(hit2D);

                // Create and save the new recob::Hit with the correct WireID
                hitPtrVec.emplace_back(recob::HitCreator(*hit2D->getHit(), hit3D.getWireIDs()[idx]).copy());

                // Recover a pointer to it...
                recob::Hit* newHit = &hitPtrVec.back();

                // Create a mapping from this hit to an art Ptr representing it
                recobHitToPtrMap[newHit] = ptrMaker(hitPtrVec.size()-1);

                // And set the pointer to this hit in the ClusterHit2D object
                const_cast<reco::ClusterHit2D*>(hit2D)->setHit(newHit);
            }
        }
    }

    size_t numNewHits = hitPtrVec.size();

    if (m_enableMonitoring)
    {
        theClockBuildNewHits.stop();

        m_timeVector[BUILDNEWHITS] = theClockBuildNewHits.accumulated_real_time();
    }

    mf::LogDebug("SnippetHit3D") << ">>>>> New output recob::Hit size: " << numNewHits << " (vs " << m_clusterHit2DMasterList.size() << " input)" << std::endl;

    return;
}

void SnippetHit3DBuilderSBN::makeWireAssns(const art::Event& evt, art::Assns<recob::Wire, recob::Hit>& wireAssns, RecobHitToPtrMap& recobHitPtrMap) const
{
    // Let's make sure the input associations container is empty
    wireAssns = art::Assns<recob::Wire, recob::Hit>();

    // First task is to recover all of the previous wire <--> hit associations and map them by channel number
    // Create the temporary container
    std::unordered_map<raw::ChannelID_t, art::Ptr<recob::Wire>> channelToWireMap;

    // Go through the list of input sources and fill out the map
    for(const auto& inputTag : m_hitFinderTagVec)
    {
        art::ValidHandle<std::vector<recob::Hit>> hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(inputTag);

        art::FindOneP<recob::Wire> hitToWireAssns(hitHandle, evt, inputTag);

        if (hitToWireAssns.isValid())
        {
            for(size_t wireIdx = 0; wireIdx < hitToWireAssns.size(); wireIdx++)
            {
                art::Ptr<recob::Wire> wire = hitToWireAssns.at(wireIdx);

                channelToWireMap[wire->Channel()] = wire;
            }
        }
    }

    // Now fill the container
    for(const auto& hitPtrPair : recobHitPtrMap)
    {
        raw::ChannelID_t channel = hitPtrPair.first->Channel();

        std::unordered_map<raw::ChannelID_t, art::Ptr<recob::Wire>>::iterator chanWireItr = channelToWireMap.find(channel);

        if (!(chanWireItr != channelToWireMap.end()))
        {
            //mf::LogDebug("SnippetHit3D") << "** Did not find channel to wire match! Skipping..." << std::endl;
            continue;
        }

        wireAssns.addSingle(chanWireItr->second, hitPtrPair.second);
    }

    return;
}

void SnippetHit3DBuilderSBN::makeRawDigitAssns(const art::Event& evt, art::Assns<raw::RawDigit, recob::Hit>& rawDigitAssns, RecobHitToPtrMap& recobHitPtrMap) const
{
    // Let's make sure the input associations container is empty
    rawDigitAssns = art::Assns<raw::RawDigit, recob::Hit>();

    // First task is to recover all of the previous wire <--> hit associations and map them by channel number
    // Create the temporary container
    std::unordered_map<raw::ChannelID_t, art::Ptr<raw::RawDigit>> channelToRawDigitMap;

    // Go through the list of input sources and fill out the map
    for(const auto& inputTag : m_hitFinderTagVec)
    {
        art::ValidHandle<std::vector<recob::Hit>> hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(inputTag);

        art::FindOneP<raw::RawDigit> hitToRawDigitAssns(hitHandle, evt, inputTag);

        if (hitToRawDigitAssns.isValid())
        {
            for(size_t rawDigitIdx = 0; rawDigitIdx < hitToRawDigitAssns.size(); rawDigitIdx++)
            {
                art::Ptr<raw::RawDigit> rawDigit = hitToRawDigitAssns.at(rawDigitIdx);

                channelToRawDigitMap[rawDigit->Channel()] = rawDigit;
            }
        }
    }

    // Now fill the container
    for(const auto& hitPtrPair : recobHitPtrMap)
    {
        raw::ChannelID_t channel = hitPtrPair.first->Channel();

        std::unordered_map<raw::ChannelID_t, art::Ptr<raw::RawDigit>>::iterator chanRawDigitItr = channelToRawDigitMap.find(channel);

        if (chanRawDigitItr == channelToRawDigitMap.end())
        {
            //mf::LogDebug("SnippetHit3D") << "** Did not find channel to wire match! Skipping..." << std::endl;
           continue;
        }

        rawDigitAssns.addSingle(chanRawDigitItr->second, hitPtrPair.second);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

DEFINE_ART_CLASS_TOOL(SnippetHit3DBuilderSBN)
} // namespace lar_cluster3d
