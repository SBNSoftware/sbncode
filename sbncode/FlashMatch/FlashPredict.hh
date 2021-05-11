#ifndef SBN_FLASHMATCH_FLASHPREDICT_HH
#define SBN_FLASHMATCH_FLASHPREDICT_HH

// save diagnostic state
//#pragma GCC diagnostic push

// turn off the specific warning. Can also use "-Wall"
//#pragma GCC diagnostic ignored "-Wconversion"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
// #include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TTree.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TSpline.h"

#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderTypes.h"
#include "sbncode/OpDet/PDMapAlg.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"

// turn the warnings back on
//#pragma GCC diagnostic pop

#include <algorithm>
#include <iterator>
#include <list>
#include <memory>
#include <set>
#include <string>

class FlashPredict;
class FlashPredict : public art::EDProducer {


public:
  explicit FlashPredict(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  FlashPredict(FlashPredict const&) = delete;
  FlashPredict(FlashPredict&&) = delete;
  FlashPredict& operator=(FlashPredict const&) = delete;
  FlashPredict& operator=(FlashPredict&&) = delete;
  // Required functions.
  void produce(art::Event& evt) override;
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;


private:
  // Declare member data here.
  //  ::flashmatch::FlashMatchManager m_flashMatchManager; ///< The flash match manager
  // art::InputTag fFlashProducer;
  // art::InputTag fT0Producer; // producer for ACPT in-time anab::T0 <-> recob::Track assocaition
  void initTree(void);
  void loadMetrics(void);
  bool computeChargeMetrics(const flashmatch::QCluster_t& qClusters);
  bool computeFlashMetrics(const std::set<unsigned>& tpcWithHits);
  bool computeScore(const std::set<unsigned>& tpcWithHits, const int pdgc);
  double hypoFlashX_splines() const;
  double hypoFlashX_fits();
  // ::flashmatch::Flash_t GetFlashPESpectrum(const recob::OpFlash& opflash);
  // void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap& pfParticleMap,
  //                                   const art::Ptr<recob::PFParticle>& particle,
  //                                   lar_pandora::PFParticleVector& downstreamPFParticles) const;
  // void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap& pfParticleMap,
  //                                   const lar_pandora::PFParticleVector& parentPFParticles,
  //                                   lar_pandora::PFParticleVector& downstreamPFParticles) const;
  void AddDaughters(const std::map<size_t, size_t>& pfpMap,
                    const art::Ptr<recob::PFParticle>& pfp_ptr,
                    const art::ValidHandle<std::vector<recob::PFParticle>>& pfp_h,
                    std::vector<art::Ptr<recob::PFParticle>>& pfp_v) const;
  double scoreTerm(const double m, const double n,
                   const double mean, const double spread) const;
  double scoreTerm(const double m,
                   const double mean, const double spread) const;
  bool pfpNeutrinoOnEvent(
    const art::ValidHandle<std::vector<recob::PFParticle>>& pfp_h) const;
  void copyOpHitsInBeamWindow(
    std::vector<recob::OpHit>& opHits,
    const art::Handle<std::vector<recob::OpHit>>& ophit_h) const;
  bool getOpHitsInFlash(std::vector<recob::OpHit>& opHits);
  bool createOpHitsTimeHist(
    const std::vector<recob::OpHit>& opHits) const;
  bool findMaxPeak(std::vector<recob::OpHit>& opHits);
  bool isPDInCryo(const int pdChannel) const;
  bool isSBNDPDRelevant(const int pdChannel,
                        const std::set<unsigned>& tpcWithHits) const;
  unsigned sbndPDinTPC(const int pdChannel) const;
  unsigned icarusPDinTPC(const int pdChannel) const;
  double flashXGl(const double hypo_x, const double flash_x) const;
  double driftDistance(const double x) const;
  unsigned driftVolume(const double x) const;
  // bool isPDInCryoTPC(double pd_x, size_t itpc);
  // bool isPDInCryoTPC(int pdChannel, size_t itpc);
  // bool isChargeInCryoTPC(double qp_x, int icryo, int itpc);
  template <typename Stream>
  void printBookKeeping(Stream&& out);
  void updateBookKeeping();
  template <typename Stream>
  void printMetrics(const std::string metric,
                    const int pdgc,
                    const std::set<unsigned>& tpcWithHits,
                    const double term,
                    Stream&& out) const;

  const art::InputTag fPandoraProducer, fSpacePointProducer,
    fOpHitProducer, fOpHitARAProducer;//, fCaloProducer, fTrackProducer;
  detinfo::DetectorClocksData const fClockData;
  const double fTickPeriod;
  const double fBeamWindowStart, fBeamWindowEnd;
  const double fLightWindowStart, fLightWindowEnd;
  const unsigned fTimeBins;
  const bool fSelectNeutrino, fUseUncoatedPMT, fUseOppVolMetric;//, fUseCalo;
  const bool fUseARAPUCAS;
  const std::string fInputFilename;
  const bool fNoAvailableMetrics, fMakeTree;
  const double fMinFlashPE, fMinOpHPE, fPEscale,
    fChargeToNPhotonsShower, fChargeToNPhotonsTrack;
  std::string fDetector; // SBND or ICARUS
  bool fSBND, fICARUS;
  std::unique_ptr<opdet::PDMapAlg> fPDMapAlgPtr;
  const int fCryostat;  // =0 or =1 to match ICARUS reco chain selection
  // geo::CryostatID fCryostat;  // TODO: use this type instead
  std::unique_ptr<geo::CryostatGeo> fGeoCryo;
  const int fNBins;
  const double fDriftDistance;
  size_t fNTPC;
  unsigned fDriftVolumes;
  unsigned fTPCPerDriftVolume;
  const unsigned fVUVToVIS;
  const double fTermThreshold;
  std::list<double> fWiresX_gl;
  // std::vector<double> fPMTChannelCorrection;

  unsigned fPeakCounter = 0;
  std::vector<recob::OpHit>::iterator fOpH_beg, fOpH_end;
  const art::ServiceHandle<geo::Geometry> geometry;

  // root stuff
  TTree* _flashmatch_nuslice_tree;
  std::unique_ptr<TH1D> fOpHitsTimeHist;
  // double rrMax, peMax;
  // TSpline3 rr_m_InvSpl, rr_h_InvSpl, rr_l_InvSpl;
  // TSpline3 pe_m_InvSpl, pe_h_InvSpl, pe_l_InvSpl;
  struct Fits {
    double min, max;
    std::unique_ptr<TF1> f;
  };
  std::array<Fits, 3> rrFits;
  std::array<Fits, 3> peFits;
  const std::array<std::string, 3> suffixes{"l", "h", "m"};
  // const std::string kPolFit = "pol3";
  const double kEps = 1e-4;

  // std::vector<double> _pe_reco_v, _pe_hypo_v;

  // Tree variables
  double _charge_x_gl, _charge_x,
    _charge_y, _charge_z, _charge_q;
  double _flash_x, _flash_x_gl, _flash_y, _flash_z,
    _flash_r, _flash_pe, _flash_unpe, _flash_ratio;
  // TODO: why not charge_time?
  double _flash_time;
  double _score, _scr_y, _scr_z, _scr_rr, _scr_ratio;
  double _hypo_x, _hypo_x_rr, _hypo_x_ratio;//, _hypo_x_fit;
  unsigned _evt, _run, _sub, _countPE; //_slices;

  std::vector<double> dy_means, dz_means, rr_means, pe_means;
  std::vector<double> dy_spreads, dz_spreads, rr_spreads, pe_spreads;

  struct ChargeDigest {
     size_t pId;
     int pfpPDGC;
     art::Ptr<recob::PFParticle> pfp_ptr;
     flashmatch::QCluster_t qClusters;
     std::set<unsigned> tpcWithHits;
    ChargeDigest() = default;
    ChargeDigest(const size_t pId_, const int pfpPDGC_,
                 const art::Ptr<recob::PFParticle>& pfp_ptr_,
                 const flashmatch::QCluster_t& qClusters_,
                 const std::set<unsigned>& tpcWithHits_) :
      pId(pId_), pfpPDGC(pfpPDGC_), pfp_ptr(pfp_ptr_),
      qClusters(qClusters_), tpcWithHits(tpcWithHits_)
      {}
  };

  const bool kNoScr = false;
  const double kNoScrTime = -9999.;
  const double kNoScrPE = -9999.;
  const int kQNoOpHScr = -1;
  const int kNoChrgScr = -2;
  const int k0VUVPEScr = -3;
  const int kNoOpHInEvt = -11;
  const int kNoPFPInEvt = -12;
  // const int kNoSlcInEvt = -13;
  struct BookKeeping {
    int job_bookkeeping, events_processed;
    unsigned events, nopfpneutrino,// noslice,
      nullophittime, nonvalidophit;

    int pfp_bookkeeping, scored_pfp;
    unsigned pfp_to_score, no_charge, no_oph_hits,
      no_flash_pe;
  };
  BookKeeping bk;
};


#endif //SBN_FLASHMATCH_FLASHPREDICT_HH
