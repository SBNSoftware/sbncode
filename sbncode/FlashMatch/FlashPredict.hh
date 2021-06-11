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
#include <limits>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <tuple>

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
  using ChargeDigestMap = std::map<double, ChargeDigest, std::greater<double>>;

  using OpHitIt = std::vector<recob::OpHit>::iterator;
  struct SimpleFlash {
    unsigned flashId, flashUId;
    OpHitIt opH_beg, opH_end;
    double maxpeak_time;
    SimpleFlash(unsigned flashId_, unsigned flashUId_,
                OpHitIt opH_beg_, OpHitIt opH_end_,
                double maxpeak_time_) :
      flashId(flashId_), flashUId(flashUId_),
      opH_beg(opH_beg_), opH_end(opH_end_),
      maxpeak_time(maxpeak_time_)
      {}
  };

  struct FlashMetrics {
    double x, x_gl, y, z, r, pe, unpe, ratio, time;
    double hypo, hypo_rr, hypo_ratio;
    bool metric_ok;
    FlashMetrics(double x_, double x_gl_, double y_, double z_, double r_,
                 double pe_, double unpe_, double ratio_, double time_,
                 double hypo_, double hypo_rr_, double hypo_ratio_,
                 bool metric_ok_) :
      x(x_), x_gl(x_gl_), y(y_), z(z_), r(r_), pe(pe_), unpe(unpe_),
      ratio(ratio_), time(time_), hypo(hypo_), hypo_rr(hypo_rr_),
      hypo_ratio(hypo_ratio_), metric_ok(metric_ok_)
      {}
    // faulty flashes constructor
    FlashMetrics() : metric_ok(false) {}
  };


private:
  // aliases for the objects that are stored
  using sFM    = sbn::SimpleFlashMatch;
  using Charge = sbn::SimpleFlashMatch::Charge;
  using Flash  = sbn::SimpleFlashMatch::Flash;
  using Score  = sbn::SimpleFlashMatch::Score;

  // Declare member data here.
  //  ::flashmatch::FlashMatchManager m_flashMatchManager; ///< The flash match manager
  // art::InputTag fFlashProducer;
  // art::InputTag fT0Producer; // producer for ACPT in-time anab::T0 <-> recob::Track assocaition
  void initTree(void);
  void loadMetrics(void);
  bool computeChargeMetrics(const flashmatch::QCluster_t& qClusters);
  bool computeFlashMetrics(const SimpleFlash& simpleFlash);
  Score computeScore(const std::set<unsigned>& tpcWithHits,
                    const int pdgc) const;
  // double hypoFlashX_splines() const;
  std::tuple<double, double, double, double> hypoFlashX_fits(
    double flash_r, double flash_ratio) const;
  // ::flashmatch::Flash_t GetFlashPESpectrum(const recob::OpFlash& opflash);
  // void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap& pfParticleMap,
  //                                   const art::Ptr<recob::PFParticle>& particle,
  //                                   lar_pandora::PFParticleVector& downstreamPFParticles) const;
  // void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap& pfParticleMap,
  //                                   const lar_pandora::PFParticleVector& parentPFParticles,
  //                                   lar_pandora::PFParticleVector& downstreamPFParticles) const;
  ChargeDigestMap makeChargeDigest(
    const art::Event& evt,
    const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h);
  void addDaughters(const std::unordered_map<size_t, size_t>& pfpMap,
                    const art::Ptr<recob::PFParticle>& pfp_ptr,
                    const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h,
                    std::vector<art::Ptr<recob::PFParticle>>& pfp_v) const;
  void updateFlashMetrics(const FlashMetrics& flashMetrics);
  void storeFlashMetrics(
    const unsigned flashUId,
    std::map<unsigned, FlashMetrics>& flashMetricsMap) const;
  inline double scoreTerm(const double m, const double n,
                   const double mean, const double spread) const;
  inline double scoreTerm(const double m,
                   const double mean, const double spread) const;
  inline bool pfpNeutrinoOnEvent(
    const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h) const;
  void copyOpHitsInBeamWindow(
    std::vector<recob::OpHit>& opHits,
    const art::Handle<std::vector<recob::OpHit>>& ophit_h) const;
  std::vector<SimpleFlash> makeSimpleFlashes(//SBND overload
    std::vector<recob::OpHit>& opHits,
    std::vector<recob::OpHit>& opHitsLeft,
    std::vector<recob::OpHit>& opHitsRght) const;
  std::vector<SimpleFlash> makeSimpleFlashes(//ICARUS overload
    std::vector<recob::OpHit>& opHits) const;
  bool createOpHitsTimeHist(//SBND overload
    const std::vector<recob::OpHit>& opHits,
    std::vector<recob::OpHit>& opHitsLeft,
    std::vector<recob::OpHit>& opHitsRght,
    std::unique_ptr<TH1D>& opHitsTimeHist,
    std::unique_ptr<TH1D>& opHitsTimeHistLeft,
    std::unique_ptr<TH1D>& opHitsTimeHistRght) const;
  bool createOpHitsTimeHist(//ICARUS overload
    const std::vector<recob::OpHit>& opHits,
    std::unique_ptr<TH1D>& opHitsTimeHist) const;
  bool findSimpleFlashes(
    std::vector<SimpleFlash>& simpleFlashes,
    std::vector<recob::OpHit>& opHits,
    const unsigned volumeId,
    std::unique_ptr<TH1D>& opHitsTimeHist) const;
  inline std::string detectorName(const std::string detName) const;
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
  const double fFlashStart, fFlashEnd;
  const unsigned fTimeBins;
  const bool fSelectNeutrino;
  const bool fOnlyCollectionWires;
  const bool fUseUncoatedPMT, fUseOppVolMetric;//, fUseCalo;
  const bool fUseARAPUCAS;
  const std::string fInputFilename;
  const bool fNoAvailableMetrics, fMakeTree;
  const double fChargeToNPhotonsShower, fChargeToNPhotonsTrack;
  const double fMinHitQ, fMinSpacePointQ, fMinParticleQ, fMinSliceQ;
  const unsigned fMaxFlashes;
  const double fMinOpHPE, fMinFlashPE;
  const art::ServiceHandle<geo::Geometry> fGeometry;
  const std::string fDetector; // SBND or ICARUS
  const bool fSBND, fICARUS;
  const std::unique_ptr<opdet::PDMapAlg> fPDMapAlgPtr;
  const int fCryostat;  // =0 or =1 to match ICARUS reco chain selection
  // geo::CryostatID fCryostat;  // TODO: use this type instead
  const std::unique_ptr<geo::CryostatGeo> fGeoCryo;
  const int fNBins;
  const double fDriftDistance;
  const size_t fNTPC;
  unsigned fDriftVolumes;
  unsigned fTPCPerDriftVolume;
  const unsigned fOpDetNormalizer;
  const double fTermThreshold;
  std::list<double> fWiresX_gl;
  // std::vector<double> fPMTChannelCorrection;

  const unsigned kLeftOpHs = 10;
  const unsigned kRghtOpHs = 20;
  const unsigned kBothOpHs = 30;

  // root stuff
  TTree* _flashmatch_nuslice_tree;
  // double rrMax, peMax;
  // TSpline3 rr_m_InvSpl, rr_h_InvSpl, rr_l_InvSpl;
  // TSpline3 pe_m_InvSpl, pe_h_InvSpl, pe_l_InvSpl;
  struct Fits {
    double min, max;
    std::unique_ptr<TF1> f;
  };
  std::array<Fits, 3> fRRFits;
  std::array<Fits, 3> fRatioFits;
  const std::array<std::string, 3> kSuffixes{"l", "h", "m"};// low, high, medium
  // const std::string kPolFit = "pol3";
  const double kEps = 1e-4;

  // std::vector<double> _pe_reco_v, _pe_hypo_v;

  // Tree variables
  double _charge_x_gl, _charge_x,
    _charge_y, _charge_z, _charge_q;
  double _flash_x, _flash_x_gl, _flash_y, _flash_z,
    _flash_r, _flash_pe, _flash_unpe, _flash_ratio, _flash_time,
    _hypo_x, _hypo_x_err, _hypo_x_rr, _hypo_x_ratio;
  double _score, _scr_y, _scr_z, _scr_rr, _scr_ratio;
  unsigned _evt, _run, _sub; //_slices;

  std::vector<double> dy_means, dz_means, rr_means, pe_means;
  std::vector<double> dy_spreads, dz_spreads, rr_spreads, pe_spreads;

  const bool kNoScr = false;
  const double kNoScrTime = -9999.;
  const double kNoScrQ  = -9999.;
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
