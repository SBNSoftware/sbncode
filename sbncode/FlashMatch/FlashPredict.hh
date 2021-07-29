#ifndef SBN_FLASHMATCH_FLASHPREDICT_HH
#define SBN_FLASHMATCH_FLASHPREDICT_HH

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
// #include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/OpHit.h"
// #include "lardataobj/RecoBase/OpFlash.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
// #include "nusimdata/SimulationBase/MCParticle.h" //uncomment for creating metrics

#include "TTree.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

// #include "sbncode/CAFMaker/RecoUtils/RecoUtils.h" //uncomment for creating metrics
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderTypes.h"
#include "sbncode/OpDet/PDMapAlg.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"

// #include "nusimdata/SimulationBase/MCParticle.h" //uncomment for creating metrics
// #include "nusimdata/SimulationBase/MCTruth.h" //uncomment for creating metrics

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
    double mcT0 = -9999;
    ChargeDigest() = default;
    ChargeDigest(const size_t pId_, const int pfpPDGC_,
                 const art::Ptr<recob::PFParticle>& pfp_ptr_,
                 const flashmatch::QCluster_t& qClusters_,
                 const std::set<unsigned>& tpcWithHits_) :
      pId(pId_), pfpPDGC(pfpPDGC_), pfp_ptr(pfp_ptr_),
      qClusters(qClusters_), tpcWithHits(tpcWithHits_)
      {}
    ChargeDigest(const size_t pId_, const int pfpPDGC_,
                 const art::Ptr<recob::PFParticle>& pfp_ptr_,
                 const flashmatch::QCluster_t& qClusters_,
                 const std::set<unsigned>& tpcWithHits_,
                 const double mcT0_) :
      pId(pId_), pfpPDGC(pfpPDGC_), pfp_ptr(pfp_ptr_),
      qClusters(qClusters_), tpcWithHits(tpcWithHits_), mcT0(mcT0_)
      {}
  };
  using ChargeDigestMap = std::map<double, ChargeDigest, std::greater<double>>;

  using OpHitIt = std::vector<recob::OpHit>::iterator;
  struct SimpleFlash {
    unsigned flashId, ophsInVolume;
    OpHitIt opH_beg, opH_end;
    double maxpeak_time;
    SimpleFlash(unsigned flashId_, unsigned ophsInVolume_,
                OpHitIt opH_beg_, OpHitIt opH_end_,
                double maxpeak_time_) :
      flashId(flashId_), ophsInVolume(ophsInVolume_),
      opH_beg(opH_beg_), opH_end(opH_end_),
      maxpeak_time(maxpeak_time_)
      {}
  };

  struct ChargeMetrics {
    double x, x_gl, y, z, q;
    bool metric_ok;
    ChargeMetrics(double x_, double x_gl_, double y_, double z_,
                  double q_, bool metric_ok_) :
      x(x_), x_gl(x_gl_), y(y_), z(z_), q(q_), metric_ok(metric_ok_)
      {}
    // faulty charges constructor
    ChargeMetrics() : metric_ok(false) {}
    std::string dumpMetrics() const
      {
        std::ostringstream stream;
        stream
          << "  x:       " << x << "\n"
          << "  x_gl:    " << x_gl << "\n"
          << "  y:       " << y << "\n"
          << "  z:       " << z << "\n"
          << "  q:       " << q << "\n"
          << "  metric_ok: " << metric_ok << "\n";
        return stream.str();
      }
  };

  struct FlashMetrics {
    double x, x_gl, y, yb, z, zb;
    double rr, pe, unpe, ratio, time;
    double h_x, h_xerr, h_xrr, h_xratio;
    double y_skew, z_skew;
    bool metric_ok;
    FlashMetrics(double x_, double x_gl_, double y_, double yb_, double z_, double zb_,
                 double rr_, double pe_, double unpe_, double ratio_, double time_,
                 double h_x_, double h_xerr_, double h_xrr_, double h_xratio_,
                 double y_skew_, double z_skew_,
                 bool metric_ok_) :
      x(x_), x_gl(x_gl_), y(y_), yb(yb_), z(z_), zb(zb_),
      rr(rr_), pe(pe_), unpe(unpe_), ratio(ratio_), time(time_),
      h_x(h_x_), h_xerr(h_xerr_), h_xrr(h_xrr_), h_xratio(h_xratio_),
      y_skew(y_skew_), z_skew(z_skew_),
      metric_ok(metric_ok_)
      {}
    // faulty flashes constructor
    FlashMetrics() :
      x(0.), x_gl(0.), y(0.), yb(0.), z(0.), zb(0.),
      rr(0.), pe(0.), unpe(0.), ratio(0.), time(0.),
      h_x(0.), h_xerr(0.), h_xrr(0.), h_xratio(0.),
      y_skew(0.), z_skew(0.),
      metric_ok(false) {}
    std::string dumpMetrics() const
      {
        std::ostringstream stream;
        stream
          << "  x:       " << x << "\n"
          << "  x_gl:    " << x_gl << "\n"
          << "  y:       " << y << "\n"
          << "  yb:      " << yb << "\n"
          << "  z:       " << z << "\n"
          << "  zb:      " << zb << "\n"
          << "  rr:      " << rr << "\n"
          << "  pe:      " << pe << "\n"
          << "  unpe:    " << unpe << "\n"
          << "  ratio:   " << ratio << "\n"
          << "  time:    " << time << "\n"
          << "  h_x:     " << h_x << "\n"
          << "  h_xerr:  " << h_xerr << "\n"
          << "  h_xrr:   " << h_xrr << "\n"
          << "  h_xratio:" << h_xratio << "\n"
          << "  y_skew:  " << y_skew << "\n"
          << "  z_skew:  " << z_skew << "\n"
          << "  metric_ok: " << std::boolalpha << metric_ok << "\n";
        return stream.str();
      }
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
  void initTree(void);
  void loadMetrics(void);
  double cheatMCT0(const std::vector<art::Ptr<recob::Hit>>& hits,
                   const std::vector<art::Ptr<simb::MCParticle>>& mcParticles);
  ChargeMetrics computeChargeMetrics(
    const flashmatch::QCluster_t& qClusters) const;
  FlashMetrics computeFlashMetrics(const SimpleFlash& simpleFlash) const;
  Score computeScore(const ChargeMetrics& charge,
                     const FlashMetrics& flash,
                     const std::set<unsigned>& tpcWithHits,
                     const int pdgc) const;
  std::tuple<double, double, double, double> hypoFlashX_fits(
    double flash_rr, double flash_ratio) const;
  std::tuple<double, double, double, double> hypoFlashX_H2(
    double flash_rr, double flash_ratio) const;
  std::tuple<double, double> xEstimateAndRMS(
    double metric_value, const TH2D* metric_h2) const;
  ChargeDigestMap makeChargeDigest(
    const art::Event& evt,
    const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h);
  void addDaughters(const std::unordered_map<size_t, size_t>& pfpMap,
                    const art::Ptr<recob::PFParticle>& pfp_ptr,
                    const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h,
                    std::vector<art::Ptr<recob::PFParticle>>& pfp_v) const;
  void updateChargeMetrics(const ChargeMetrics& chargeMetrics);
  void updateFlashMetrics(const FlashMetrics& flashMetrics);
  void updateScore(const Score& score);
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
    std::vector<recob::OpHit>& opHitsRght,
    std::vector<recob::OpHit>& opHitsLeft) const;
  std::vector<SimpleFlash> makeSimpleFlashes(//ICARUS overload
    std::vector<recob::OpHit>& opHits) const;
  bool createOpHitsTimeHist(//SBND overload
    const std::vector<recob::OpHit>& opHits,
    std::vector<recob::OpHit>& opHitsRght,
    std::vector<recob::OpHit>& opHitsLeft,
    std::unique_ptr<TH1D>& opHitsTimeHist,
    std::unique_ptr<TH1D>& opHitsTimeHistRght,
    std::unique_ptr<TH1D>& opHitsTimeHistLeft) const;
  unsigned createOpHitsTimeHist(//ICARUS overload
    const std::vector<recob::OpHit>& opHits,
    std::unique_ptr<TH1D>& opHitsTimeHist) const;
  bool findSimpleFlashes(
    std::vector<SimpleFlash>& simpleFlashes,
    std::vector<recob::OpHit>& opHits,
    const unsigned ophsInVolume,
    std::unique_ptr<TH1D>& opHitsTimeHist) const;
  inline std::string detectorName(const std::string detName) const;
  bool isPDInCryo(const int pdChannel) const;
  bool isSBNDPDRelevant(const int pdChannel,
                        const std::set<unsigned>& tpcWithHits) const;
  unsigned sbndPDinTPC(const int pdChannel) const;
  unsigned icarusPDinTPC(const int pdChannel) const;
  double wallXWithMaxPE(const OpHitIt opH_beg,
                        const OpHitIt opH_end) const;
  std::list<double> wiresXGl() const;
  double driftDistance() const;
  double flashXGl(const double hypo_x, const double flash_x) const;
  double foldXGl(const double x_gl) const;
  unsigned driftVolume(const double x) const;
  template <typename Stream>
  void printBookKeeping(Stream&& out);
  void updateBookKeeping();
  template <typename Stream>
  void printMetrics(const std::string metric,
                    const ChargeMetrics& charge,
                    const FlashMetrics& flash,
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
  const bool fForceConcurrence;
  const bool fUseUncoatedPMT, fUseOppVolMetric;//, fUseCalo;
  const bool fUseARAPUCAS;
  const bool fStoreCheatMCT0;
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
  const size_t fNTPC;
  const int fCryostat;  // =0 or =1 to match ICARUS reco chain selection
  // geo::CryostatID fCryostat;  // TODO: use this type instead
  const std::unique_ptr<geo::CryostatGeo> fGeoCryo;
  const std::list<double> fWiresX_gl;
  const double fDriftDistance;
  const int fXBins;
  const double fXBinWidth;
  const std::string fRR_TF1_fit, fRatio_TF1_fit;
  unsigned fYBins,fZBins;
  double fYLow, fYHigh, fZLow, fZHigh;
  double fYSkewThreshold, fZSkewThreshold;
  const double fYBiasSlope, fZBiasSlope;
  unsigned fDriftVolumes;
  unsigned fTPCPerDriftVolume;
  const unsigned fOpDetNormalizer;
  const double fTermThreshold;

  static constexpr unsigned kRght = 0;
  static constexpr unsigned kLeft = 1;

  static constexpr unsigned kActivityInRght = 100;
  static constexpr unsigned kActivityInLeft = 200;
  static constexpr unsigned kActivityInBoth = 300;

  // root stuff
  TTree* _flashmatch_nuslice_tree;
  struct Fits {
    double min, max;
    std::unique_ptr<TF1> f;
  };
  TH2D* fRRH2; TH2D* fRatioH2;
  static constexpr unsigned kMinEntriesInProjection = 100;
  std::array<Fits, 3> fRRFits;
  std::array<Fits, 3> fRatioFits;
  const std::array<std::string, 3> kSuffixes{"l", "h", "m"};// low, high, medium
  static constexpr double kEps = 1e-4;

  // Tree variables
  double _charge_x_gl, _charge_x,
    _charge_y, _charge_z, _charge_q;
  double _flash_x, _flash_x_gl, _flash_y, _flash_yb, _flash_z, _flash_zb,
    _flash_rr, _flash_pe, _flash_unpe, _flash_ratio, _flash_time,
    _hypo_x, _hypo_x_err, _hypo_x_rr, _hypo_x_ratio,
    _y_skew, _z_skew;
  double _score, _scr_y, _scr_z, _scr_rr, _scr_ratio;
  unsigned _evt, _run, _sub;
  unsigned _slices = -1; unsigned _true_nus = -1;
  double _mcT0 = -9999.;

  std::vector<double> fdYMeans, fdZMeans, fRRMeans, fRatioMeans;
  std::vector<double> fdYSpreads, fdZSpreads, fRRSpreads, fRatioSpreads;

  static constexpr bool kNoScr = false;
  static constexpr double kNoScrTime = -9999.;
  static constexpr double kNoScrQ  = -9999.;
  static constexpr double kNoScrPE = -9999.;
  static constexpr int kQNoOpHScr = -1;
  static constexpr int kNoChrgScr = -2;
  static constexpr int k0VUVPEScr = -3;
  static constexpr int kNoOpHInEvt = -11;
  static constexpr int kNoPFPInEvt = -12;
  static constexpr int kNoSlcInEvt = -13;
  struct BookKeeping {
    int job_bookkeeping, events_processed;
    unsigned events, nopfpneutrino, noslice,
      nullophittime, nonvalidophit;

    int pfp_bookkeeping, scored_pfp;
    unsigned pfp_to_score, no_charge, no_oph_hits,
      no_flash_pe;
  };
  BookKeeping bk;
};


#endif //SBN_FLASHMATCH_FLASHPREDICT_HH
