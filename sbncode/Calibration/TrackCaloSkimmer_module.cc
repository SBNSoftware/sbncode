///////////////////////////////////////////////////////////////////////
// Class:       TrackCaloSkimmer
// Plugin Type: analyzer (art v3_06_03)
// File:        TrackCaloSkimmer_module.cc
//
// Generated at Mon May 17 09:46:34 2021 by Gray Putnam using cetskelgen
// from cetlib version v3_11_01.
//
// Module for creating a skim of track calorimetry reconstruction for use
// with calibrations in ICARUS.
////////////////////////////////////////////////////////////////////////

#include "TrackCaloSkimmer.h"
#include "art/Utilities/make_tool.h"
#include "sbncode/CAFMaker/FillTrue.h"
#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"

//global functions, data for fitting
const size_t MAX_N_FIT_DATA = 30;
static int N_FIT_DATA = -1;
static double FIT_RR[MAX_N_FIT_DATA];
static double FIT_DQDX[MAX_N_FIT_DATA];

void ConstResiduals(
  int &npar, 
  double *g, 
  double &result, 
  double *par, 
  int flag) {
  double ret = 0;
  double C = *par;
  for (int i = 0; i < N_FIT_DATA; i++) {
    double diff = FIT_DQDX[i] - C;
    ret += diff*diff;
  }
  result = sqrt(ret);
}

void ExpResiduals(
  int &npar, 
  double *g, 
  double &result, 
  double *par, 
  int flag) {
  double ret = 0;
  double A = par[0];
  double R = par[1];
  for (int i = 0; i < N_FIT_DATA; i++) {
    double diff = FIT_DQDX[i] - A*exp(-FIT_RR[i]/R);
    ret += diff*diff;
  }
  result = sqrt(ret);
}

sbn::TrackCaloSkimmer::~TrackCaloSkimmer() {
  delete fTrack;
}

sbn::TrackCaloSkimmer::TrackCaloSkimmer(
  fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fFitExp(2),
    fFitConst(1)
{
  // Grab config
  fproducerCRTTag = p.get< std::vector<art::InputTag> > ("CRTTagproducer", {"CRTT0Tagging"});  
  fPFPproducer = p.get< art::InputTag > ("PFPproducer","pandoraGausCryo0");
  fT0producers = p.get< std::vector<art::InputTag> > ("T0producers", {"pandoraGausCryo0"} );
  fCALOproducer = p.get< art::InputTag > ("CALOproducer");
  fTRKproducer = p.get< art::InputTag > ("TRKproducer" );
  fTRKHMproducer = p.get< art::InputTag > ("TRKHMproducer", "");

  fMCSproducerUF = p.get< art::InputTag > ("MCSproducerUF"); 
  fMCSproducerUI = p.get< art::InputTag > ("MCSproducerUI"); 
  fMCSproducerIF3D = p.get< art::InputTag > ("MCSproducerIF3D"); 
  fMCSproducerIF2DI1 = p.get< art::InputTag > ("MCSproducerIF2DI1"); 
  fMCSproducerIF2DI2 = p.get< art::InputTag > ("MCSproducerIF2DI2"); 
  fMCSproducerIF2DC = p.get< art::InputTag > ("MCSproducerIF2DC"); 
  fMCSproducerII3D = p.get< art::InputTag > ("MCSproducerII3D"); 
  fMCSproducerII2DI1 = p.get< art::InputTag > ("MCSproducerII2DI1"); 
  fMCSproducerII2DI2 = p.get< art::InputTag > ("MCSproducerII2DI2"); 
  fMCSproducerII2DC = p.get< art::InputTag > ("MCSproducerII2DC"); 

  fRangeInputtag = p.get< art::InputTag > ("RangeInputtag"); 
  fHITproducer = p.get< art::InputTag > ("HITproducer" );
  fG4producer = p.get< std::string > ("G4producer" );
  fG4Droppedproducer = p.get< std::string > ("G4Droppedproducer" );
  fSimChannelproducer = p.get< std::string > ("SimChannelproducer" );
  fRequireT0 = p.get<bool>("RequireT0", false);
  fDoTailFit = p.get<bool>("DoTailFit", true);
  fVerbose = p.get<bool>("Verbose", false);
  //fVerbose = true;
  fSilenceMissingDataProducts = p.get<bool>("SilenceMissingDataProducts", false);
  fHitRawDigitsTickCollectWidth = p.get<double>("HitRawDigitsTickCollectWidth", 50.);
  fHitRawDigitsWireCollectWidth = p.get<int>("HitRawDigitsWireCollectWidth", 5);
  fTailFitResidualRange = p.get<double>("TailFitResidualRange", 5.);
  //if (fTailFitResidualRange > 10.) std::cout << "sbn::TrackCaloSkimmer: Bad tail fit residual range config :(" << fTailFitResidualRange << "). Fits will not be meaningful.\n";
  fFillTrackEndHits = p.get<bool>("FillTrackEndHits", true);
  fTrackEndHitWireBox = p.get<float>("TrackEndHitWireBox", 60); // 20 cm
  fTrackEndHitTimeBox = p.get<float>("TrackEndHitTimeBox", 300); // about 20cm
  fRawDigitproducers = p.get<std::vector<art::InputTag>>("RawDigitproducers", {});
  std::vector<fhicl::ParameterSet> selection_tool_configs(p.get<std::vector<fhicl::ParameterSet>>("SelectionTools"), {});
  for (const fhicl::ParameterSet &p: selection_tool_configs) {
    fSelectionTools.push_back(art::make_tool<sbn::ITCSSelectionTool>(p));
  }

  // Setup meta info
  fMeta.iproc = -1;
  fMeta.ifile = -1;
  const char *process_str = std::getenv("PROCESS");
  if (process_str) {
    try {
      fMeta.iproc = std::stoi(process_str);
    }
    catch (...) {}
  }

  // Output stuff
  fTrack = new sbn::TrackInfo();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("TrackCaloSkim", "Calo Tree");
  fTree->Branch("trk", &fTrack);
}

void sbn::TrackCaloSkimmer::analyze(
  art::Event const& e) {
  unsigned evt = e.event();
  unsigned sub = e.subRun();
  unsigned run = e.run();

  fMeta.evt = evt;
  fMeta.subrun = sub;
  fMeta.run = run;
  fMeta.time = e.time().value();

  // Services
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const dprop =
  art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  // Identify which detector: can only detect either sbnd or icarus

  std::string gdml = geometry->GDMLFile();
  gdml = basename(gdml.c_str()); 

  for(unsigned int i = 0; i <gdml.size(); ++i) gdml[i] = std::tolower(gdml[i]); 

  EDet det = kNOTDEFINED;

  const bool hasSBND = ((gdml.find("sbnd") != std::string::npos) ||
			(geometry->DetectorName().find("sbnd") != std::string::npos));

  const bool hasICARUS = ((gdml.find("icarus") != std::string::npos) ||
			(geometry->DetectorName().find("icarus") != std::string::npos));

  if(hasSBND == hasICARUS) { 
    //std::cout << "TrackCaloSkimmer: Unable to automatically determine either SBND or ICARUS!" << std::endl;
    abort();
  }
 
  if(hasSBND) det = kSBND;
  if(hasICARUS) det = kICARUS;

  // Setup the volumes
  std::vector<std::vector<geo::BoxBoundedGeo>> TPCVols;
  std::vector<geo::BoxBoundedGeo> AVs;

  // First the TPC
  for (auto const &cryo: geometry->Iterate<geo::CryostatGeo>()) {
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryo.ID())) {
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
    }
    TPCVols.push_back(std::move(this_tpc_volumes));
  }

  for (const std::vector<geo::BoxBoundedGeo> &tpcs: TPCVols) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    AVs.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
  }

  // Truth information
  std::vector<art::Ptr<simb::MCParticle>> mcparticles;
  if (fG4producer.size()) {
    art::ValidHandle<std::vector<simb::MCParticle>> mcparticle_handle = e.getValidHandle<std::vector<simb::MCParticle>>(fG4producer);
    art::fill_ptr_vector(mcparticles, mcparticle_handle);
  }

  std::vector<art::Ptr<simb::MCParticle>> droppedmcparticles;
  if (fG4Droppedproducer.size()) {
    art::ValidHandle<std::vector<simb::MCParticle>> droppedmcparticle_handle = e.getValidHandle<std::vector<simb::MCParticle>>(fG4Droppedproducer);
    art::fill_ptr_vector(droppedmcparticles, droppedmcparticle_handle);
  }

  std::vector<art::Ptr<sim::SimChannel>> simchannels;
  if (fSimChannelproducer.size()) {
    art::ValidHandle<std::vector<sim::SimChannel>> simchannel_handle = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelproducer);
    art::fill_ptr_vector(simchannels, simchannel_handle);
  }
  // Reconstructed Information
  std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
  try {
    art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::fill_ptr_vector(PFParticleList, pfparticles);
  }
  catch(...) {
      //<< "PFP's with tag: " << fPFPproducer << " not present.\n";
      // Real data may have missing products -- just ignore the event
      if (fSilenceMissingDataProducts) return;
      else throw;
  }

  std::vector<art::FindManyP<anab::T0>> fmCRTTaggedT0;
  for (unsigned i_label = 0; i_label < fproducerCRTTag.size(); i_label++) {
    //std::cout << " ilabel " << i_label << fproducerCRTTag[i_label] << std::endl;
    fmCRTTaggedT0.emplace_back(PFParticleList, e, fproducerCRTTag[i_label]);
  }

  std::vector<art::FindManyP<anab::T0>> fmT0;
  for (unsigned i_label = 0; i_label < fT0producers.size(); i_label++) {
    //std::cout << " old t0 ilabel " << i_label << fT0producers[i_label] << std::endl;
    fmT0.emplace_back(PFParticleList, e, fT0producers[i_label]);
  }
  art::FindManyP<recob::SpacePoint> PFParticleSPs(PFParticleList, e, fPFPproducer);

  // Now we don't need to guard access to further data. If this is an empty event it should be caught by PFP's or Hit's
  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(fTRKproducer); 

  // Track - associated data
  art::FindManyP<recob::Track> fmTracks(PFParticleList, e, fTRKproducer);

  art::InputTag thm_label = fTRKHMproducer.empty() ? fTRKproducer : fTRKHMproducer;
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(tracks, e, thm_label);
  art::FindManyP<anab::Calorimetry> fmCalo(tracks, e, fCALOproducer);
  // Collect raw digits for saving hits
  std::vector<art::Ptr<raw::RawDigit>> rawdigitlist;
  for (const art::InputTag &t: fRawDigitproducers) {
    try {
      art::ValidHandle<std::vector<raw::RawDigit>> thisdigits = e.getValidHandle<std::vector<raw::RawDigit>>(t);
      art::fill_ptr_vector(rawdigitlist, thisdigits);
    }
    catch(...) {
      if (!fSilenceMissingDataProducts) throw;
      else {} // Allow Raw Digits to not be present
    }
  }

  //mcs uboone algorithm, final cut
  art::ValidHandle<std::vector<recob::MCSFitResult>> MCS_handleUF = e.getValidHandle<std::vector<recob::MCSFitResult>>(fMCSproducerUF);
  std::vector<art::Ptr<recob::MCSFitResult>> MCSresultsUF;
  art::fill_ptr_vector(MCSresultsUF, MCS_handleUF);

  //mcs uboone algorithm, initial cut
  art::ValidHandle<std::vector<recob::MCSFitResult>> MCS_handleUI = e.getValidHandle<std::vector<recob::MCSFitResult>>(fMCSproducerUI);
  std::vector<art::Ptr<recob::MCSFitResult>> MCSresultsUI;
  art::fill_ptr_vector(MCSresultsUI, MCS_handleUI);

  //mcs icarus algorithm, final cut, 3d
  art::ValidHandle<std::vector<recob::MCSFitResultGS>> MCS_handleIF3D = e.getValidHandle<std::vector<recob::MCSFitResultGS>>(fMCSproducerIF3D);
  std::vector<art::Ptr<recob::MCSFitResultGS>> MCSresultsIF3D;
  art::fill_ptr_vector(MCSresultsIF3D, MCS_handleIF3D);

  //mcs icarus algorithm, final cut, 2d induction1
  art::ValidHandle<std::vector<recob::MCSFitResultGS>> MCS_handleIF2DI1 = e.getValidHandle<std::vector<recob::MCSFitResultGS>>(fMCSproducerIF2DI1);
  std::vector<art::Ptr<recob::MCSFitResultGS>> MCSresultsIF2DI1;
  art::fill_ptr_vector(MCSresultsIF2DI1, MCS_handleIF2DI1);

  //mcs icarus algorithm, final cut, 2d induction2
  art::ValidHandle<std::vector<recob::MCSFitResultGS>> MCS_handleIF2DI2 = e.getValidHandle<std::vector<recob::MCSFitResultGS>>(fMCSproducerIF2DI2);
  std::vector<art::Ptr<recob::MCSFitResultGS>> MCSresultsIF2DI2;
  art::fill_ptr_vector(MCSresultsIF2DI2, MCS_handleIF2DI2);

  //mcs icarus algorithm, final cut, 2d collection
  art::ValidHandle<std::vector<recob::MCSFitResultGS>> MCS_handleIF2DC = e.getValidHandle<std::vector<recob::MCSFitResultGS>>(fMCSproducerIF2DC);
  std::vector<art::Ptr<recob::MCSFitResultGS>> MCSresultsIF2DC;
  art::fill_ptr_vector(MCSresultsIF2DC, MCS_handleIF2DC);

  //mcs icarus algorithm, initial cut, 3d
  art::ValidHandle<std::vector<recob::MCSFitResultGS>> MCS_handleII3D = e.getValidHandle<std::vector<recob::MCSFitResultGS>>(fMCSproducerII3D);
  std::vector<art::Ptr<recob::MCSFitResultGS>> MCSresultsII3D;
  art::fill_ptr_vector(MCSresultsII3D, MCS_handleII3D);

  //mcs icarus algorithm, initial cut, 2d induction1
  art::ValidHandle<std::vector<recob::MCSFitResultGS>> MCS_handleII2DI1 = e.getValidHandle<std::vector<recob::MCSFitResultGS>>(fMCSproducerII2DI1);
  std::vector<art::Ptr<recob::MCSFitResultGS>> MCSresultsII2DI1;
  art::fill_ptr_vector(MCSresultsII2DI1, MCS_handleII2DI1);

  //mcs icarus algorithm, initial cut, 2d induction2
  art::ValidHandle<std::vector<recob::MCSFitResultGS>> MCS_handleII2DI2 = e.getValidHandle<std::vector<recob::MCSFitResultGS>>(fMCSproducerII2DI2);
  std::vector<art::Ptr<recob::MCSFitResultGS>> MCSresultsII2DI2;
  art::fill_ptr_vector(MCSresultsII2DI2, MCS_handleII2DI2);

  //mcs icarus algorithm, initial cut, 2d collection
  art::ValidHandle<std::vector<recob::MCSFitResultGS>> MCS_handleII2DC = e.getValidHandle<std::vector<recob::MCSFitResultGS>>(fMCSproducerII2DC);
  std::vector<art::Ptr<recob::MCSFitResultGS>> MCSresultsII2DC;
  art::fill_ptr_vector(MCSresultsII2DC, MCS_handleII2DC);

  //range momentum
  auto const& range_handle = e.getValidHandle<std::vector<sbn::RangeP>>(fRangeInputtag);
  std::vector<art::Ptr<sbn::RangeP>> rangePs;
  art::fill_ptr_vector(rangePs, range_handle);

  // The raw digit list is not sorted, so make it into a map on the WireID
  std::map<geo::WireID, art::Ptr<raw::RawDigit>> rawdigits;
  for (const art::Ptr<raw::RawDigit> &d: rawdigitlist) {
    
    std::vector<geo::WireID> wids;
    // Handle bad channel ID
    try {
      wids = geometry->ChannelToWire(d->Channel());
    }
    catch(...) {
      continue;
    }

    // Ignore channel with no mapped wire
    if (wids.size() == 0) continue;

    // Ignore wires that are already mapped
    if (rawdigits.count(wids[0])) continue;

    rawdigits[wids[0]] = d;
  } 
  
  // Collect all hits
  art::ValidHandle<std::vector<recob::Hit>> allhit_handle = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
  std::vector<art::Ptr<recob::Hit>> allHits;
  art::fill_ptr_vector(allHits, allhit_handle);

  // And lookup the SP's
  art::FindManyP<recob::SpacePoint> allHitSPs(allHits, e, fPFPproducer);

  // Use helper functions from CAFMaker/FillTrue
  std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map;
  std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map;
  const cheat::BackTrackerService *bt = NULL;

  if (simchannels.size()) {
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    id_to_ide_map = caf::PrepSimChannels(simchannels, *geometry);
    id_to_truehit_map = caf::PrepTrueHits(allHits, clock_data, *bt_serv.get());
    bt = bt_serv.get();
  }

  // Build global track info
  std::vector<GlobalTrackInfo> track_infos;
  for (const recob::Track &t: *tracks) 
    track_infos.push_back({
      t.Start(), t.End(), t.StartDirection(), t.EndDirection(), t.ID()
    });

  for (art::Ptr<recob::PFParticle> p_pfp: PFParticleList) {
    const recob::PFParticle &pfp = *p_pfp;
    bool isCRTTagged=false;
    const std::vector<art::Ptr<recob::Track>> thisTrack = fmTracks.at(p_pfp.key());
    if (thisTrack.size() != 1) continue;

    art::Ptr<recob::Track> trkPtr = thisTrack.at(0);

    std::vector<art::Ptr<anab::Calorimetry>> emptyCaloVector;
    const std::vector<art::Ptr<anab::Calorimetry>> &calo = fmCalo.isValid() ? fmCalo.at(trkPtr.key()) : emptyCaloVector;

    std::vector<art::Ptr<recob::Hit>> emptyHitVector;
    const std::vector<art::Ptr<recob::Hit>> &trkHits = fmtrkHits.isValid() ? fmtrkHits.at(trkPtr.key()) : emptyHitVector;

    art::Ptr<recob::MCSFitResult> mcsResultUF = MCSresultsUF.at(trkPtr.key());
    art::Ptr<recob::MCSFitResult> mcsResultUI = MCSresultsUI.at(trkPtr.key());
    art::Ptr<recob::MCSFitResultGS> mcsResultIF3D = MCSresultsIF3D.at(trkPtr.key());
    art::Ptr<recob::MCSFitResultGS> mcsResultIF2DI1 = MCSresultsIF2DI1.at(trkPtr.key());
    art::Ptr<recob::MCSFitResultGS> mcsResultIF2DI2 = MCSresultsIF2DI2.at(trkPtr.key());
    art::Ptr<recob::MCSFitResultGS> mcsResultIF2DC = MCSresultsIF2DC.at(trkPtr.key());
    art::Ptr<recob::MCSFitResultGS> mcsResultII3D = MCSresultsII3D.at(trkPtr.key());
    art::Ptr<recob::MCSFitResultGS> mcsResultII2DI1 = MCSresultsII2DI1.at(trkPtr.key());
    art::Ptr<recob::MCSFitResultGS> mcsResultII2DI2 = MCSresultsII2DI2.at(trkPtr.key());
    art::Ptr<recob::MCSFitResultGS> mcsResultII2DC = MCSresultsII2DC.at(trkPtr.key());

    art::Ptr<sbn::RangeP> rangeP = rangePs.at(trkPtr.key());
    art::FindManyP<recob::SpacePoint> fmtrkHitSPs(trkHits, e, fPFPproducer);
    std::vector<const recob::TrackHitMeta*> emptyTHMVector;
    const std::vector<const recob::TrackHitMeta*> &trkHitMetas = fmtrkHits.isValid() ? fmtrkHits.data(trkPtr.key()) : emptyTHMVector;
    art::Ptr<recob::SpacePoint> nullSP;
    std::vector<art::Ptr<recob::SpacePoint>> trkHitSPs;

    if (fmtrkHitSPs.isValid()) {
      for (unsigned i_hit = 0; i_hit < trkHits.size(); i_hit++) {
        const std::vector<art::Ptr<recob::SpacePoint>> &h_sp = fmtrkHitSPs.at(i_hit);
        if (h_sp.size()) {
          trkHitSPs.push_back(h_sp.at(0));
        }
        else {
          trkHitSPs.push_back(nullSP);
        }
      }
    }

    std::cout << "track with ID = " << trkPtr->ID() << std::endl;

    int whicht0 = -1;
    float t0 = std::numeric_limits<float>::signaling_NaN();
    std::cout << "fmT0.size() = " << fmT0.size() << std::endl;
    for (unsigned i_t0 = 0; i_t0 < fmT0.size(); i_t0++) {
      if (fmT0[i_t0].isValid() && fmT0[i_t0].at(p_pfp.key()).size()) {
        t0 = fmT0[i_t0].at(p_pfp.key()).at(0)->Time();
        //whicht0 = i_t0;
        break;
      }
    }

    int whichCRT = -1;
    float confidence = -1;
    float crtt0 = -1;
    float CRTT0 = std::numeric_limits<float>::signaling_NaN();
    std::cout << "fmCRTTaggedT0.size() = " << fmCRTTaggedT0.size() << std::endl;
    for (unsigned i_t0 = 0; i_t0 < fmCRTTaggedT0.size(); i_t0++) {
	    if (fmCRTTaggedT0[i_t0].isValid() && fmCRTTaggedT0[i_t0].at(p_pfp.key()).size()) {
        if (fmCRTTaggedT0[i_t0].at(p_pfp.key()).at(0)->TriggerBits() != 0) continue;
        if (fmCRTTaggedT0[i_t0].at(p_pfp.key()).at(0)->TriggerConfidence() > 100.) continue;
        CRTT0 = fmCRTTaggedT0[i_t0].at(p_pfp.key()).at(0)->Time();
        isCRTTagged = true;
        whichCRT = fmCRTTaggedT0[i_t0].at(p_pfp.key()).at(0)->TriggerBits();
        confidence = fmCRTTaggedT0[i_t0].at(p_pfp.key()).at(0)->TriggerConfidence();
        crtt0 = CRTT0;
        whicht0 = i_t0;
        break;
      }
    }

    std::cout << " " << std::endl;
    //std::cout << " simchannels.size " << simchannels.size() << " iscrttagged " << isCRTTagged << " CRTT0" << CRTT0 << std::endl;
    if ((!simchannels.size() && !isCRTTagged)) {
      crtt0 = -1;
      continue; }

    std::cout << "processing new track with ID = " << trkPtr->ID() << " and time = " << crtt0 << std::endl;
    std::cout << "whicht0 = " << whicht0 << std::endl;

    // Reset the track object
    *fTrack = sbn::TrackInfo();

    // Reset other persistent info
    fSnippetCount.clear();
    fWiresToSave.clear();

    // Fill the track!
    FillTrack(
      *trkPtr, pfp, t0, trkHits, trkHitMetas, trkHitSPs, calo, 
      mcsResultUF, 
      mcsResultUI, 
      mcsResultIF3D,
      mcsResultIF2DI1,
      mcsResultIF2DI2,
      mcsResultIF2DC,
      mcsResultII3D,
      mcsResultII2DI1,
      mcsResultII2DI2,
      mcsResultII2DC,
      rangeP, rawdigits, track_infos, geometry, clock_data, bt, det, 
      whichCRT, crtt0, confidence);
    fTrack->whicht0 = whicht0;
    FillTrackDaughterRays(*trkPtr, pfp, PFParticleList, PFParticleSPs);

    if (fFillTrackEndHits) FillTrackEndHits(geometry, dprop, *trkPtr, allHits, allHitSPs);

    // Fill the truth information if configured
    if (simchannels.size()) FillTrackTruth(clock_data, trkHits, mcparticles, droppedmcparticles ,AVs, TPCVols, id_to_ide_map, id_to_truehit_map, dprop, geometry); 
    if (simchannels.size()) FillDroppedTrackTruth(clock_data, trkHits, mcparticles, droppedmcparticles, AVs, TPCVols, id_to_ide_map, id_to_truehit_map, dprop, geometry); 
    bool select = false;
    if (!fSelectionTools.size()) select = true;
    if (!fSelectionTools.size()) std::cout << "selection tools not empty" << std::endl;

    // Take the OR of each selection tool
    int i_select = 0;
    for (const std::unique_ptr<sbn::ITCSSelectionTool> &t: fSelectionTools) {
      if (t->DoSelect(*fTrack)) {
        select = true;
        fTrack->selected = i_select;
        fTrack->nprescale = t->GetPrescale();
        break;
      }
      i_select ++;
    }
    select = true;
    if (select) {
      std::cout << "i_select = " << i_select << std::endl;
      fTree->Fill();
      std::cout << " after filling " << std::endl;
    }
  }
}

//returns the minimum hit time for hits in either TPC E (TPCE == true) or TPC W (TPCE == false)
float HitMinTime(
  const std::vector<sbn::TrackHitInfo> &hits, 
  bool TPCE, 
  sbn::EDet det) {
  double min = -1;
  bool hit_is_TPCE = -1;

  for (const sbn::TrackHitInfo &h: hits) {
    
    // In ICARUS, TPC E is 0, 1 and TPC W is 2, 3
    if(det == sbn::kICARUS) hit_is_TPCE = h.h.tpc <= 1;
    
    // In SBND, TPC 0 and 1
    if(det == sbn::kSBND) hit_is_TPCE = h.h.tpc <= 0;
    
    if (h.oncalo && hit_is_TPCE == TPCE) {
      if (min < 0. || h.h.time < min) min = h.h.time;
    } 
  }
  return min;
}

//returns the maximum hit time for hits in either TPC E (TPCE == true) or TPC W (TPCE == false)
float HitMaxTime(
  const std::vector<sbn::TrackHitInfo> &hits, 
  bool TPCE, 
  sbn::EDet det) {
  double max = -1;
  bool hit_is_TPCE = -1;

  for (const sbn::TrackHitInfo &h: hits) {
    
    // In ICARUS, TPC E is 0, 1 and TPC W is 2, 3
    if(det == sbn::kICARUS) hit_is_TPCE = h.h.tpc <= 1;
    
    // In SBND, TPC 0 and 1
    if(det == sbn::kSBND) hit_is_TPCE = h.h.tpc <= 0;
    
    if (h.oncalo && hit_is_TPCE == TPCE) {
      if (max < 0. || h.h.time > max) max = h.h.time;
    } 
  }
  return max;
}

sbn::Vector3D ConvertTVector(
  const TVector3 &tv) {
  sbn::Vector3D v;
  v.x = tv.X();
  v.y = tv.Y();
  v.z = tv.Z();
  return v;
}

//turn a particle position to a space-charge induced position
geo::Point_t TrajectoryToWirePosition(
  const geo::Point_t &loc, 
  const geo::TPCID &tpc) {
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  art::ServiceHandle<geo::Geometry const> geom;

  geo::Point_t ret = loc;

  // Returned X is the drift -- multiply by the drift direction to undo this
  int corr = geom->TPC(tpc).DriftDir().X();
  
  if (sce && sce->EnableSimSpatialSCE()) {
    geo::Vector_t offset = sce->GetPosOffsets(ret);
  
    ret.SetX(ret.X() + corr * offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
}

//turn a space-charge induced position to a trajectory Position
geo::Point_t WireToTrajectoryPosition(
  const geo::Point_t &loc, 
  const geo::TPCID &tpc) {
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Point_t ret = loc;

  if (sce && sce->EnableSimSpatialSCE()) {
    geo::Vector_t offset = sce->GetCalPosOffsets(ret, tpc.TPC);

    ret.SetX(ret.X() + offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
  
}

//collect MCParticle information
sbn::TrueParticle TrueParticleInfo(
  const simb::MCParticle &particle, 
  const std::vector<geo::BoxBoundedGeo> &active_volumes, 
  const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
  const std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE *>>> &id_to_ide_map,
  const std::map<int, std::vector<art::Ptr<recob::Hit>>> &id_to_truehit_map, 
  const detinfo::DetectorPropertiesData &dprop,
  const geo::GeometryCore *geo) {

  std::vector<std::pair<geo::WireID, const sim::IDE *>> empty;
  const std::vector<std::pair<geo::WireID, const sim::IDE *>> &particle_ides = id_to_ide_map.count(particle.TrackId()) ? id_to_ide_map.at(particle.TrackId()) : empty;
  
  std::vector<art::Ptr<recob::Hit>> emptyHits;
  const std::vector<art::Ptr<recob::Hit>> &particle_hits = id_to_truehit_map.count(particle.TrackId()) ? id_to_truehit_map.at(particle.TrackId()) : emptyHits;

  sbn::TrueParticle trueparticle;

  trueparticle.length = 0.;
  trueparticle.crosses_tpc = false;
  trueparticle.wallin = (int)caf::kWallNone;
  trueparticle.wallout = (int)caf::kWallNone;
  trueparticle.plane0VisE = 0.;
  trueparticle.plane1VisE = 0.;
  trueparticle.plane2VisE = 0.;
  trueparticle.plane0nhit = 0;
  trueparticle.plane1nhit = 0;
  trueparticle.plane2nhit = 0;
  for (auto const &ide_pair: particle_ides) {
    const geo::WireID &w = ide_pair.first;
    const sim::IDE *ide = ide_pair.second;
    
    if (w.Plane == 0) {
      trueparticle.plane0VisE += ide->energy / 1000. /* MeV -> GeV*/;
    }
    else if (w.Plane == 1) {
      trueparticle.plane1VisE += ide->energy / 1000. /* MeV -> GeV*/;
    }
    else if (w.Plane == 2) {
      trueparticle.plane2VisE += ide->energy / 1000. /* MeV -> GeV*/;
    }
  }

  for (const art::Ptr<recob::Hit> h: particle_hits) {
    const geo::WireID &w = h->WireID();

    if (w.Plane == 0) {
      trueparticle.plane0nhit ++;
    }
    else if (w.Plane == 1) {
      trueparticle.plane1nhit ++;
    }
    else if (w.Plane == 2) {
      trueparticle.plane2nhit ++;
    }

  }

  // if no trajectory points, then assume outside AV
  trueparticle.cont_tpc = particle.NumberTrajectoryPoints() > 0;
  trueparticle.contained = particle.NumberTrajectoryPoints() > 0;

  // Get the entry and exit points
  int entry_point = -1;
  
  int cryostat_index = -1;
  int tpc_index = -1;

  for (unsigned j = 0; j < particle.NumberTrajectoryPoints(); j++) {
    for (unsigned i = 0; i < active_volumes.size(); i++) {
      if (active_volumes.at(i).ContainsPosition(particle.Position(j).Vect())) {
        entry_point = j;
        cryostat_index = i;
        break;
      }
    }
    if (entry_point != -1) break;
  }

  // get the wall
  if (entry_point > 0) {
    trueparticle.wallin = (int)caf::GetWallCross(active_volumes.at(cryostat_index), particle.Position(entry_point).Vect(), particle.Position(entry_point-1).Vect());
  }

  int exit_point = -1;

  // now setup the cryostat the particle is in
  std::vector<geo::BoxBoundedGeo> volumes;
  if (entry_point >= 0) {
    volumes = tpc_volumes.at(cryostat_index);
    for (unsigned i = 0; i < volumes.size(); i++) {
      if (volumes[i].ContainsPosition(particle.Position(entry_point).Vect())) {
        tpc_index = i;
        trueparticle.cont_tpc = entry_point == 0;
        break;
      }
    }
    trueparticle.contained = entry_point == 0;
  }
  // if we couldn't find the initial point, set not contained
  else {
    trueparticle.contained = false;
  }

  if (tpc_index < 0) {
    trueparticle.cont_tpc = false;
  }

  // Get the length and determine if any point leaves the active volume
  // Use every trajectory point if possible
  if (entry_point >= 0) {
    // particle trajectory
    const simb::MCTrajectory &trajectory = particle.Trajectory();
    TVector3 pos = trajectory.Position(entry_point).Vect();
    for (unsigned i = entry_point+1; i < particle.NumberTrajectoryPoints(); i++) {
      TVector3 this_point = trajectory.Position(i).Vect();
      // get the exit point
      // update if particle is contained
      // check if particle has crossed TPC
      if (!trueparticle.crosses_tpc) {
        for (unsigned j = 0; j < volumes.size(); j++) {
          if (volumes[j].ContainsPosition(this_point) && tpc_index >= 0 && j != ((unsigned)tpc_index)) {
            trueparticle.crosses_tpc = true;
            break;
          }
        }
      }
      // check if particle has left tpc
      if (trueparticle.cont_tpc) {
        trueparticle.cont_tpc = volumes[tpc_index].ContainsPosition(this_point);
      }

      if (trueparticle.contained) {
        trueparticle.contained = active_volumes.at(cryostat_index).ContainsPosition(this_point);
      }

     trueparticle.length += (this_point - pos).Mag();

      if (!active_volumes.at(cryostat_index).ContainsPosition(this_point) && active_volumes.at(cryostat_index).ContainsPosition(pos)) {
        exit_point = i-1;
      }

      pos = trajectory.Position(i).Vect();
    }
  }
  if (exit_point < 0 && entry_point >= 0) {
    exit_point = particle.NumberTrajectoryPoints() - 1;
  }
  if (exit_point >= 0 && ((unsigned)exit_point) < particle.NumberTrajectoryPoints() - 1) {
    trueparticle.wallout = (int)caf::GetWallCross(active_volumes.at(cryostat_index), particle.Position(exit_point).Vect(), particle.Position(exit_point+1).Vect());
  }

  // other truth information
  trueparticle.pdg = particle.PdgCode();
  if (particle.PdgCode() == 11) std::cout << " truth electron " << std::endl;
  
  trueparticle.gen = ConvertTVector(particle.NumberTrajectoryPoints() ? particle.Position().Vect() : TVector3(-9999, -9999, -9999));
  trueparticle.genT = particle.NumberTrajectoryPoints() ? particle.Position().T() / 1000. /* ns -> us*/: -9999;
  trueparticle.genp = ConvertTVector(particle.NumberTrajectoryPoints() ? particle.Momentum().Vect(): TVector3(-9999, -9999, -9999));
  trueparticle.genE = particle.NumberTrajectoryPoints() ? particle.Momentum().E(): -9999;
  
  trueparticle.start = ConvertTVector((entry_point >= 0) ? particle.Position(entry_point).Vect(): TVector3(-9999, -9999, -9999));
  trueparticle.startT = (entry_point >= 0) ? particle.Position(entry_point).T() / 1000. /* ns-> us*/: -9999;
  trueparticle.end = ConvertTVector((exit_point >= 0) ? particle.Position(exit_point).Vect(): TVector3(-9999, -9999, -9999));
  trueparticle.endT = (exit_point >= 0) ? particle.Position(exit_point).T() / 1000. /* ns -> us */ : -9999;
  
  trueparticle.startp = ConvertTVector((entry_point >= 0) ? particle.Momentum(entry_point).Vect() : TVector3(-9999, -9999, -9999));
  trueparticle.startE = (entry_point >= 0) ? particle.Momentum(entry_point).E() : -9999.;
  trueparticle.endp = ConvertTVector((exit_point >= 0) ? particle.Momentum(exit_point).Vect() : TVector3(-9999, -9999, -9999));
  trueparticle.endE = (exit_point >= 0) ? particle.Momentum(exit_point).E() : -9999.;
  
  trueparticle.start_process = (int)caf::GetG4ProcessID(particle.Process());
  trueparticle.end_process = (int)caf::GetG4ProcessID(particle.EndProcess());
  if (particle.PdgCode() == 11)  std::cout << " truth process " << trueparticle.start_process << std::endl;
  if (particle.PdgCode() == 11) std::cout << " truth genE " << trueparticle.genE << std::endl;
  
  trueparticle.G4ID = particle.TrackId();
  trueparticle.parent = particle.Mother();

  // Organize deposition info into per-wire true "Hits" -- key is the Channel Number
  std::map<unsigned, sbn::TrueHit> truehits; 

  for (auto const &ide_pair: particle_ides) {
    const geo::WireID &w = ide_pair.first;
    unsigned c = geo->PlaneWireToChannel(w);
    const sim::IDE *ide = ide_pair.second;

    // Set stuff
    truehits[c].cryo = w.Cryostat;
    truehits[c].tpc = w.TPC;
    truehits[c].plane = w.Plane;
    truehits[c].wire = w.Wire;
    truehits[c].channel = c;

    // Average stuff using charge-weighting
    float old_elec = truehits[c].nelec;
    float new_elec = old_elec + ide->numElectrons;
    truehits[c].p.x = (truehits[c].p.x*old_elec + ide->x*ide->numElectrons) / new_elec;
    truehits[c].p.y = (truehits[c].p.y*old_elec + ide->y*ide->numElectrons) / new_elec;
    truehits[c].p.z = (truehits[c].p.z*old_elec + ide->z*ide->numElectrons) / new_elec;

    // Also get the position with space charge un-done
    geo::Point_t ide_p(ide->x, ide->y, ide->z);
    geo::Point_t ide_p_scecorr = WireToTrajectoryPosition(ide_p, w);

    truehits[c].p_scecorr.x = (truehits[c].p_scecorr.x*old_elec + ide_p_scecorr.x()*ide->numElectrons) / new_elec; 
    truehits[c].p_scecorr.y = (truehits[c].p_scecorr.y*old_elec + ide_p_scecorr.y()*ide->numElectrons) / new_elec; 
    truehits[c].p_scecorr.z = (truehits[c].p_scecorr.z*old_elec + ide_p_scecorr.z()*ide->numElectrons) / new_elec; 
    
    // Sum stuff
    truehits[c].nelec += ide->numElectrons;
    truehits[c].e += ide->energy;
    truehits[c].ndep += 1;
  }

  // Compute widths
  for (auto const &ide_pair: particle_ides) {
    const geo::WireID &w = ide_pair.first;
    unsigned c = geo->PlaneWireToChannel(w);
    const sim::IDE *ide = ide_pair.second;

    geo::Point_t ide_p(ide->x, ide->y, ide->z);
    geo::Point_t ide_p_scecorr = WireToTrajectoryPosition(ide_p, w);

    // Average stuff using charge-weighting
    float this_elec = ide->numElectrons;

    truehits[c].p_width.x += (ide_p.x() - truehits[c].p.x) * (ide_p.x() - truehits[c].p.x) * this_elec / truehits[c].nelec;
    truehits[c].p_width.y += (ide_p.y() - truehits[c].p.y) * (ide_p.y() - truehits[c].p.y) * this_elec / truehits[c].nelec;
    truehits[c].p_width.z += (ide_p.z() - truehits[c].p.z) * (ide_p.z() - truehits[c].p.z) * this_elec / truehits[c].nelec;

    truehits[c].p_scecorr_width.x += (ide_p_scecorr.x() - truehits[c].p_scecorr.x) * (ide_p_scecorr.x() - truehits[c].p_scecorr.x) * this_elec / truehits[c].nelec;
    truehits[c].p_scecorr_width.y += (ide_p_scecorr.y() - truehits[c].p_scecorr.y) * (ide_p_scecorr.y() - truehits[c].p_scecorr.y) * this_elec / truehits[c].nelec;
    truehits[c].p_scecorr_width.z += (ide_p_scecorr.z() - truehits[c].p_scecorr.z) * (ide_p_scecorr.z() - truehits[c].p_scecorr.z) * this_elec / truehits[c].nelec;
  }

  // Convert to vector
  std::vector<sbn::TrueHit> truehits_v;
  for (auto const &p: truehits) {
    truehits_v.push_back(p.second);
  }

  // Compute the time of each hit
  for (sbn::TrueHit &h: truehits_v) {
    h.time = dprop.ConvertXToTicks(h.p.x, h.plane, h.tpc, h.cryo);

    double xdrift = abs(h.p.x - geo->Plane(geo::PlaneID(h.cryo, h.tpc, 0)).GetCenter().X());
    h.tdrift = xdrift / dprop.DriftVelocity(); 
  }

  // Compute the pitch of each hit and order it in the trajectory
  for (sbn::TrueHit &h: truehits_v) {
    // Use the SCE-undone hit since this matches to the Trajectory
    TVector3 h_p(h.p_scecorr.x, h.p_scecorr.y, h.p_scecorr.z);

    TVector3 direction;
    float closest_dist = -1.;
    int traj_index = -1;
    for (unsigned i_traj = 0; i_traj < particle.NumberTrajectoryPoints(); i_traj++) {
      if (closest_dist < 0. || (particle.Position(i_traj).Vect() - h_p).Mag() < closest_dist) {
        direction = particle.Momentum(i_traj).Vect().Unit();
        closest_dist = (particle.Position(i_traj).Vect() - h_p).Mag();
        traj_index = i_traj;
      }
    }

    // If we got a direction, get the pitch
    if (closest_dist >= 0. && direction.Mag() > 1e-4) {
      geo::PlaneID plane(h.cryo, h.tpc, h.plane);
      float angletovert = geo->WireAngleToVertical(geo->View(plane), plane) - 0.5*::util::pi<>();
      float cosgamma = abs(cos(angletovert) * direction.Z() + sin(angletovert) * direction.Y());
      float pitch = geo->WirePitch(plane) / cosgamma;
      h.pitch = pitch;
    }
    else {
      h.pitch = -1.;
    }
    // And the pitch induced by SCE
    if (closest_dist >= 0. && direction.Mag() > 1e-4) {
      geo::PlaneID plane(h.cryo, h.tpc, h.plane);
      float angletovert = geo->WireAngleToVertical(geo->View(plane), plane) - 0.5*::util::pi<>();

      TVector3 loc_mdx_v = h_p - direction * (geo->WirePitch(geo->View(plane)) / 2.);
      TVector3 loc_pdx_v = h_p + direction * (geo->WirePitch(geo->View(plane)) / 2.);

      // Convert types for helper functions
      geo::Point_t loc_mdx(loc_mdx_v.X(), loc_mdx_v.Y(), loc_mdx_v.Z());
      geo::Point_t loc_pdx(loc_pdx_v.X(), loc_pdx_v.Y(), loc_pdx_v.Z());
      geo::Point_t h_p_point(h_p.X(), h_p.Y(), h_p.Z());

      loc_mdx = TrajectoryToWirePosition(loc_mdx, plane);
      loc_pdx = TrajectoryToWirePosition(loc_pdx, plane);
      
      // Direction at wires
      geo::Vector_t dir = (loc_pdx - loc_mdx) /  (loc_mdx - loc_pdx).r(); 

      // Pitch at wires
      double cosgamma = std::abs(std::sin(angletovert)*dir.Y() + std::cos(angletovert)*dir.Z());
      double pitch;
      if (cosgamma) {
        pitch = geo->WirePitch(geo->View(plane))/cosgamma;
      }
      else {
        pitch = 0.;
      }

      // Now bring that back to the particle trajectory
      geo::Point_t loc_w = TrajectoryToWirePosition(h_p_point, plane);
      
      geo::Point_t locw_pdx_traj = WireToTrajectoryPosition(loc_w + pitch*dir, plane);
      geo::Point_t loc = WireToTrajectoryPosition(loc_w, plane);
      
      h.pitch_sce = (locw_pdx_traj - loc).R();
    }
    else {
      h.pitch_sce = -1.;
    }

    // And the trajectory location
    h.itraj = traj_index;

    // And the residual range of the hit
    h.rr = 0.;
    if (traj_index >= 0) {
      for (int i_traj = traj_index+1; i_traj < (int)particle.NumberTrajectoryPoints(); i_traj++) {
        h.rr += (particle.Position(i_traj).Vect() - particle.Position(i_traj-1).Vect()).Mag();
      }

      // Also account for the distance from the Hit point to the matched trajectory point
      double hit_distance_along_particle = (h_p - particle.Position(traj_index).Vect()).Dot(particle.Momentum(traj_index).Vect().Unit());
      h.rr += -hit_distance_along_particle;

    }
  }

  // Order the hits by their location along the trajectory, start to end
  std::sort(truehits_v.begin(), truehits_v.end(), 
    [](auto const &lhs, auto const &rhs) {
      return lhs.itraj < rhs.itraj;
  });

  // Save depositions into the True Particle
  for (sbn::TrueHit &h: truehits_v) {
    if (h.plane == 0) {
      trueparticle.truehits0.push_back(h);
    }
    else if (h.plane == 1) {
      trueparticle.truehits1.push_back(h);
    }
    else if (h.plane == 2) {
      trueparticle.truehits2.push_back(h);
    }
  }

  // Save the true trajectory
  for (unsigned i_traj = 0; i_traj < particle.NumberTrajectoryPoints(); i_traj++) {
    // Get trajectory point
    TVector3 traj = particle.Position(i_traj).Vect();
    geo::Point_t traj_p(traj.X(), traj.Y(), traj.Z());

    // lookup TPC
    geo::TPCID tpc; // invalid by default
    for (auto const &cryo: geo->Iterate<geo::CryostatGeo>()) {
      for (auto const& TPC : geo->Iterate<geo::TPCGeo>(cryo.ID())) {
        if (TPC.ActiveBoundingBox().ContainsPosition(traj_p)) {
          tpc = TPC.ID();
          break;
        }
      }
      if (tpc.isValid) break;
    }

    // add in space-charge-deflected position if applicable
    geo::Point_t traj_p_sce = tpc.isValid ? TrajectoryToWirePosition(traj_p, tpc) : traj_p;

    sbn::Vector3D traj_v;
    traj_v.x = traj_p.x();
    traj_v.y = traj_p.y();
    traj_v.z = traj_p.z();

    sbn::Vector3D traj_v_sce;
    traj_v_sce.x = traj_p_sce.x();
    traj_v_sce.y = traj_p_sce.y();
    traj_v_sce.z = traj_p_sce.z();

    trueparticle.traj.push_back(traj_v);
    trueparticle.traj_sce.push_back(traj_v_sce);
  }

  return trueparticle;
}

void sbn::TrackCaloSkimmer::FillTrackEndHits(
  const geo::GeometryCore *geometry,
  const detinfo::DetectorPropertiesData &dprop,
  const recob::Track &track,
  const std::vector<art::Ptr<recob::Hit>> &allHits,
  const art::FindManyP<recob::SpacePoint> &allHitSPs) {

  (void) dprop; // TODO: use??

  geo::TPCID tpc_end = geometry->FindTPCAtPosition(track.End());
  if (!tpc_end) return;

  geo::PlaneID plane_end(tpc_end, 2 /* collection */);

  float end_w = geometry->WireCoordinate(track.End(), plane_end);

  float end_t = -1000.;
  float closest_wire_dist = -1.;
  // Get the hit closest to the end to get the end time 
  for (const TrackHitInfo &h: fTrack->hits2) {
    if (h.oncalo && (closest_wire_dist < 0. || abs(h.h.wire - end_w) < closest_wire_dist)) {
      closest_wire_dist = abs(h.h.wire - end_w);
      end_t = h.h.time;
    }
  }

  for (const art::Ptr<recob::Hit> &hit: allHits) {
    geo::PlaneID h_p = hit->WireID();
    if (h_p != plane_end) continue;

    // Inside the box?
    float h_w = (float)hit->WireID().Wire;
    float h_t = hit->PeakTime();

    if (abs(h_w - end_w) < fTrackEndHitWireBox &&
        abs(h_t - end_t) < fTrackEndHitTimeBox) {

      // HitInfo to save
      sbn::HitInfo hinfo;
      
      // information from the hit object
      hinfo.integral = hit->Integral();
      hinfo.sumadc = hit->SummedADC();
      hinfo.width = hit->RMS();
      hinfo.time = hit->PeakTime();
      hinfo.mult = hit->Multiplicity();
      hinfo.wire = hit->WireID().Wire;
      hinfo.plane = hit->WireID().Plane;
      hinfo.tpc = hit->WireID().TPC;
      hinfo.end = hit->EndTick();
      hinfo.start = hit->StartTick();
      hinfo.id = hit.key();

      const std::vector<art::Ptr<recob::SpacePoint>> &h_sp = allHitSPs.at(hit.key());
      if (h_sp.size()) {
        const recob::SpacePoint &sp = *h_sp[0];
        hinfo.sp.x = sp.position().x();
        hinfo.sp.y = sp.position().y();
        hinfo.sp.z = sp.position().z();

        hinfo.hasSP = true;
      }
      else {
        hinfo.hasSP = false;
      }

      fTrack->endhits.push_back(hinfo);
      
    } 
  }

}

void sbn::TrackCaloSkimmer::FillTrackTruth(
  const detinfo::DetectorClocksData &clock_data,
  const std::vector<art::Ptr<recob::Hit>> &trkHits,
  const std::vector<art::Ptr<simb::MCParticle>> &mcparticles,
  const std::vector<art::Ptr<simb::MCParticle>> &droppedmcparticles,
  const std::vector<geo::BoxBoundedGeo> &active_volumes,
  const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
  const std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map,
  const std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map,
  const detinfo::DetectorPropertiesData &dprop,
  const geo::GeometryCore *geo) {

  // Lookup the true-particle match -- use utils in CAF
  std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, trkHits, true);
  float total_energy = CAFRecoUtils::TotalHitEnergy(clock_data, trkHits); 
  fTrack->truth.depE = total_energy / 1000. /* MeV -> GeV */;
  
  // sort highest energy match to lowest
  std::sort(matches.begin(), matches.end(),
      [](const auto &a, const auto &b) {
        return a.second > b.second;
      }
  );

  // Save the best match
  if (matches.size()) {
    std::pair<int, float> bestmatch = matches[0];

    fTrack->truth.pur = bestmatch.second / total_energy;

    std::cout << " mcparticles size " << mcparticles.size() << std::endl;
    std::cout << " dropped size " << droppedmcparticles.size() << std::endl;

    for (const art::Ptr<simb::MCParticle> &p_mcp: mcparticles) {
      //std::cout << " nondropped trackid " << p_mcp->TrackId() << " pdg: " << p_mcp->PdgCode() << " process: " << p_mcp->Process() << std::endl;
      //std::cout << " nondropped trackid " << p_mcp->TrackId() << " energy: " << p_mcp->E()  << std::endl;
      if (p_mcp->TrackId() == bestmatch.first) {
        std::cout << "Matched! Track ID: " << p_mcp->TrackId() << " pdg: " << p_mcp->PdgCode() << " process: " << p_mcp->EndProcess() << std::endl;
        fTrack->truth.p = TrueParticleInfo(*p_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        fTrack->truth.eff = fTrack->truth.depE / (fTrack->truth.p.plane0VisE + fTrack->truth.p.plane1VisE + fTrack->truth.p.plane2VisE);

        // Lookup any Michel
        for (const art::Ptr<simb::MCParticle> &d_mcp: mcparticles) {
          if (d_mcp->Mother() == p_mcp->TrackId() && // correct parent
              (d_mcp->Process() == "Decay" || d_mcp->Process() == "muMinusCaptureAtRest") && // correct process
              abs(d_mcp->PdgCode()) == 11) { // correct PDG code

            fTrack->truth.michel = TrueParticleInfo(*d_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
            break;
          }
        }
        std::vector<std::pair<const simb::MCParticle*, float>> ionization_electrons;

        // Now look for ionization electrons from this particle in the dropped list
        for (const art::Ptr<simb::MCParticle> &c_mcp: droppedmcparticles) {
          
          if (!(c_mcp->PdgCode() == 11 && c_mcp->Process() == "muIoni")) continue;
          if (c_mcp->Mother() == p_mcp->TrackId()) {
            
            std::cout << "Dropped Matched! Track ID: " << c_mcp->TrackId() << " pdg: " << c_mcp->PdgCode() << " process: " << c_mcp->Process() << "  Energy: " << c_mcp->E() * 1000.0 << " keV" << std::endl; //da vedere se MeV o kev
            
            float energy = c_mcp->E(); // in GeV
            ionization_electrons.emplace_back(c_mcp.get(), energy);

            //sbn::TrueParticle deltaray = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
            
            //fTrack->truth.deltarays.push_back(deltaray);
          
          }
        }

        // Sort ionization electrons by energy
        std::sort(ionization_electrons.begin(), ionization_electrons.end(),
          [](const auto &a, const auto &b) {
            return a.second > b.second; // descending order
          });

        // Fill the deltaray information
        //fTrack->truth.Ndeltarays = ionization_electrons.size();
          
        if (ionization_electrons.size() > 0) {
          const simb::MCParticle* c_mcp = ionization_electrons[0].first;
          fTrack->truth.deltaray_1 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }
        
        if (ionization_electrons.size() > 1) {
          const simb::MCParticle* c_mcp = ionization_electrons[1].first;
          fTrack->truth.deltaray_2 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }
        
        if (ionization_electrons.size() > 2) {
          const simb::MCParticle* c_mcp = ionization_electrons[2].first;
          fTrack->truth.deltaray_3 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }

        if (ionization_electrons.size() > 3) {
          const simb::MCParticle* c_mcp = ionization_electrons[3].first;
          fTrack->truth.deltaray_4 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }

        if (ionization_electrons.size() > 4) {
          const simb::MCParticle* c_mcp = ionization_electrons[4].first;
          fTrack->truth.deltaray_5 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }
        
        if (ionization_electrons.size() > 5) {
          const simb::MCParticle* c_mcp = ionization_electrons[5].first;
          fTrack->truth.deltaray_6 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }
        
        if (ionization_electrons.size() > 6) {
          const simb::MCParticle* c_mcp = ionization_electrons[6].first;
          fTrack->truth.deltaray_7 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }
        
        if (ionization_electrons.size() > 7) {
          const simb::MCParticle* c_mcp = ionization_electrons[7].first;
          fTrack->truth.deltaray_8 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }

        if (ionization_electrons.size() > 8) {
          const simb::MCParticle* c_mcp = ionization_electrons[8].first;
          fTrack->truth.deltaray_9 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }

        if (ionization_electrons.size() > 9) {
          const simb::MCParticle* c_mcp = ionization_electrons[9].first;
          fTrack->truth.deltaray_10 = TrueParticleInfo(*c_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        }

        break;
      } 
    }
  }
}

void sbn::TrackCaloSkimmer::FillTrackDaughterRays(const recob::Track &trk, const recob::PFParticle &pfp, const std::vector<art::Ptr<recob::PFParticle>> &PFParticleList, const art::FindManyP<recob::SpacePoint> &PFParticleSPs) {
  for (unsigned d: pfp.Daughters()) {
    const recob::PFParticle &d_pfp = *PFParticleList[d];

    fTrack->daughter_pdg.push_back(d_pfp.PdgCode());

    unsigned nsp = 0;
    float min_distance = -1.;
    for (const art::Ptr<recob::SpacePoint> &sp: PFParticleSPs.at(d)) {
      if (min_distance < 0. || (sp->position() - trk.End()).r() < min_distance) {
        min_distance = (sp->position() - trk.End()).r();
      }
      nsp++;
    }

    fTrack->daughter_sp_toend_dist.push_back(min_distance);
    fTrack->daughter_nsp.push_back(nsp);
  }

}

void sbn::TrackCaloSkimmer::FillTrack(const recob::Track &track, 
  const recob::PFParticle &pfp, float t0, 
  const std::vector<art::Ptr<recob::Hit>> &hits,
  const std::vector<const recob::TrackHitMeta*> &thms,
  const std::vector<art::Ptr<recob::SpacePoint>> &sps,
  const std::vector<art::Ptr<anab::Calorimetry>> &calo,
	const art::Ptr<recob::MCSFitResult> &mcsUF,
  const art::Ptr<recob::MCSFitResult> &mcsUI,
	const art::Ptr<recob::MCSFitResultGS> &mcsIF3D,
  const art::Ptr<recob::MCSFitResultGS> &mcsIF2DI1,
  const art::Ptr<recob::MCSFitResultGS> &mcsIF2DI2,
  const art::Ptr<recob::MCSFitResultGS> &mcsIF2DC,
	const art::Ptr<recob::MCSFitResultGS> &mcsII3D,
  const art::Ptr<recob::MCSFitResultGS> &mcsII2DI1,
  const art::Ptr<recob::MCSFitResultGS> &mcsII2DI2,
  const art::Ptr<recob::MCSFitResultGS> &mcsII2DC,
  const art::Ptr<sbn::RangeP> &rangeP,
  const std::map<geo::WireID, art::Ptr<raw::RawDigit>> &rawdigits,
  const std::vector<GlobalTrackInfo> &tracks,
  const geo::GeometryCore *geo,
  const detinfo::DetectorClocksData &clock_data,
  const cheat::BackTrackerService *bt_serv,
  const sbn::EDet det,
  const int whichcrt,
  const float crtt0,
  const float confidence) {
  std::cout << " filltrack " << std::endl;

  fTrack->whichCRT = whichcrt;
  fTrack->CRT_T0 = crtt0;
  fTrack->CRTconfidence = confidence;

  // Fill top level stuff
  fTrack->meta = fMeta;
  fTrack->t0 = t0;
  fTrack->id = track.ID();
  fTrack->clear_cosmic_muon = pfp.Parent() == recob::PFParticle::kPFParticlePrimary;
  std::cout << " filltrack " << std::endl;

  fTrack->length = track.Length();
  fTrack->start.x = track.Start().X();
  fTrack->start.y = track.Start().Y();
  fTrack->start.z = track.Start().Z();
  fTrack->end.x = track.End().X();
  fTrack->end.y = track.End().Y();
  fTrack->end.z = track.End().Z();
  fTrack->dir.x = track.StartDirection().X();
  fTrack->dir.y = track.StartDirection().Y();
  fTrack->dir.z = track.StartDirection().Z();
  std::cout << " filltrack " << std::endl;

  fTrack->range_p = rangeP->range_p;
  
  fTrack->mcs_pbest_UF = mcsUF->fwdMomentum();
  fTrack->mcs_perr_UF  = mcsUF->fwdMomUncertainty();
  fTrack->mcs_pbest_UI = mcsUI->fwdMomentum();
  fTrack->mcs_perr_UI  = mcsUI->fwdMomUncertainty();
  std::cout << " filltrack " << std::endl;

  for (size_t ja = 0; ja < mcsUF->scatterAngles().size(); ja++)
    fTrack->mcs_angles_UF.push_back(mcsUF->scatterAngles().at(ja));
  for (size_t ja = 0; ja < mcsUI->scatterAngles().size(); ja++)
    fTrack->mcs_angles_UI.push_back(mcsUI->scatterAngles().at(ja));
  std::cout << " filltrack " << std::endl;

  fTrack->mcs_stop = mcsIF3D->GeoStopCheck();
  
  fTrack->mcs_isdelta_3D = mcsIF3D->IsDelta();
  fTrack->mcs_isdelta_2DI1 = mcsIF2DI1->IsDelta();
  fTrack->mcs_isdelta_2DI2 = mcsIF2DI2->IsDelta();
  fTrack->mcs_isdelta_2DC = mcsIF2DC->IsDelta();

  fTrack->mcs_isdelta_indexC = mcsIF2DC->IsDeltaIndex();
  fTrack->mcs_isdelta_indexI2 = mcsIF2DI2->IsDeltaIndex();
  fTrack->mcs_isdelta_indexI1 = mcsIF2DI1->IsDeltaIndex();
  std::cout << " filltrack " << std::endl;

  for (unsigned int jd = 0; jd < fTrack->mcs_isdelta_2DI2.size(); jd++) 
    if (fTrack->mcs_isdelta_2DI2.at(jd) > 0)
      std::cout << " candidate isdelta 2di2 " << fTrack->mcs_isdelta_2DI2.at(jd) << std::endl; 
  for (unsigned int jd = 0; jd < fTrack->mcs_isdelta_2DI1.size(); jd++) 
    if (fTrack->mcs_isdelta_2DI1.at(jd) > 0)
      std::cout << " candidate isdelta 2di1 " << fTrack->mcs_isdelta_2DI1.at(jd) << std::endl; 
  for (unsigned int jd = 0; jd < fTrack->mcs_isdelta_2DC.size(); jd++) 
    if (fTrack->mcs_isdelta_2DC.at(jd) > 0)
      std::cout << " candidate isdelta 2dC " << fTrack->mcs_isdelta_2DC.at(jd) << std::endl; 
  for (unsigned int jd = 0; jd < fTrack->mcs_isdelta_indexI2.size(); jd++) 
    if (fTrack->mcs_isdelta_indexI2.at(jd) > 0)
      std::cout << " candidate isdelta indexI2 " << fTrack->mcs_isdelta_indexI2.at(jd) << std::endl; 

  fTrack->mcs_sigma3p_3D   = mcsIF3D->sigma3P();
  fTrack->mcs_sigma3p_2DI1 = mcsIF2DI1->sigma3P();
  fTrack->mcs_sigma3p_2DI2 = mcsIF2DI2->sigma3P();
  fTrack->mcs_sigma3p_2DC  = mcsIF2DC->sigma3P();
  
  fTrack->mcs_c2atrange_IF3D = mcsIF3D->C2AtRange();
  fTrack->mcs_c2atrange_IF2DI1 = mcsIF2DI1->C2AtRange();
  fTrack->mcs_c2atrange_IF2DI2 = mcsIF2DI2->C2AtRange();
  fTrack->mcs_c2atrange_IF2DC = mcsIF2DC->C2AtRange();
  fTrack->mcs_c2atrange_II3D = mcsII3D->C2AtRange();
  fTrack->mcs_c2atrange_II2DI1 = mcsII2DI1->C2AtRange();
  fTrack->mcs_c2atrange_II2DI2 = mcsII2DI2->C2AtRange();
  fTrack->mcs_c2atrange_II2DC = mcsII2DC->C2AtRange();
  
  fTrack->mcs_tailsatrange_IF3D = mcsIF3D->TailsAtRange();
  fTrack->mcs_tailsatrange_IF2DI1 = mcsIF2DI1->TailsAtRange();
  fTrack->mcs_tailsatrange_IF2DI2 = mcsIF2DI2->TailsAtRange();
  fTrack->mcs_tailsatrange_IF2DC = mcsIF2DC->TailsAtRange();
  fTrack->mcs_tailsatrange_II3D = mcsII3D->TailsAtRange();
  fTrack->mcs_tailsatrange_II2DI1 = mcsII2DI1->TailsAtRange();
  fTrack->mcs_tailsatrange_II2DI2 = mcsII2DI2->TailsAtRange();
  fTrack->mcs_tailsatrange_II2DC = mcsII2DC->TailsAtRange();

  fTrack->mcs_l1DI1 = mcsIF2DI1->length1D();
  fTrack->mcs_l2DI1 = mcsIF2DI1->length2D();
  fTrack->mcs_l3DI1 = mcsIF2DI1->length3D();
  fTrack->mcs_l1DI2 = mcsIF2DI2->length1D();
  fTrack->mcs_l2DI2 = mcsIF2DI2->length2D();
  fTrack->mcs_l3DI2 = mcsIF2DI2->length3D();
  fTrack->mcs_l1DC  = mcsIF2DC->length1D();
  fTrack->mcs_l2DC  = mcsIF2DC->length2D();
  fTrack->mcs_l3DC  = mcsIF2DC->length3D();

  fTrack->mcs_pbest_IF3D  = mcsIF3D->bestMomentum();
  fTrack->mcs_perr_IF3D   = mcsIF3D->errorMomentum();
  fTrack->mcs_pmin_IF3D   = mcsIF3D->minMomentum();
  fTrack->mcs_pmax_IF3D   = mcsIF3D->maxMomentum();
  fTrack->mcs_alpha_IF3D  = mcsIF3D->alphaFit();
  fTrack->mcs_dalpha_IF3D = mcsIF3D->alphaFitErr();
  fTrack->mcs_beta_IF3D   = mcsIF3D->betaFit();
  fTrack->mcs_dbeta_IF3D  = mcsIF3D->betaFitErr();
  for (size_t ja = 0; ja < mcsIF3D->testMomentum().size(); ja++)
    fTrack->mcs_ptest_IF3D.push_back(mcsIF3D->testMomentum().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->C2Function().size(); ja++)
    fTrack->mcs_c2function_IF3D.push_back(mcsIF3D->C2Function().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->TailsSize().size(); ja++)
    fTrack->mcs_tailssize_IF3D.push_back(mcsIF3D->TailsSize().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->segmentLengths().size(); ja++)
    fTrack->mcs_seglens_IF3D.push_back(mcsIF3D->segmentLengths().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->segmentCumLengths().size(); ja++)
    fTrack->mcs_cumseglens_IF3D.push_back(mcsIF3D->segmentCumLengths().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->segmentHits().size(); ja++)
    fTrack->mcs_seghits_IF3D.push_back(mcsIF3D->segmentHits().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->segmentCumHits().size(); ja++)
    fTrack->mcs_cumseghits_IF3D.push_back(mcsIF3D->segmentCumHits().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->expectedLinAngles().size(); ja++)
    fTrack->mcs_dthetalinexp_IF3D.push_back(mcsIF3D->expectedLinAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->measuredLinAngles().size(); ja++)
    fTrack->mcs_dthetalin_IF3D.push_back(mcsIF3D->measuredLinAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->expectedPolyAngles().size(); ja++)
    fTrack->mcs_dthetapolyexp_IF3D.push_back(mcsIF3D->expectedPolyAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF3D->measuredPolyAngles().size(); ja++)
    fTrack->mcs_dthetapoly_IF3D.push_back(mcsIF3D->measuredPolyAngles().at(ja));

  fTrack->mcs_pbest_IF2DI1  = mcsIF2DI1->bestMomentum();
  fTrack->mcs_perr_IF2DI1   = mcsIF2DI1->errorMomentum();
  fTrack->mcs_pmin_IF2DI1   = mcsIF2DI1->minMomentum();
  fTrack->mcs_pmax_IF2DI1   = mcsIF2DI1->maxMomentum();
  fTrack->mcs_alpha_IF2DI1  = mcsIF2DI1->alphaFit();
  fTrack->mcs_dalpha_IF2DI1 = mcsIF2DI1->alphaFitErr();
  fTrack->mcs_beta_IF2DI1   = mcsIF2DI1->betaFit();
  fTrack->mcs_dbeta_IF2DI1  = mcsIF2DI1->betaFitErr();
  for (size_t ja = 0; ja < mcsIF2DI1->testMomentum().size(); ja++)
    fTrack->mcs_ptest_IF2DI1.push_back(mcsIF2DI1->testMomentum().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->C2Function().size(); ja++)
    fTrack->mcs_c2function_IF2DI1.push_back(mcsIF2DI1->C2Function().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->TailsSize().size(); ja++)
    fTrack->mcs_tailssize_IF2DI1.push_back(mcsIF2DI1->TailsSize().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->segmentLengths().size(); ja++)
    fTrack->mcs_seglens_IF2DI1.push_back(mcsIF2DI1->segmentLengths().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->segmentCumLengths().size(); ja++)
    fTrack->mcs_cumseglens_IF2DI1.push_back(mcsIF2DI1->segmentCumLengths().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->segmentHits().size(); ja++)
    fTrack->mcs_seghits_IF2DI1.push_back(mcsIF2DI1->segmentHits().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->segmentCumHits().size(); ja++)
    fTrack->mcs_cumseghits_IF2DI1.push_back(mcsIF2DI1->segmentCumHits().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->expectedLinAngles().size(); ja++)
    fTrack->mcs_dthetalinexp_IF2DI1.push_back(mcsIF2DI1->expectedLinAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->measuredLinAngles().size(); ja++)
    fTrack->mcs_dthetalin_IF2DI1.push_back(mcsIF2DI1->measuredLinAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->expectedPolyAngles().size(); ja++)
    fTrack->mcs_dthetapolyexp_IF2DI1.push_back(mcsIF2DI1->expectedPolyAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI1->measuredPolyAngles().size(); ja++)
    fTrack->mcs_dthetapoly_IF2DI1.push_back(mcsIF2DI1->measuredPolyAngles().at(ja));

  fTrack->mcs_pbest_IF2DI2  = mcsIF2DI2->bestMomentum();
  fTrack->mcs_perr_IF2DI2   = mcsIF2DI2->errorMomentum();
  fTrack->mcs_pmin_IF2DI2   = mcsIF2DI2->minMomentum();
  fTrack->mcs_pmax_IF2DI2   = mcsIF2DI2->maxMomentum();
  fTrack->mcs_alpha_IF2DI2  = mcsIF2DI2->alphaFit();
  fTrack->mcs_dalpha_IF2DI2 = mcsIF2DI2->alphaFitErr();
  fTrack->mcs_beta_IF2DI2   = mcsIF2DI2->betaFit();
  fTrack->mcs_dbeta_IF2DI2  = mcsIF2DI2->betaFitErr();
  for (size_t ja = 0; ja < mcsIF2DI2->testMomentum().size(); ja++)
    fTrack->mcs_ptest_IF2DI2.push_back(mcsIF2DI2->testMomentum().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->C2Function().size(); ja++)
    fTrack->mcs_c2function_IF2DI2.push_back(mcsIF2DI2->C2Function().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->TailsSize().size(); ja++)
    fTrack->mcs_tailssize_IF2DI2.push_back(mcsIF2DI2->TailsSize().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->segmentLengths().size(); ja++)
    fTrack->mcs_seglens_IF2DI2.push_back(mcsIF2DI2->segmentLengths().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->segmentCumLengths().size(); ja++)
    fTrack->mcs_cumseglens_IF2DI2.push_back(mcsIF2DI2->segmentCumLengths().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->segmentHits().size(); ja++)
    fTrack->mcs_seghits_IF2DI2.push_back(mcsIF2DI2->segmentHits().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->segmentCumHits().size(); ja++)
    fTrack->mcs_cumseghits_IF2DI2.push_back(mcsIF2DI2->segmentCumHits().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->expectedLinAngles().size(); ja++)
    fTrack->mcs_dthetalinexp_IF2DI2.push_back(mcsIF2DI2->expectedLinAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->measuredLinAngles().size(); ja++)
    fTrack->mcs_dthetalin_IF2DI2.push_back(mcsIF2DI2->measuredLinAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->expectedPolyAngles().size(); ja++)
    fTrack->mcs_dthetapolyexp_IF2DI2.push_back(mcsIF2DI2->expectedPolyAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DI2->measuredPolyAngles().size(); ja++)
    fTrack->mcs_dthetapoly_IF2DI2.push_back(mcsIF2DI2->measuredPolyAngles().at(ja));

  fTrack->mcs_pbest_IF2DC  = mcsIF2DC->bestMomentum();
  fTrack->mcs_perr_IF2DC   = mcsIF2DC->errorMomentum();
  fTrack->mcs_pmin_IF2DC   = mcsIF2DC->minMomentum();
  fTrack->mcs_pmax_IF2DC   = mcsIF2DC->maxMomentum();
  fTrack->mcs_alpha_IF2DC  = mcsIF2DC->alphaFit();
  fTrack->mcs_dalpha_IF2DC = mcsIF2DC->alphaFitErr();
  fTrack->mcs_beta_IF2DC   = mcsIF2DC->betaFit();
  fTrack->mcs_dbeta_IF2DC  = mcsIF2DC->betaFitErr();
  for (size_t ja = 0; ja < mcsIF2DC->testMomentum().size(); ja++)
    fTrack->mcs_ptest_IF2DC.push_back(mcsIF2DC->testMomentum().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->C2Function().size(); ja++)
    fTrack->mcs_c2function_IF2DC.push_back(mcsIF2DC->C2Function().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->TailsSize().size(); ja++)
    fTrack->mcs_tailssize_IF2DC.push_back(mcsIF2DC->TailsSize().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->segmentLengths().size(); ja++)
    fTrack->mcs_seglens_IF2DC.push_back(mcsIF2DC->segmentLengths().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->segmentCumLengths().size(); ja++)
    fTrack->mcs_cumseglens_IF2DC.push_back(mcsIF2DC->segmentCumLengths().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->segmentHits().size(); ja++)
    fTrack->mcs_seghits_IF2DC.push_back(mcsIF2DC->segmentHits().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->segmentCumHits().size(); ja++)
    fTrack->mcs_cumseghits_IF2DC.push_back(mcsIF2DC->segmentCumHits().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->expectedLinAngles().size(); ja++)
    fTrack->mcs_dthetalinexp_IF2DC.push_back(mcsIF2DC->expectedLinAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->measuredLinAngles().size(); ja++)
    fTrack->mcs_dthetalin_IF2DC.push_back(mcsIF2DC->measuredLinAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->expectedPolyAngles().size(); ja++)
    fTrack->mcs_dthetapolyexp_IF2DC.push_back(mcsIF2DC->expectedPolyAngles().at(ja));
  for (size_t ja = 0; ja < mcsIF2DC->measuredPolyAngles().size(); ja++)
    fTrack->mcs_dthetapoly_IF2DC.push_back(mcsIF2DC->measuredPolyAngles().at(ja));

  fTrack->mcs_pbest_II3D  = mcsII3D->bestMomentum();
  fTrack->mcs_perr_II3D   = mcsII3D->errorMomentum();
  fTrack->mcs_pmin_II3D   = mcsII3D->minMomentum();
  fTrack->mcs_pmax_II3D   = mcsII3D->maxMomentum();
  fTrack->mcs_alpha_II3D  = mcsII3D->alphaFit();
  fTrack->mcs_dalpha_II3D = mcsII3D->alphaFitErr();
  fTrack->mcs_beta_II3D   = mcsII3D->betaFit();
  fTrack->mcs_dbeta_II3D  = mcsII3D->betaFitErr();
  for (size_t ja = 0; ja < mcsII3D->testMomentum().size(); ja++)
    fTrack->mcs_ptest_II3D.push_back(mcsII3D->testMomentum().at(ja));
  for (size_t ja = 0; ja < mcsII3D->C2Function().size(); ja++)
    fTrack->mcs_c2function_II3D.push_back(mcsII3D->C2Function().at(ja));
  for (size_t ja = 0; ja < mcsII3D->TailsSize().size(); ja++)
    fTrack->mcs_tailssize_II3D.push_back(mcsII3D->TailsSize().at(ja));

  fTrack->mcs_pbest_II2DI1  = mcsII2DI1->bestMomentum();
  fTrack->mcs_perr_II2DI1   = mcsII2DI1->errorMomentum();
  fTrack->mcs_pmin_II2DI1   = mcsII2DI1->minMomentum();
  fTrack->mcs_pmax_II2DI1   = mcsII2DI1->maxMomentum();
  fTrack->mcs_alpha_II2DI1  = mcsII2DI1->alphaFit();
  fTrack->mcs_dalpha_II2DI1 = mcsII2DI1->alphaFitErr();
  fTrack->mcs_beta_II2DI1   = mcsII2DI1->betaFit();
  fTrack->mcs_dbeta_II2DI1  = mcsII2DI1->betaFitErr();
  for (size_t ja = 0; ja < mcsII2DI1->testMomentum().size(); ja++)
    fTrack->mcs_ptest_II2DI1.push_back(mcsII2DI1->testMomentum().at(ja));
  for (size_t ja = 0; ja < mcsII2DI1->C2Function().size(); ja++)
    fTrack->mcs_c2function_II2DI1.push_back(mcsII2DI1->C2Function().at(ja));
  for (size_t ja = 0; ja < mcsII2DI1->TailsSize().size(); ja++)
    fTrack->mcs_tailssize_II2DI1.push_back(mcsII2DI1->TailsSize().at(ja));

  fTrack->mcs_pbest_II2DI2  = mcsII2DI2->bestMomentum();
  fTrack->mcs_perr_II2DI2   = mcsII2DI2->errorMomentum();
  fTrack->mcs_pmin_II2DI2   = mcsII2DI2->minMomentum();
  fTrack->mcs_pmax_II2DI2   = mcsII2DI2->maxMomentum();
  fTrack->mcs_alpha_II2DI2  = mcsII2DI2->alphaFit();
  fTrack->mcs_dalpha_II2DI2 = mcsII2DI2->alphaFitErr();
  fTrack->mcs_beta_II2DI2   = mcsII2DI2->betaFit();
  fTrack->mcs_dbeta_II2DI2  = mcsII2DI2->betaFitErr();
  for (size_t ja = 0; ja < mcsII2DI2->testMomentum().size(); ja++)
    fTrack->mcs_ptest_II2DI2.push_back(mcsII2DI2->testMomentum().at(ja));
  for (size_t ja = 0; ja < mcsII2DI2->C2Function().size(); ja++)
    fTrack->mcs_c2function_II2DI2.push_back(mcsII2DI2->C2Function().at(ja));
  for (size_t ja = 0; ja < mcsII2DI2->TailsSize().size(); ja++)
    fTrack->mcs_tailssize_II2DI2.push_back(mcsII2DI2->TailsSize().at(ja));

  fTrack->mcs_pbest_II2DC  = mcsII2DC->bestMomentum();
  fTrack->mcs_perr_II2DC   = mcsII2DC->errorMomentum();
  fTrack->mcs_pmin_II2DC   = mcsII2DC->minMomentum();
  fTrack->mcs_pmax_II2DC   = mcsII2DC->maxMomentum();
  fTrack->mcs_alpha_II2DC  = mcsII2DC->alphaFit();
  fTrack->mcs_dalpha_II2DC = mcsII2DC->alphaFitErr();
  fTrack->mcs_beta_II2DC   = mcsII2DC->betaFit();
  fTrack->mcs_dbeta_II2DC  = mcsII2DC->betaFitErr();
  for (size_t ja = 0; ja < mcsII2DC->testMomentum().size(); ja++)
    fTrack->mcs_ptest_II2DC.push_back(mcsII2DC->testMomentum().at(ja));
  for (size_t ja = 0; ja < mcsII2DC->C2Function().size(); ja++)
    fTrack->mcs_c2function_II2DC.push_back(mcsII2DC->C2Function().at(ja));
  for (size_t ja = 0; ja < mcsII2DC->TailsSize().size(); ja++)
    fTrack->mcs_tailssize_II2DC.push_back(mcsII2DC->TailsSize().at(ja));
    
  if (hits.size() > 0) {
    fTrack->cryostat = hits[0]->WireID().Cryostat;
  }

  // Fill each hit
  for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
    sbn::TrackHitInfo hinfo = MakeHit(*hits[i_hit], hits[i_hit].key(), *thms[i_hit], track, sps[i_hit], calo, geo, clock_data, bt_serv);
    if (hinfo.h.plane == 0) {
      fTrack->hits0.push_back(hinfo);
    }
    else if (hinfo.h.plane == 1) {
      fTrack->hits1.push_back(hinfo);
    } 
    else if (hinfo.h.plane == 2) {
      fTrack->hits2.push_back(hinfo);
    }
  }
  std::cout << " skimmer collection count " << fTrack->hits2.size() << std::endl;

  // Hit summary info
  fTrack->hit_min_time_p0_tpcE = HitMinTime(fTrack->hits0, true, det);
  fTrack->hit_max_time_p0_tpcE = HitMaxTime(fTrack->hits0, true, det);
  fTrack->hit_min_time_p0_tpcW = HitMinTime(fTrack->hits0, false, det);
  fTrack->hit_max_time_p0_tpcW = HitMaxTime(fTrack->hits0, false, det);
  fTrack->hit_min_time_p1_tpcE = HitMinTime(fTrack->hits1, true, det);
  fTrack->hit_max_time_p1_tpcE = HitMaxTime(fTrack->hits1, true, det);
  fTrack->hit_min_time_p1_tpcW = HitMinTime(fTrack->hits1, false, det);
  fTrack->hit_max_time_p1_tpcW = HitMaxTime(fTrack->hits1, false, det);
  fTrack->hit_min_time_p2_tpcE = HitMinTime(fTrack->hits2, true, det);
  fTrack->hit_max_time_p2_tpcE = HitMaxTime(fTrack->hits2, true, det);
  fTrack->hit_min_time_p2_tpcW = HitMinTime(fTrack->hits2, false, det);
  fTrack->hit_max_time_p2_tpcW = HitMaxTime(fTrack->hits2, false, det);

  // Save information on a fit to the end of the track
  if (fDoTailFit) DoTailFit();

  // Save the Wire ADC values we need to
  for (auto const &w_pair: fWiresToSave) {
    geo::WireID wire = w_pair.first;

    if (rawdigits.count(wire)) {
      const raw::RawDigit &thisdigit = *rawdigits.at(wire);
      int min_tick = std::max(0, w_pair.second.first);
      int max_tick = std::min((int)thisdigit.NADC(), w_pair.second.second);

      // collect the adcs
      std::vector<short> adcs;
      for (int t = min_tick; t < max_tick; t++) {
        adcs.push_back(thisdigit.ADC(t));
      }

      WireInfo winfo;
      winfo.wire = wire.Wire;
      winfo.plane = wire.Plane;
      winfo.tpc = wire.TPC;
      winfo.channel = geo->PlaneWireToChannel(wire);
      winfo.tdc0 = min_tick;
      winfo.adcs = adcs;

      if (winfo.plane == 0) {
        fTrack->wires0.push_back(winfo);
      }
      else if (winfo.plane == 1) {
        fTrack->wires1.push_back(winfo);
      }
      else if (winfo.plane == 2) {
        fTrack->wires2.push_back(winfo);
      }
    }
  }
  // Sort the ADC values by wire
  std::sort(fTrack->wires0.begin(), fTrack->wires0.end(), [](auto const &lhs, auto const &rhs) {return lhs.wire < rhs.wire;});
  std::sort(fTrack->wires1.begin(), fTrack->wires1.end(), [](auto const &lhs, auto const &rhs) {return lhs.wire < rhs.wire;});
  std::sort(fTrack->wires2.begin(), fTrack->wires2.end(), [](auto const &lhs, auto const &rhs) {return lhs.wire < rhs.wire;});

  // get information on nearby tracks
  for (const GlobalTrackInfo &othr: tracks) {
    if (othr.ID == track.ID()) continue;

    if ((track.End() - othr.start).r() < 50. || (track.End() - othr.end).r() < 50.) {
      fTrack->tracks_near_end_dist.push_back(std::min((track.End() - othr.start).r(), (track.End() - othr.end).r())); 
      fTrack->tracks_near_end_costh.push_back(
        (track.End() - othr.start).r() < (track.End() - othr.end).r() ? 
          track.EndDirection().Dot(othr.dir) : track.EndDirection().Dot(othr.enddir)); 
    }
  }

  for (const GlobalTrackInfo &othr: tracks) {
    if (othr.ID == track.ID()) continue;

    if ((track.Start() - othr.start).r() < 50. || (track.Start() - othr.end).r() < 50.) {
      fTrack->tracks_near_start_dist.push_back(std::min((track.Start() - othr.start).r(), (track.Start() - othr.end).r())); 
      fTrack->tracks_near_start_costh.push_back(
        (track.Start() - othr.start).r() < (track.Start() - othr.end).r() ? 
          track.StartDirection().Dot(othr.dir) : track.StartDirection().Dot(othr.enddir)); 
    }
  }
  std::cout << " end filltrack " << std::endl;
}

void sbn::TrackCaloSkimmer::DoTailFit() {
  // Try fitting the constant and exponentials to the tail of dQ/dx v. RR on the collection plane
  std::vector<double> fit_rr;
  std::vector<double> fit_dqdx;

  for (const TrackHitInfo &h: fTrack->hits2) {
    if (h.oncalo && h.rr > 0. && h.rr < fTailFitResidualRange) {
      fit_rr.push_back(h.rr);
      fit_dqdx.push_back(h.dqdx);
    }
  }

  // TODO: should we throw an exception here??
  // 
  // If there is too much data to fit in the array, throw exception
  if (fit_rr.size() > MAX_N_FIT_DATA) {
    throw cet::exception("sbn::TrackCaloSkimmer::DoTailFit: More fitting points required (" 
      + std::to_string(fit_rr.size()) + ") than available in fit array (" + std::to_string(MAX_N_FIT_DATA) + ").\n");
  }

  // Copy the fit data to the global array
  for (unsigned i = 0; i < fit_rr.size() && i < MAX_N_FIT_DATA; i++) {
    FIT_RR[i] = fit_rr[i];
    FIT_DQDX[i] = fit_dqdx[i];
  }
  N_FIT_DATA = std::min(fit_rr.size(), MAX_N_FIT_DATA);

  fTrack->n_fit_point = N_FIT_DATA;
  if (fTrack->n_fit_point > 2) { // need more points than params
    // Fit the Exponential
    fFitExp.SetFCN(ExpResiduals);
    fFitExp.SetParameter(0, "A", *std::max_element(fit_dqdx.begin(), fit_dqdx.end()), 200, 0, 5000);
    fFitExp.SetParameter(1, "R", 10., 0.5, 0, 1000);
    fFitExp.ExecuteCommand("MIGRAD", 0, 0);

    double A = fFitExp.GetParameter(0);
    double R = fFitExp.GetParameter(1);

    int nparam;
    double param[2] {A, R};
    double residuals = -1;
    ExpResiduals(nparam, NULL, residuals, param, 0);

    //fTrack->exp_fit_A = A; 
    //fTrack->exp_fit_R = R;
    //fTrack->exp_fit_residuals = residuals;

    // Fit the Constant
    fFitConst.SetFCN(ConstResiduals);
    fFitConst.SetParameter(0, "C", std::accumulate(fit_dqdx.begin(), fit_dqdx.end(), 0.), 200, 0, 5000);
    fFitConst.ExecuteCommand("MIGRAD", 0, 0);

    //double C = fFitConst.GetParameter(0);
    //double cresiduals = -1;
    //ConstResiduals(nparam, NULL, cresiduals, &C, 0);
    //fTrack->const_fit_C = C;    
    //fTrack->const_fit_residuals = cresiduals;


  }
}

sbn::TrackHitInfo sbn::TrackCaloSkimmer::MakeHit(const recob::Hit &hit,
  unsigned hkey,
  const recob::TrackHitMeta &thm,
  const recob::Track &trk,
  const art::Ptr<recob::SpacePoint> &sp,
  const std::vector<art::Ptr<anab::Calorimetry>> &calo,
  const geo::GeometryCore *geo,
  const detinfo::DetectorClocksData &dclock,
  const cheat::BackTrackerService *bt_serv) {

  // TrackHitInfo to save
  sbn::TrackHitInfo hinfo;

  // information from the hit object
  if (trk.HasValidPoint(thm.Index())) {
    hinfo.h.integral = hit.Integral();
    hinfo.h.sumadc = hit.SummedADC();
    hinfo.h.width = hit.RMS();
    hinfo.h.time = hit.PeakTime();
    hinfo.h.mult = hit.Multiplicity();
    hinfo.h.wire = hit.WireID().Wire;
    hinfo.h.plane = hit.WireID().Plane;
    if (hinfo.h.plane == 2) std::cout << " adding skimmer index " << thm.Index() << std::endl;
    hinfo.h.channel = geo->PlaneWireToChannel(hit.WireID());
    hinfo.h.tpc = hit.WireID().TPC;
    hinfo.h.end = hit.EndTick();
    hinfo.h.start = hit.StartTick();
    hinfo.h.id = (int)hkey;
  }

  // Do back-tracking on each hit
  if (bt_serv) {
    // The default BackTracking function goes from (peak - width, peak + width).
    //
    // This time range does not match well hits with a non-Gaussian shape where
    // the Gaussian-fit-width does not replicate the width of the pulse. 
    //
    // Instead, we use the Hit (start, end) time range. This is also more relevant
    // for (e.g.) the SummedADC charge extraction method.
    //
    // Don't use this:
    // std::vector<sim::TrackIDE> ides = bt_serv->HitToTrackIDEs(dclock, hit);
    //
    // Use this:
    std::vector<sim::TrackIDE> ides = bt_serv->ChannelToTrackIDEs(dclock, hit.Channel(), hit.StartTick(), hit.EndTick());

    hinfo.h.truth.e = 0.;
    hinfo.h.truth.nelec = 0.;

    for (const sim::TrackIDE &ide: ides) {
      hinfo.h.truth.e += ide.energy;
      hinfo.h.truth.nelec += ide.numElectrons;
    }
  }
  else {
    hinfo.h.truth.e = -1.;
    hinfo.h.truth.nelec = -1.;
  }

  // look up the snippet
  sbn::TrackCaloSkimmer::Snippet snippet {hit.WireID(), hit.StartTick(), hit.EndTick()};
  if (!fSnippetCount.count(snippet)) {
    fSnippetCount[snippet] = 0;
    hinfo.i_snippet = 0;
  }
  else {
    fSnippetCount[snippet] ++;
    hinfo.i_snippet = fSnippetCount[snippet];
  }

  // Which wires to save
  int min_tick = (int)std::floor(hit.PeakTime() - fHitRawDigitsTickCollectWidth);
  int max_tick = (int)std::ceil(hit.PeakTime() + fHitRawDigitsTickCollectWidth);
  for (int wire = hinfo.h.wire - fHitRawDigitsWireCollectWidth; wire <= hinfo.h.wire + fHitRawDigitsWireCollectWidth; wire++) {
    geo::WireID w(hit.WireID(), wire);

    if (fWiresToSave.count(w)) {
      fWiresToSave.at(w).first = std::min(fWiresToSave.at(w).first, min_tick);
      fWiresToSave.at(w).second = std::max(fWiresToSave.at(w).second, max_tick);
    }
    else {
      fWiresToSave[w] = {min_tick, max_tick};
    }
  }

  // Information from the TrackHitMeta
  bool badhit = (thm.Index() == std::numeric_limits<unsigned int>::max()) ||
                    (!trk.HasValidPoint(thm.Index()));

  hinfo.ontraj = !badhit;

  // Save trajectory information if we can
  if (!badhit) {
    geo::Point_t loc = trk.LocationAtPoint(thm.Index());
    hinfo.tp.x = loc.X();
    hinfo.tp.y = loc.Y();
    hinfo.tp.z = loc.Z();

    geo::Vector_t dir = trk.DirectionAtPoint(thm.Index());
    hinfo.dir.x = dir.X();
    hinfo.dir.y = dir.Y();
    hinfo.dir.z = dir.Z();

    // And determine if the Hit is on a Calorimetry object
    for (const art::Ptr<anab::Calorimetry> &c: calo) {
      if (c->PlaneID().Plane != hinfo.h.plane) continue;

      // Found the plane! Now find the hit:
      for (unsigned i_calo = 0; i_calo < c->dQdx().size(); i_calo++) {
        if (c->TpIndices()[i_calo] == hkey) { // "TpIndices" match to the hit key
          // Fill the calo information associated with the hit
          hinfo.oncalo = true;
          hinfo.pitch = c->TrkPitchVec()[i_calo];
          hinfo.dqdx = c->dQdx()[i_calo];
          hinfo.rr = c->ResidualRange()[i_calo];
          break;
        }
      }
      break;
    }
  }

  // Save SpacePoint information
  if (sp && trk.HasValidPoint(thm.Index())) {
    hinfo.h.sp.x = sp->position().x();
    hinfo.h.sp.y = sp->position().y();
    hinfo.h.sp.z = sp->position().z();
    if(isnan(sp->position().x())) std::cout << " x nan " << thm.Index() << std::endl;
    if(isnan(sp->position().y())) std::cout << " y nan " << thm.Index() << std::endl;
    if(isnan(sp->position().z())) std::cout << " z nan " << thm.Index() << std::endl;

    hinfo.h.hasSP = true;
  }
  else {
    hinfo.h.hasSP = false;
  }

  return hinfo;
}

void sbn::TrackCaloSkimmer::FillDroppedTrackTruth(
  const detinfo::DetectorClocksData &clock_data,
  const std::vector<art::Ptr<recob::Hit>> &trkHits,
  const std::vector<art::Ptr<simb::MCParticle>> &mcparticles,
  const std::vector<art::Ptr<simb::MCParticle>> &droppedmcparticles,
  const std::vector<geo::BoxBoundedGeo> &active_volumes,
  const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
  const std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map,
  const std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map,
  const detinfo::DetectorPropertiesData &dprop,
  const geo::GeometryCore *geo) {

  // Lookup the true-particle match -- use utils in CAF
  std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, trkHits, true);
  float total_energy = CAFRecoUtils::TotalHitEnergy(clock_data, trkHits); 
  fTrack->truth.depE = total_energy / 1000. /* MeV -> GeV */;
  
  // sort highest energy match to lowest
  std::sort(matches.begin(), matches.end(),
      [](const auto &a, const auto &b) {
        return a.second > b.second;
      }
  );
std::cout << " dropped size " << droppedmcparticles.size() << std::endl;

  // loop on matches
  for (unsigned int j=0;j<matches.size();j++) {
    std::pair<int, float> match = matches[j];


    fTrack->truth.pur = match.second / total_energy;

    for (const art::Ptr<simb::MCParticle> &p_mcp: droppedmcparticles) {
     //std::cout << " nondropped trackid " << p_mcp->TrackId() << " pdg: " << p_mcp->PdgCode() << " process: " << p_mcp->Process() << std::endl;
     //std::cout << " nondropped trackid " << p_mcp->TrackId() << " energy: " << p_mcp->E()  << std::endl;

      if (p_mcp->TrackId() == match.first) {
        std::cout << "Matched reco track " << j << " Track ID: " << p_mcp->TrackId() << " pdg: " << p_mcp->PdgCode() << " process: " << p_mcp->EndProcess() << std::endl;
        fTrack->truth.p = TrueParticleInfo(*p_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
        fTrack->truth.eff = fTrack->truth.depE / (fTrack->truth.p.plane0VisE + fTrack->truth.p.plane1VisE + fTrack->truth.p.plane2VisE);

        // Lookup any delta
        for (const art::Ptr<simb::MCParticle> &d_mcp: mcparticles) {
          if (d_mcp->Mother() == p_mcp->TrackId() && // correct parent
              (d_mcp->Process() == "Decay" || d_mcp->Process() == "muIoni") && // correct process
              abs(d_mcp->PdgCode()) == 11) { // correct PDG code

            fTrack->truth.michel = TrueParticleInfo(*d_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo);
            break;
          }
        }
      }
    }

    
  }
}

DEFINE_ART_MODULE(sbn::TrackCaloSkimmer)
