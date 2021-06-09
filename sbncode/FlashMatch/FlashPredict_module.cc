///////////////////////////////////////////////////////////////////////
// Class:       FlashPredict
// Plugin Type: producer (art v3_04_00)
// File:        FlashPredict_module.cc
//
// Created: February-2020  Iker Loïc de Icaza Astiz (icaza@fnal.gov)
//
////////////////////////////////////////////////////////////////////////
#include "sbncode/FlashMatch/FlashPredict.hh"


FlashPredict::FlashPredict(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fPandoraProducer(p.get<std::string>("PandoraProducer"))
  , fSpacePointProducer(p.get<std::string>("SpacePointProducer"))
  , fOpHitProducer(p.get<std::string>("OpHitProducer"))
  , fOpHitARAProducer(p.get<std::string>("OpHitARAProducer", ""))
  // , fCaloProducer(p.get<std::string>("CaloProducer"))
  // , fTrackProducer(p.get<std::string>("TrackProducer"))
  , fClockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob())
  , fTickPeriod(fClockData.OpticalClock().TickPeriod()) // us
  , fBeamWindowStart(p.get<double>("BeamWindowStart")) //us // TODO: should come from service
  , fBeamWindowEnd(p.get<double>("BeamWindowEnd"))// us // TODO: should come from service
  , fFlashStart(p.get<double>("FlashStart")) // in us w.r.t. flash time
  , fFlashEnd(p.get<double>("FlashEnd"))  // in us w.r.t flash time
  , fTimeBins(unsigned(1/fTickPeriod * (fBeamWindowEnd - fBeamWindowStart)))
  , fSelectNeutrino(p.get<bool>("SelectNeutrino", true))
  , fOnlyCollectionWires(p.get<bool>("OnlyCollectionWires", true))
  , fUseUncoatedPMT(p.get<bool>("UseUncoatedPMT", false))
  , fUseOppVolMetric(p.get<bool>("UseOppVolMetric", false))
  , fUseARAPUCAS(p.get<bool>("UseARAPUCAS", false))
  // , fUseCalo(p.get<bool>("UseCalo", false))
  , fInputFilename(p.get<std::string>("InputFileName")) // root file with score metrics
  , fNoAvailableMetrics(p.get<bool>("NoAvailableMetrics", false))
  , fMakeTree(p.get<bool>("MakeTree", false))
  , fChargeToNPhotonsShower(p.get<double>("ChargeToNPhotonsShower", 1.0))  // ~40000/1600
  , fChargeToNPhotonsTrack(p.get<double>("ChargeToNPhotonsTrack", 1.0))  // ~40000/1600
  , fMinHitQ(p.get<double>("MinHitQ", 0.0))
  , fMinSpacePointQ(p.get<double>("MinSpacePointQ", 0.0))
  , fMinParticleQ(p.get<double>("MinParticleQ", 0.0))
  , fMinSliceQ(p.get<double>("MinSliceQ", 0.0))
  , fMaxFlashes(p.get<unsigned>("MaxFlashes", 1))
  , fMinOpHPE(p.get<double>("MinOpHPE", 0.0))
  , fMinFlashPE(p.get<double>("MinFlashPE", 0.0))
  , fDetector(detectorName(fGeometry->DetectorName()))
  , fSBND((fDetector == "SBND") ? true : false )
  , fICARUS((fDetector == "ICARUS") ? true : false )
  , fPDMapAlgPtr(art::make_tool<opdet::PDMapAlg>(p.get<fhicl::ParameterSet>("PDMapAlg")))
  , fCryostat(p.get<int>("Cryostat", 0)) //set =0 or =1 for ICARUS to match reco chain selection
  , fGeoCryo(std::make_unique<geo::CryostatGeo>(fGeometry->Cryostat(fCryostat)))
  , fNBins(p.get<int>("n_bins"))
  , fDriftDistance(p.get<double>("DriftDistance"))// rounded up for binning
  , fNTPC(fGeometry->NTPC())
  , fOpDetNormalizer((fSBND) ? 4 : 1)
  , fTermThreshold(p.get<double>("ThresholdTerm", 30.))
{
  produces< std::vector<sbn::SimpleFlashMatch> >();
  produces< art::Assns <recob::PFParticle, sbn::SimpleFlashMatch> >();

  for (size_t t = 0; t < fNTPC; t++) {
    const geo::TPCGeo& tpcg = fGeoCryo->TPC(t);
    fWiresX_gl.push_back(tpcg.LastPlane().GetCenter().X());
  }
  fWiresX_gl.unique([](double l, double r) { return std::abs(l - r) < 0.00001;});

  if(fSBND && !fICARUS) {
    if(fUseOppVolMetric) {
      throw cet::exception("FlashPredict")
        << "UseOppVolMetric: " << std::boolalpha << fUseOppVolMetric << "\n"
        << "Not supported on SBND. Stopping.";
    }
    fTPCPerDriftVolume = 1;
    fDriftVolumes = fNTPC/fTPCPerDriftVolume;
  }
  else if(fICARUS && !fSBND) {
    if(fUseUncoatedPMT || fUseARAPUCAS) {
      throw cet::exception("FlashPredict")
        << "UseUncoatedPMT: " << std::boolalpha << fUseUncoatedPMT << ",\n"
        << "UseARAPUCAS:    " << std::boolalpha << fUseARAPUCAS << "\n"
        << "Not supported on ICARUS. Stopping.\n";
    }
    fTPCPerDriftVolume = 2;
    fDriftVolumes = fNTPC/fTPCPerDriftVolume;
  }
  else {
      throw cet::exception("FlashPredict")
        << "Detector: " << fDetector
        << ", not supported. Stopping.\n";
  }

  // TODO no point on having fCryostat as parameter, user whatever comes from geometry
  if (fSBND && fCryostat != 0) {
    throw cet::exception("FlashPredict")
      << "SBND has only one cryostat. \n"
      << "Check Detector and Cryostat parameter." << std::endl;
  }
  else if (fICARUS && fCryostat > 1) {
    throw cet::exception("FlashPredict")
      << "ICARUS has only two cryostats. \n"
      << "Check Detector and Cryostat parameter." << std::endl;
  }

  if (fMakeTree) initTree();

  loadMetrics();

  consumes<std::vector<recob::PFParticle>>(fPandoraProducer);
  // consumes<std::vector<recob::Slice>>(fPandoraProducer);
  consumes<art::Assns<recob::SpacePoint, recob::PFParticle>>(fPandoraProducer);
  consumes<std::vector<recob::SpacePoint>>(fSpacePointProducer);
  consumes<art::Assns<recob::Hit, recob::SpacePoint>>(fSpacePointProducer);
  consumes<std::vector<recob::OpHit>>(fOpHitProducer);
  if(fUseARAPUCAS && !fOpHitARAProducer.empty())
    consumes<std::vector<recob::OpHit>>(fOpHitARAProducer);
} // FlashPredict::FlashPredict(fhicl::ParameterSet const& p)


void FlashPredict::produce(art::Event& evt)
{
  // sFM is an alias for sbn::SimpleFlashMatch
  std::unique_ptr< std::vector<sFM> >
    sFM_v(new std::vector<sFM>);
  std::unique_ptr< art::Assns <recob::PFParticle, sFM> >
    pfp_sFM_assn_v(new art::Assns<recob::PFParticle, sFM>);

  // reset TTree variables
  _evt = evt.event();
  _sub = evt.subRun();
  _run = evt.run();
  // _slices        = 0;
  _flash_time    = -9999.;
  _flash_pe      = -9999.;
  _flash_unpe    = -9999.;
  _flash_r       = -9999.;
  _flash_ratio   = -9999.;
  _hypo_x        = -9999.;
  _hypo_x_err    = -9999.;
  _hypo_x_rr     = -9999.;
  _hypo_x_ratio  = -9999.;
  // _hypo_x_fit    = -9999.;
  _score         = -9999.;
  bk.events++;

  // grab PFParticles in event
  const auto pfps_h =
    evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraProducer);
  if (fSelectNeutrino &&
      !pfpNeutrinoOnEvent(pfps_h)) {
    mf::LogInfo("FlashPredict")
      << "No pfp neutrino on event. Skipping...";
    bk.nopfpneutrino++;
    updateBookKeeping();
    for(size_t pId=0; pId<pfps_h->size(); pId++) {
      if(!pfps_h->at(pId).IsPrimary()) continue;
      const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
      sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                           Flash(kNoScrPE), Score(kNoPFPInEvt)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
    }
    evt.put(std::move(sFM_v));
    evt.put(std::move(pfp_sFM_assn_v));
    return;
  }

  // auto const& slice_h = evt.getValidHandle<std::vector<recob::Slice>>(fPandoraProducer);
  // _slices = slice_h.product()->size();
  // if (_slices == 0) {
  //   mf::LogWarning("FlashPredict")
  //     << "No recob:Slice on event. Skipping...";
  //   bk.noslice++;
  //   updateBookKeeping();
  //   for(size_t pId=0; pId<pfps_h->size(); pId++) {
  //     if(!pfps_h->at(pId).IsPrimary()) continue;
  //     const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
  //     sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
  //                          Flash(kNoScrPE), Score(kNoSlcInEvt)));
  //     util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
  //   }
  //   evt.put(std::move(sFM_v));
  //   evt.put(std::move(pfp_sFM_assn_v));
  //   return;
  // }

  // load OpHits previously created
  art::Handle<std::vector<recob::OpHit>> ophits_h;
  evt.getByLabel(fOpHitProducer, ophits_h);
  if(!ophits_h.isValid()) {
    mf::LogError("FlashPredict")
      << "No optical hits from producer module "
      << fOpHitProducer;
    bk.nonvalidophit++;
    updateBookKeeping();
    for(size_t pId=0; pId<pfps_h->size(); pId++) {
      if(!pfps_h->at(pId).IsPrimary()) continue;
      const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
      sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                           Flash(kNoScrPE), Score(kNoOpHInEvt)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
    }
    evt.put(std::move(sFM_v));
    evt.put(std::move(pfp_sFM_assn_v));
    return;
  }

  std::vector<recob::OpHit> opHits(ophits_h->size());
  copyOpHitsInBeamWindow(opHits, ophits_h);

  if(fUseARAPUCAS && !fOpHitARAProducer.empty()){
    art::Handle<std::vector<recob::OpHit>> ophits_ara_h;
    evt.getByLabel(fOpHitARAProducer, ophits_ara_h);
    if(!ophits_ara_h.isValid()) {
      mf::LogError("FlashPredict")
        << "Non valid ophits from ARAPUCAS"
        << "\nfUseARAPUCAS: " << std::boolalpha << fUseARAPUCAS
        << "\nfOpHitARAProducer: " << fOpHitARAProducer;
    }
    else{
      std::vector<recob::OpHit> opHitsARA(ophits_ara_h->size());
      copyOpHitsInBeamWindow(opHitsARA, ophits_ara_h);
      opHits.insert(opHits.end(),
                    opHitsARA.begin(), opHitsARA.end());
    }
  }


  std::vector<recob::OpHit> opHitsLeft, opHitsRght;
  const std::vector<SimpleFlash> simpleFlashes = (fSBND) ?
    makeSimpleFlashes(opHits, opHitsLeft, opHitsRght) : makeSimpleFlashes(opHits);
  if(simpleFlashes.empty()){
    mf::LogWarning("FlashPredict")
      << "\nNo OpHits in beam window,"
      << "\nor the integral is less than " << fMinFlashPE
      << "\nSkipping...";
    bk.nullophittime++;
    updateBookKeeping();
    for(size_t pId=0; pId<pfps_h->size(); pId++) {
      if(!pfps_h->at(pId).IsPrimary()) continue;
      const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
      sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                           Flash(kNoScrPE), Score(kNoOpHInEvt)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
    }
    evt.put(std::move(sFM_v));
    evt.put(std::move(pfp_sFM_assn_v));
    return;
  }

  // std::set<unsigned> tpcWithOpH;
  // if(fSBND) {// no point for ICARUS
  //   for(auto oph=fOpH_beg; oph!=fOpH_end; ++oph){
  //     tpcWithOpH.insert(sbndPDinTPC(oph->OpChannel()));
  //     if(tpcWithOpH.size() == fNTPC) break;
  //   }
  // }


  ChargeDigestMap chargeDigestMap = makeChargeDigest(evt, pfps_h);

  std::map<unsigned, FlashMetrics> flashMetricsMap;
  for(auto& chargeDigest : chargeDigestMap) {
    //const size_t pId = chargeDigest.second.pId;
    const int pfpPDGC = chargeDigest.second.pfpPDGC;
    const auto& pfp_ptr = chargeDigest.second.pfp_ptr;
    const auto& qClusters = chargeDigest.second.qClusters;
    const auto& tpcWithHits = chargeDigest.second.tpcWithHits;

    // if(fSBND){// because SBND has an opaque cathode
    //   std::set<unsigned> tpcWithHitsOpH;
    //   std::set_intersection(tpcWithHits.begin(), tpcWithHits.end(),
    //                         tpcWithOpH.begin(), tpcWithOpH.end(),
    //                         std::inserter(tpcWithHitsOpH, tpcWithHitsOpH.begin()));
    //   if (tpcWithHitsOpH.size() == 0) {
    //     mf::LogWarning("FlashPredict") << "No OpHits where there's charge. Skipping...";
    //     bk.no_oph_hits++;
    //     mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
    //     sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
    //                          Flash(kNoScrPE), Score(kQNoOpHScr)));
    //     util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
    //     continue;
    //   }
    // }

    if(!computeChargeMetrics(qClusters)){
      mf::LogWarning("FlashPredict") << "Clusters with No Charge. Skipping...";
      bk.no_charge++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                           Flash(kNoScrPE), Score(kNoChrgScr)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      continue;
    }


    // TODO: finish rest of items for multiflash support:
    for(auto& simpleFlash : simpleFlashes) {
      if(simpleFlash.flashUId != kBothOpHs) continue;// TODO
      // check if metrics already computed, if so skip computing and grab them
      if(!computeFlashMetrics(simpleFlash)){
        printMetrics("ERROR", pfpPDGC, tpcWithHits, 0, mf::LogError("FlashPredict"));
        bk.no_flash_pe++;
        mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
        sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                             Flash(kNoScrPE), Score(k0VUVPEScr)));
        util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
        continue;
      }

      // TODO:
      // update metrics map outside computeflashmetrics

      // TODO:
      // computeScore is also inside the loop, the assigned score and
      // variables that are writen to the root tree and the artobject
      // should be those of the lowest score
      if(computeScore(tpcWithHits, pfpPDGC)){
        if (fMakeTree) {_flashmatch_nuslice_tree->Fill();}
        bk.scored_pfp++;
        mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
        Charge charge{_charge_q, TVector3(_charge_x_gl, _charge_y, _charge_z)};
        Flash flash{_flash_pe, TVector3(_flash_x_gl, _flash_y, _flash_z)};
        Score score{_score, _scr_y, _scr_z, _scr_rr, _scr_ratio};
        sFM_v->push_back(sFM(true, _flash_time, charge, flash, score));
        util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      }
      else{
        mf::LogError("FlashPredict") << "ERROR: score <= 0. Dumping info."
                                     << "\n_score:     " << _score
                                     << "\n_scr_y:     " << _scr_y
                                     << "\n_scr_z:     " << _scr_z
                                     << "\n_scr_rr:    " << _scr_rr
                                     << "\n_scr_ratio: " << _scr_ratio;
        printMetrics("ERROR", pfpPDGC, tpcWithHits, 0, mf::LogError("FlashPredict"));
      }
    }
  } // chargeDigestMap: PFparticles that pass criteria
  bk.events_processed++;
  updateBookKeeping();

  evt.put(std::move(sFM_v));
  evt.put(std::move(pfp_sFM_assn_v));

}// end of producer module


void FlashPredict::initTree(void)
{
  art::ServiceHandle<art::TFileService> tfs;
  _flashmatch_nuslice_tree = tfs->make<TTree>("nuslicetree", "nu FlashPredict tree");
  _flashmatch_nuslice_tree->Branch("evt", &_evt, "evt/I");
  _flashmatch_nuslice_tree->Branch("run", &_run, "run/I");
  _flashmatch_nuslice_tree->Branch("sub", &_sub, "sub/I");
  // _flashmatch_nuslice_tree->Branch("slices", &_slices, "slices/I");
  _flashmatch_nuslice_tree->Branch("flash_time", &_flash_time, "flash_time/D");
  _flashmatch_nuslice_tree->Branch("flash_x_gl", &_flash_x_gl, "flash_x_gl/D");
  _flashmatch_nuslice_tree->Branch("flash_x", &_flash_x, "flash_x/D");
  _flashmatch_nuslice_tree->Branch("flash_y", &_flash_y, "flash_y/D");
  _flashmatch_nuslice_tree->Branch("flash_z", &_flash_z, "flash_z/D");
  _flashmatch_nuslice_tree->Branch("flash_r", &_flash_r, "flash_r/D");
  _flashmatch_nuslice_tree->Branch("flash_pe", &_flash_pe, "flash_pe/D");
  _flashmatch_nuslice_tree->Branch("flash_unpe", &_flash_unpe, "flash_unpe/D");
  _flashmatch_nuslice_tree->Branch("flash_ratio", &_flash_ratio, "flash_ratio/D");
  // TODO: add charge_time?
  _flashmatch_nuslice_tree->Branch("charge_x_gl", &_charge_x_gl, "charge_x_gl/D");
  _flashmatch_nuslice_tree->Branch("charge_x", &_charge_x, "charge_x/D");
  _flashmatch_nuslice_tree->Branch("charge_y", &_charge_y, "charge_y/D");
  _flashmatch_nuslice_tree->Branch("charge_z", &_charge_z, "charge_z/D");
  _flashmatch_nuslice_tree->Branch("charge_q", &_charge_q, "charge_q/D");
  _flashmatch_nuslice_tree->Branch("hypo_x", &_hypo_x, "hypo_x/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_err", &_hypo_x_err, "hypo_x_err/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_rr", &_hypo_x_rr, "hypo_x_rr/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_ratio", &_hypo_x_ratio, "hypo_x_ratio/D");
  // _flashmatch_nuslice_tree->Branch("hypo_x_fit", &_hypo_x_fit, "hypo_x_fit/D");
  _flashmatch_nuslice_tree->Branch("score", &_score, "score/D");
  _flashmatch_nuslice_tree->Branch("scr_y", &_scr_y, "scr_y/D");
  _flashmatch_nuslice_tree->Branch("scr_z", &_scr_z, "scr_z/D");
  _flashmatch_nuslice_tree->Branch("scr_rr", &_scr_rr, "scr_rr/D");
  _flashmatch_nuslice_tree->Branch("scr_ratio", &_scr_ratio, "scr_ratio/D");
}


void FlashPredict::loadMetrics()
{
  // TODO: fill histos with less repetition and range for loops
  // read histograms and fill vectors for match score calculation
  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  if(!sp.find_file(fInputFilename, fname)) {
    mf::LogError("FlashPredict")
      << "Could not find the light-charge match root file '"
      << fInputFilename << "' on FW_SEARCH_PATH\n";
    throw cet::exception("FlashPredict")
      << "Could not find the light-charge match root file '"
      << fInputFilename << "' on FW_SEARCH_PATH\n";
  }
  mf::LogInfo("FlashPredict") << "Opening file with metrics: " << fname;
  TFile *infile = new TFile(fname.c_str(), "READ");
  auto metricsInFile = infile->GetListOfKeys();
  if(!metricsInFile->Contains("dy_h1") ||
     !metricsInFile->Contains("dz_h1") ||
     !metricsInFile->Contains("rr_h1") ||
     !metricsInFile->Contains("pe_h1") ||
     !metricsInFile->Contains("rr_fit_l") ||
     !metricsInFile->Contains("rr_fit_m") ||
     !metricsInFile->Contains("rr_fit_h") ||
     !metricsInFile->Contains("pe_fit_l") ||
     !metricsInFile->Contains("pe_fit_m") ||
     !metricsInFile->Contains("pe_fit_h"))
  {
    mf::LogError("FlashPredict")
      << "The metrics file '" << fname << "' lacks at least one metric.";
    throw cet::exception("FlashPredict")
      << "The metrics file '" << fname << "'lacks at least one metric.";
  }
  //
  TH1 *temphisto = (TH1*)infile->Get("dy_h1");
  int bins = 0;
  bins = temphisto->GetNbinsX();
  if (bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for dy " << bins << " bins " << std::endl;
    bins = 1;
    dy_means.push_back(0);
    dy_spreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= bins; ++ib) {
      dy_means.push_back(temphisto->GetBinContent(ib));
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in dy" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      dy_spreads.push_back(tt);
    }
  }
  //
  temphisto = (TH1*)infile->Get("dz_h1");
  bins = temphisto->GetNbinsX();
  if (bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for dz " << bins << " bins " << std::endl;
    bins = 1;
    dz_means.push_back(0);
    dz_spreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= bins; ++ib) {
      dz_means.push_back(temphisto->GetBinContent(ib));
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in dz" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      dz_spreads.push_back(tt);
    }
  }
  //
  temphisto = (TH1*)infile->Get("rr_h1");
  bins = temphisto->GetNbinsX();
  if (bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for rr " << bins << " bins " << std::endl;
    bins = 1;
    rr_means.push_back(0);
    rr_spreads.push_back(0.001);
  }
  else {
    // std::vector<double> x, yL, yM, yH;
    for (int ib = 1; ib <= bins; ++ib) {
      double me = temphisto->GetBinContent(ib);
      rr_means.push_back(me);
      // x.push_back(temphisto->GetBinCenter(ib));
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in rr" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      rr_spreads.push_back(tt);
      // yL.push_back(me - tt);
      // yM.push_back(me);
      // yH.push_back(me + tt);
    }
    // rrMax = *std::max_element(yH.begin(), yH.end());
    // rr_l_InvSpl = TSpline3("rr_l_InvSpl", yL.data(), x.data(), yL.size());
    // rr_m_InvSpl = TSpline3("rr_m_InvSpl", yM.data(), x.data(), yM.size());
    // rr_h_InvSpl = TSpline3("rr_h_InvSpl", yH.data(), x.data(), yH.size());
    unsigned s = 0;
    for(auto& rrF : fRRFits){
      std::string nold = "rr_fit_" + kSuffixes[s];
      std::string nnew = "rrFit_" + kSuffixes[s];
      TF1* tempF1 = (TF1*)infile->Get(nold.c_str());
      auto params = tempF1->GetParameters();
      rrF.f = std::make_unique<TF1>(nnew.c_str(), "pol3", 0., fDriftDistance);
      rrF.f->SetParameters(params);
      rrF.min = rrF.f->GetMinimum(0., fDriftDistance, kEps);
      rrF.max = rrF.f->GetMaximum(0., fDriftDistance, kEps);
      s++;
    }
  }
  //
  temphisto = (TH1*)infile->Get("pe_h1");
  bins = temphisto->GetNbinsX();
  if (bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for pe " << bins << " bins " << std::endl;
    bins = 1;
    pe_means.push_back(0);
    pe_spreads.push_back(0.001);
  }
  else {
    // std::vector<double> x, yL, yM, yH;
    for (int ib = 1; ib <= bins; ++ib) {
      double me = temphisto->GetBinContent(ib);
      pe_means.push_back(me);
      // x.push_back(temphisto->GetBinCenter(ib));
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in pe" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      pe_spreads.push_back(tt);
      // yL.push_back(me - tt);
      // yM.push_back(me);
      // yH.push_back(me + tt);
    }
    // peMax = *std::max_element(yH.begin(), yH.end());
    // pe_l_InvSpl = TSpline3("pe_l_InvSpl", yL.data(), x.data(), yL.size());
    // pe_m_InvSpl = TSpline3("pe_m_InvSpl", yM.data(), x.data(), yM.size());
    // pe_h_InvSpl = TSpline3("pe_h_InvSpl", yH.data(), x.data(), yH.size());
    unsigned s = 0;
    for(auto& peF : fRatioFits){
      std::string nold = "pe_fit_" + kSuffixes[s];
      std::string nnew = "peFit_" + kSuffixes[s];
      TF1* tempF1 = (TF1*)infile->Get(nold.c_str());
      auto params = tempF1->GetParameters();
      peF.f = std::make_unique<TF1>(nnew.c_str(), "pol3", 0., fDriftDistance);
      peF.f->SetParameters(params);
      peF.min = peF.f->GetMinimum(0., fDriftDistance, kEps);
      peF.max = peF.f->GetMaximum(0., fDriftDistance, kEps);
      s++;
    }
  }

  infile->Close();
  delete infile;
  mf::LogInfo("FlashPredict") << "Finish loading metrics";
}


bool FlashPredict::computeChargeMetrics(const flashmatch::QCluster_t& qClusters)
{
  double xave = 0.; double yave = 0.;
  double zave = 0.; double norm = 0.;
  for (auto& qp : qClusters) {
    xave += qp.q * qp.x;
    yave += qp.q * qp.y;
    zave += qp.q * qp.z;
    norm += qp.q;
  }
  if (norm > 0.) {
    _charge_q = norm;
    _charge_x_gl = xave / norm;
    _charge_x = driftDistance(_charge_x_gl);
    _charge_y = yave / norm;
    _charge_z = zave / norm;
    return true;
  }
  return false;
}


bool FlashPredict::computeFlashMetrics(const SimpleFlash& simpleFlash)
{
  // unsigned flashMetricsId = fFlashCounter;
  // unsigned tens = 10;
  // for(unsigned t : tpcWithHits) {
  //   flashMetricsId += tens * t;
  //   tens *= 10;
  // }
  // if(auto it=flashMetricsMap.find(flashMetricsId); it!=flashMetricsMap.end()){
  //   mf::LogDebug("FlashPredict") << "Reusing metrics previously computed, "
  //                                << "for flashMetricsId " << flashMetricsId;
  //   auto flashMetrics = it->second;
  //   if(!flashMetrics.metric_ok) return false;
  //   updateFlashMetrics(flashMetrics);
  //   return true;
  // }

  const OpHitIt opH_beg = simpleFlash.opH_beg;
  const OpHitIt opH_end = simpleFlash.opH_end;

  // NOTE: opHMax_X (and later _flash_x) holds the X coordinate of the
  //       opdet that registered the most PEs in the flash.
  // TODO: change to peSumMax_wallX (and later _flash_x) to hold the X
  //       coordinate of the wall that registered the most PEs overall
  auto compareOpHits = [] (const recob::OpHit& oph1, const recob::OpHit& oph2)->bool
    { return oph1.PE() < oph2.PE(); };
  auto opHMax = std::max_element(opH_beg, opH_end, compareOpHits);
  double opHMax_X =
    fGeometry->OpDetGeoFromOpChannel(opHMax->OpChannel()).GetCenter().X();

  double sum = 0.;
  double sum_PE = 0.;
  double sum_PE2 = 0.;
  double sum_unPE = 0.;
  double sum_visARA_PE = 0.;
  double sum_PE2Y  = 0.; double sum_PE2Z  = 0.;
  double sum_PE2Y2 = 0.; double sum_PE2Z2 = 0.;

  for(auto oph=opH_beg; oph!=opH_end; ++oph){
    int opChannel = oph->OpChannel();
    // TODO remove this if !isSBNDPDRelevant(), eventually
    // if(fSBND && !isSBNDPDRelevant(opChannel, tpcWithHits)) continue;
    auto opDet = fGeometry->OpDetGeoFromOpChannel(opChannel);
    auto opDetXYZ = opDet.GetCenter();

    bool is_pmt_vis = false, is_ara_vis = false;
    if(fSBND){// because VIS light
      auto op_type = fPDMapAlgPtr->pdType(opChannel);
      if(op_type == "pmt_uncoated") {
        if(!fUseUncoatedPMT) continue;
        is_pmt_vis = true, is_ara_vis = false;
      }
      else if(op_type == "xarapuca_vis" || op_type == "arapuca_vis") {
        is_pmt_vis = false, is_ara_vis = true;
        // if !fUseArapucas, they weren't loaded at all
      }
    }

    double ophPE  = oph->PE();
    double ophPE2 = ophPE * ophPE;
    sum       += 1.0;
    sum_PE    += ophPE;
    sum_PE2   += ophPE2;
    sum_PE2Y  += ophPE2 * opDetXYZ.Y();
    sum_PE2Z  += ophPE2 * opDetXYZ.Z();
    sum_PE2Y2 += ophPE2 * opDetXYZ.Y() * opDetXYZ.Y();
    sum_PE2Z2 += ophPE2 * opDetXYZ.Z() * opDetXYZ.Z();

    if(fICARUS){
      // TODO: change this if block
      if(fUseOppVolMetric){
        unsigned pdVolume = icarusPDinTPC(opChannel)/fTPCPerDriftVolume;
        geo::Point_t q(_charge_x_gl, _charge_y, _charge_z);
        unsigned qVolume = driftVolume(_charge_x_gl);
        if(qVolume < fDriftVolumes){
          if (pdVolume != qVolume){
            sum_unPE += ophPE;
          }
        }
      }
      // // TODO: for this block instead, needs testing!
      // if(fUseOppVolMetric && std::abs(peSumMax_wallX-opDetXYZ.X()) < 10.){
      //   sum_unPE += ophPE;
      // }
    }
    else {// fSBND
      if(fUseUncoatedPMT && is_pmt_vis) {
        sum_unPE += ophPE;
      }
      else if(fUseARAPUCAS && is_ara_vis) {
        sum_visARA_PE += ophPE;
      }
    }
  } // for opHits

  if (sum_PE > 0.) {
    _flash_pe    = sum_PE;
    _flash_unpe  = sum_unPE;
    _flash_ratio = fOpDetNormalizer * _flash_unpe / _flash_pe;
    if(fUseARAPUCAS) {
      _flash_unpe  += sum_visARA_PE;
      _flash_ratio = (fOpDetNormalizer * sum_unPE  + sum_visARA_PE ) / _flash_pe;
    }
    _flash_y  = sum_PE2Y / sum_PE2;
    _flash_z  = sum_PE2Z / sum_PE2;
    _flash_r = std::sqrt(
      std::abs(sum_PE2Y2 + sum_PE2Z2 + sum_PE2 * (_flash_y * _flash_y + _flash_z * _flash_z)
       - 2.0 * (_flash_y * sum_PE2Y + _flash_z * sum_PE2Z) ) / sum_PE2);
    std::tie(_hypo_x, _hypo_x_err, _hypo_x_rr, _hypo_x_ratio) =
      hypoFlashX_fits(_flash_r, _flash_ratio);
    // TODO: using _hypo_x make further corrections to _flash_time to
    // account for light transport time and/or rising edge
    // _flash_time = timeCorrections(simpleFlash.maxpeak_time, _hypo_x);
    _flash_time = simpleFlash.maxpeak_time;
    _flash_x = opHMax_X; // TODO: _flash_x = peSumMax_wallX;
    _flash_x_gl = flashXGl(_hypo_x, _flash_x);
    // storeFlashMetrics(flashMetricsId, flashMetricsMap);
    return true;
  }
  else {
    // flashMetricsMap[flashMetricsId] = FlashMetrics();//faulty flash
    std::string channels;
    for(auto oph=opH_beg; oph!=opH_end; ++oph) channels += std::to_string(oph->OpChannel()) + ' ';
    // std::string tpcs;
    // for(auto itpc: tpcWithHits) tpcs += std::to_string(itpc) + ' ';
    mf::LogError("FlashPredict")
      << "Really odd that I landed here, this shouldn't had happen.\n"
      << "sum:          \t" << sum << "\n"
      << "sum_PE:       \t" << sum_PE << "\n"
      << "sum_unPE:     \t" << sum_unPE << "\n"
      // << "tpcWithHits:  \t" << tpcs << "\n"
      << "opHits size:  \t" << std::distance(opH_beg, opH_end) << "\n"
      << "channels:     \t" << channels << std::endl;
    _flash_y = 0.;
    _flash_z = 0.;
    _flash_r = 0.;
    _flash_pe = 0.;
    _flash_unpe = 0.;
    _flash_ratio = 0.;
    return false;
  }
}


bool FlashPredict::computeScore(const std::set<unsigned>& tpcWithHits,
                                const int pdgc)
{
  _score = 0.;
  unsigned tcount = 0;
  int isl = int(fNBins * (_charge_x / fDriftDistance));
  auto out = mf::LogWarning("FlashPredict");

  _scr_y = scoreTerm(_flash_y, _charge_y, dy_means[isl], dy_spreads[isl]);
  if (_scr_y > fTermThreshold) printMetrics("Y", pdgc, tpcWithHits, _scr_y, out);
  _score += _scr_y;
  tcount++;
  _scr_z = scoreTerm(_flash_z, _charge_z, dz_means[isl], dz_spreads[isl]);
  if (_scr_z > fTermThreshold) printMetrics("Z", pdgc, tpcWithHits, _scr_z, out);    
  _score += _scr_z;
  tcount++;
  _scr_rr = scoreTerm(_flash_r, rr_means[isl], rr_spreads[isl]);
  if (_scr_rr > fTermThreshold) printMetrics("R", pdgc, tpcWithHits, _scr_rr, out);
  _score += _scr_rr;
  tcount++;
  if (fUseUncoatedPMT || fUseOppVolMetric) {
    _scr_ratio = scoreTerm(_flash_ratio, pe_means[isl], pe_spreads[isl]);
    if (_scr_ratio > fTermThreshold) printMetrics("RATIO", pdgc, tpcWithHits, _scr_ratio, out);
    _score += _scr_ratio;
    tcount++;
  }
  mf::LogDebug("FlashPredict")
    << "score:\t" << _score << "using " << tcount << " terms";
  if(_score > 0.) return true;
  else return false;
}


// double FlashPredict::hypoFlashX_splines() const
// {
//   double lX = 0., mX = 0., hX = fDriftDistance;
//   if(_flash_r < rrMax){//don't bother to interpolate otherwise
//     double rr_lX = rr_l_InvSpl.Eval(_flash_r);
//     double rr_mX = rr_m_InvSpl.Eval(_flash_r);
//     double rr_hX = rr_h_InvSpl.Eval(_flash_r);
//     if(0.<rr_mX && rr_mX<fDriftDistance) mX = rr_mX;
//     else mX = (_flash_r<rr_means[0]) ? 0. : fDriftDistance;
//     hX = (0.<rr_hX && rr_hX<fDriftDistance) ? rr_hX : mX;
//     lX = (0.<rr_lX && rr_lX<fDriftDistance) ? rr_lX : mX;
//     if(lX == hX) hX = fDriftDistance, lX = fDriftDistance/3.;//TODO: hack
//   }
//   double rr_hypoX =  (lX + hX)/2.;
//   double rr_hypoXWgt =  1./std::abs(lX - hX);
//   if(!fUseUncoatedPMT && !fUseOppVolMetric)
//     return rr_hypoX;
//   lX = 0., mX = 0., hX = fDriftDistance;
//   if(_flash_ratio < peMax){//don't bother to interpolate otherwise
//     double pe_lX = pe_l_InvSpl.Eval(_flash_ratio);
//     double pe_mX = pe_m_InvSpl.Eval(_flash_ratio);
//     double pe_hX = pe_h_InvSpl.Eval(_flash_ratio);
//     if(0.<pe_mX && pe_mX<fDriftDistance) mX = pe_mX;
//     else mX = (_flash_ratio<pe_means[0]) ? 0. : fDriftDistance;
//     hX = (0.<pe_hX && pe_hX<fDriftDistance) ? pe_hX : mX;
//     lX = (0.<pe_lX && pe_lX<fDriftDistance) ? pe_lX : mX;
//     if(lX == hX) hX = fDriftDistance/2., lX = 0.;//TODO: hack
//   }
//   double pe_hypoX =  (lX + hX)/2.;
//   double pe_hypoXWgt =  1./std::abs(lX - hX);
//   return (rr_hypoX*rr_hypoXWgt + pe_hypoX*pe_hypoXWgt) / (rr_hypoXWgt + pe_hypoXWgt);
// }


std::tuple<double, double, double, double> FlashPredict::hypoFlashX_fits(
  double flash_r, double flash_ratio) const
{
  std::vector<double> rrXs;
  double rr_hypoX, rr_hypoXWgt;
  for(const auto& rrF : fRRFits){
    if(rrF.min < flash_r && flash_r < rrF.max){
      try{
        rrXs.push_back(rrF.f->GetX(flash_r, 0., fDriftDistance, kEps));
      }catch (...) {
        mf::LogWarning("FlashPredict")
          << "Couldn't find root with fRRFits.\n"
          << "min/_flash_ratio/max:"
          << rrF.min << "/" << flash_r << "/" << rrF.max;
      }
    }
  }
  if(rrXs.size() > 1){//between: [l,h], [l,m], or [h,m]
    rr_hypoX = (rrXs[0] + rrXs[1])/2.;
    double half_interval = (rrXs[1] - rrXs[0])/2.;
    rr_hypoXWgt =  1./(half_interval*half_interval);
  }
  else if(rrXs.size() == 0){//can't estimate
    rr_hypoX = 0.;
    rr_hypoXWgt = 0.;
  }
  else{//(rrXs.size() == 1)
    if(flash_r < fRRFits[2].min){//between: [l, m)
      rr_hypoX =  rrXs[0]/2.;
      double half_interval = rrXs[0]/2.;
      rr_hypoXWgt = 1./(half_interval*half_interval);
    }
    else{//between: (m, h]
      rr_hypoX = (rrXs[0] + fDriftDistance)/2.;
      double half_interval = (fDriftDistance - rrXs[0]);
      rr_hypoXWgt = 1./(half_interval*half_interval);
    }
  }

  if(!fUseUncoatedPMT && !fUseOppVolMetric)
    return {rr_hypoX, rr_hypoXWgt, rr_hypoX, 0.};

  std::vector<double> ratioXs;
  double ratio_hypoX, ratio_hypoXWgt;
  for(const auto& ratioF : fRatioFits){
    if(ratioF.min < flash_ratio && flash_ratio < ratioF.max){
      try{
        ratioXs.push_back(ratioF.f->GetX(flash_ratio, 0., fDriftDistance, kEps));
      }catch (...) {
        mf::LogWarning("FlashPredict")
          << "Couldn't find root with fRatioFits.\n"
          << "min/flash_ratio/max:"
          << ratioF.min << "/" << flash_ratio << "/" << ratioF.max;
      }
    }
  }
  if(ratioXs.size() > 1){//between: [l,h], [l,m], or [h,m]
    ratio_hypoX = (ratioXs[0] + ratioXs[1])/2.;
    double half_interval = (ratioXs[1] - ratioXs[0])/2.;
    ratio_hypoXWgt =  1./(half_interval*half_interval);
  }
  else if(ratioXs.size() == 0){//can't estimate
    ratio_hypoX = 0.;
    ratio_hypoXWgt = 0.;
  }
  else{//(ratioXs.size() == 1)
    if(flash_r < fRatioFits[2].min){//between: [l, m)
      ratio_hypoX =  ratioXs[0]/2.;
      double half_interval = ratioXs[0]/2.;
      ratio_hypoXWgt = 1./(half_interval*half_interval);
    }
    else{//between: (m, h]
      ratio_hypoX = (ratioXs[0] + fDriftDistance)/2.;
      double half_interval = (fDriftDistance - ratioXs[0]);
      ratio_hypoXWgt = 1./(half_interval*half_interval);
    }
  }

  double sum_weights = rr_hypoXWgt + ratio_hypoXWgt;
  double hypo_x =
    (rr_hypoX*rr_hypoXWgt + ratio_hypoX*ratio_hypoXWgt) / sum_weights;
  double hypo_x_err = std::sqrt(sum_weights) / sum_weights;
  return {hypo_x, hypo_x_err, rr_hypoX, ratio_hypoX};
}


// ::flashmatch::Flash_t FlashPredict::GetFlashPESpectrum(const recob::OpFlash& opflash)
// {
//   // prepare container to store flash
//   ::flashmatch::Flash_t flash;
//   flash.time = opflash.Time();
//   // geometry service
//   const art::ServiceHandle<geo::Geometry> geometry;
//   uint nOpDets(geometry->NOpDets());
//   std::vector<double> PEspectrum;
//   PEspectrum.resize(nOpDets);
//   // apply gain to OpDets
//   for (uint OpChannel=0; OpChannel<nOpDets; ++OpChannel) {
//     uint opdet = geometry->OpDetFromOpChannel(OpChannel);
//     PEspectrum[opdet] = opflash.PEs().at(OpChannel);
//   }
//   _pe_reco_v = PEspectrum;

//   // Reset variables
//   flash.x = flash.y = flash.z = 0;
//   flash.x_err = flash.y_err = flash.z_err = 0;
//   double totalPE = 0.;
//   double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;
//   for (unsigned int opdet=0; opdet<PEspectrum.size(); opdet++) {
//     double PMTxyz[3];
//     geometry->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);
//     // Add up the position, weighting with PEs
//     sumy    += PEspectrum[opdet] * PMTxyz[1];
//     sumy2   += PEspectrum[opdet] * PMTxyz[1] * PMTxyz[1];
//     sumz    += PEspectrum[opdet] * PMTxyz[2];
//     sumz2   += PEspectrum[opdet] * PMTxyz[2] * PMTxyz[2];
//     totalPE += PEspectrum[opdet];
//   }

//   flash.y = sumy / totalPE;
//   flash.z = sumz / totalPE;
//   // This is just sqrt(<x^2> - <x>^2)
//   if ((sumy2 * totalPE - sumy * sumy) > 0.)
//     flash.y_err = std::sqrt(sumy2 * totalPE - sumy * sumy) / totalPE;
//   if ((sumz2 * totalPE - sumz * sumz) > 0.)
//     flash.z_err = std::sqrt(sumz2 * totalPE - sumz * sumz) / totalPE;
//   // Set the flash properties
//   flash.pe_v.resize(nOpDets);
//   flash.pe_err_v.resize(nOpDets);
//   // Fill the flash with the PE spectrum
//   for (unsigned int i=0; i<nOpDets; ++i) {
//     const auto PE(PEspectrum.at(i));
//     flash.pe_v.at(i) = PE;
//     flash.pe_err_v.at(i) = std::sqrt(PE);
//   }
//   if (flash.pe_v.size() != nOpDets)
//     throw cet::exception("FlashNeutrinoId") << "Number of channels in beam flash doesn't match the number of OpDets!" << std::endl;
//   return flash;

// }// ::flashmatch::Flash_t FlashPredict::GetFlashPESpectrum


// void FlashPredict::CollectDownstreamPFParticles(
//   const lar_pandora::PFParticleMap& pfParticleMap,
//   const art::Ptr<recob::PFParticle>& particle,
//   lar_pandora::PFParticleVector& downstreamPFParticles) const
// {
//   if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end()){
//     downstreamPFParticles.push_back(particle);
//   }
//   for (const auto &daughterId : particle->Daughters()) {
//     const auto iter(pfParticleMap.find(daughterId));
//     if (iter == pfParticleMap.end()){
//       throw cet::exception("FlashNeutrinoId") << "Scrambled PFParticle IDs" << std::endl;
//     }
//     this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
//   }
// } // void FlashPredict::CollectDownstreamPFParticles


// void FlashPredict::CollectDownstreamPFParticles(
//   const lar_pandora::PFParticleMap& pfParticleMap,
//   const lar_pandora::PFParticleVector& parentPFParticles,
//   lar_pandora::PFParticleVector& downstreamPFParticles) const
// {
//   for (const auto &particle : parentPFParticles){
//     this->CollectDownstreamPFParticles(pfParticleMap, particle, downstreamPFParticles);
//   }
// } // void FlashPredict::CollectDownstreamPFParticles


FlashPredict::ChargeDigestMap FlashPredict::makeChargeDigest(
  const art::Event& evt,
  const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h)
{
  // grab spacepoints associated with PFParticles
  const art::FindManyP<recob::SpacePoint>
    pfp_spacepoints_assns(pfps_h, evt, fPandoraProducer);
  const auto& spacepoints_h =
    evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointProducer);
  const art::FindManyP<recob::Hit>
    spacepoint_hits_assns(spacepoints_h, evt, fSpacePointProducer);
  // grab tracks associated with PFParticles
  // auto const& track_h = evt.getValidHandle<std::vector<recob::Track> >(fTrackProducer);
  // art::FindManyP<recob::Track> pfp_track_assn_v(track_h, evt, fTrackProducer);
  // grab calorimetry info for tracks
  // auto const& calo_h =
  //   evt.getValidHandle<std::vector<anab::Calorimetry> >(fCaloProducer);
  // art::FindManyP<anab::Calorimetry>  track_calo_assn_v(calo_h,
  //                                                      evt, fCaloProducer);

  std::unordered_map<size_t, size_t> pfpMap;
  for(size_t pId=0; pId<pfps_h->size(); pId++) {
    pfpMap[pfps_h->at(pId).Self()] = pId;
  }

  ChargeDigestMap chargeDigestMap;
  // Loop over pandora pfp particles
  for(size_t pId=0; pId<pfps_h->size(); pId++) {
    if(!pfps_h->at(pId).IsPrimary()) continue;
    const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
    unsigned pfpPDGC = std::abs(pfp_ptr->PdgCode());
    if(fSelectNeutrino &&
       (pfpPDGC != 12) && (pfpPDGC != 14) && (pfpPDGC != 16)) continue;
    bk.pfp_to_score++;
    std::set<unsigned> tpcWithHits;

    std::vector<art::Ptr<recob::PFParticle>> particles_in_slice;
    addDaughters(pfpMap, pfp_ptr, pfps_h, particles_in_slice);

    double sliceQ = 0.;
    flashmatch::QCluster_t particlesClusters;
    for(const auto& particle: particles_in_slice) {
      const auto particle_key = particle.key();
      const auto& particle_spacepoints = pfp_spacepoints_assns.at(particle_key);
      double particleQ = 0.;
      flashmatch::QCluster_t spsClusters;
      for(const auto& spacepoint : particle_spacepoints) {
        const auto spacepoint_key = spacepoint.key();
        const auto& hits = spacepoint_hits_assns.at(spacepoint_key);
        const auto& pos = spacepoint->XYZ();
        double spacepointQ = 0.;
        flashmatch::QCluster_t hitsClusters;
        for(const auto& hit : hits) {
          geo::WireID wId = hit->WireID();
          if(fOnlyCollectionWires &&
             fGeometry->SignalType(wId) != geo::kCollection) continue;
          const double hitQ = hit->Integral();
          if(hitQ < fMinHitQ) continue;
          spacepointQ += hitQ;
          hitsClusters.emplace_back(pos[0], pos[1], pos[2], hitQ);
          const auto itpc = wId.TPC;
          tpcWithHits.insert(itpc);
        } // for hits associated to spacepoint
        if(spacepointQ < fMinSpacePointQ) continue;
        particleQ += spacepointQ;
        spsClusters.insert(spsClusters.end(),
                           hitsClusters.begin(), hitsClusters.end());
      } // for spacepoints in particle
      double chargeToNPhots = lar_pandora::LArPandoraHelper::IsTrack(particle) ?
        fChargeToNPhotonsTrack : fChargeToNPhotonsShower;
      particleQ *= chargeToNPhots;
      if(particleQ < fMinParticleQ) continue;
      sliceQ += particleQ;
      particlesClusters.insert(particlesClusters.end(),
                               spsClusters.begin(), spsClusters.end());
    } // for particles in slice
    if(sliceQ < fMinSliceQ) continue;
    chargeDigestMap[sliceQ] =
      ChargeDigest(pId, pfpPDGC, pfp_ptr, particlesClusters, tpcWithHits);
  } // over all slices
  return chargeDigestMap;
}


void FlashPredict::addDaughters(
  const std::unordered_map<size_t, size_t>& pfpMap,
  const art::Ptr<recob::PFParticle>& pfp_ptr,
  const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h,
  std::vector<art::Ptr<recob::PFParticle>>& mothers_daughters) const
{
  mothers_daughters.push_back(pfp_ptr);
  const auto daughters = pfp_ptr->Daughters();
  for(const auto daughter : daughters) {
    const art::Ptr<recob::PFParticle> daughter_ptr(pfps_h, pfpMap.at(daughter));
    addDaughters(pfpMap, daughter_ptr, pfps_h, mothers_daughters);
  }
  return;
}


inline
void FlashPredict::updateFlashMetrics(const FlashMetrics& flashMetrics)
{
  const auto& f = flashMetrics;
  _flash_x = f.x; _flash_x_gl = f.x_gl; _flash_y = f.y;
  _flash_z = f.z; _flash_r = f.r; _flash_pe = f.pe; _flash_unpe = f.unpe;
  _flash_ratio = f.ratio; _flash_time = f.time; _hypo_x = f.hypo;
  _hypo_x_rr  = f.hypo_rr; _hypo_x_ratio = f.hypo_ratio;
}


inline
void FlashPredict::storeFlashMetrics(
  const unsigned flashUId,
  std::map<unsigned, FlashMetrics>& flashMetricsMap) const
{
  flashMetricsMap[flashUId] = FlashMetrics(
    _flash_x, _flash_x_gl, _flash_y, _flash_z, _flash_r, _flash_pe, _flash_unpe,
    _flash_ratio, _flash_time, _hypo_x, _hypo_x_rr, _hypo_x_ratio, true);
}


inline
double FlashPredict::scoreTerm(const double m, const double n,
                               const double mean, const double spread) const
{
  return std::abs(std::abs(m - n) - mean) / spread;
}


inline
double FlashPredict::scoreTerm(const double m,
                               const double mean, const double spread) const
{
  return std::abs(m - mean) / spread;
}


inline
bool FlashPredict::pfpNeutrinoOnEvent(
  const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h) const
{
  for (auto const& p : (*pfps_h)) {
    unsigned pfpPDGC = std::abs(p.PdgCode());
    if ((pfpPDGC == 12) ||
        (pfpPDGC == 14) ||
        (pfpPDGC == 16)) {
      return true;
    }
  }
  return false;
}


void FlashPredict::copyOpHitsInBeamWindow(
  std::vector<recob::OpHit>& opHits,
  const art::Handle<std::vector<recob::OpHit>>& ophits_h) const
{
  double s = fBeamWindowStart;
  double e = fBeamWindowEnd;
  double m = fMinOpHPE;
  // copy ophits that are inside the time window and with PEs
  auto opHitInWindow =
    [s, e, m, this](const recob::OpHit& oph)-> bool
    {return ((oph.PeakTime() > s) &&
             (oph.PeakTime() < e) &&
             (oph.PE() > m) &&
             isPDInCryo(oph.OpChannel())); };
  auto it = std::copy_if(ophits_h->begin(), ophits_h->end(), opHits.begin(),
                         opHitInWindow);
  opHits.resize(std::distance(opHits.begin(), it));
}


//SBND overload
std::vector<FlashPredict::SimpleFlash> FlashPredict::makeSimpleFlashes(
  std::vector<recob::OpHit>& opHits,
  std::vector<recob::OpHit>& opHitsLeft,
  std::vector<recob::OpHit>& opHitsRght) const
{
  std::unique_ptr<TH1D> opHitsTimeHist = std::make_unique<TH1D>(
    "opHitsTimeHist", "ophittime", fTimeBins, fBeamWindowStart, fBeamWindowEnd);
  opHitsTimeHist->SetOption("HIST");
  opHitsTimeHist->SetDirectory(0);//turn off ROOT's object ownership
  std::unique_ptr<TH1D> opHitsTimeHistLeft = std::make_unique<TH1D>(
    "opHitsTimeHistLeft", "ophittimel", fTimeBins, fBeamWindowStart, fBeamWindowEnd);
  opHitsTimeHistLeft->SetOption("HIST");
  opHitsTimeHistLeft->SetDirectory(0);//turn off ROOT's object ownership
  std::unique_ptr<TH1D> opHitsTimeHistRght = std::make_unique<TH1D>(
    "opHitsTimeHistRght", "ophittimer", fTimeBins, fBeamWindowStart, fBeamWindowEnd);
  opHitsTimeHistRght->SetOption("HIST");
  opHitsTimeHistRght->SetDirectory(0);//turn off ROOT's object ownership
  if(!createOpHitsTimeHist(
       opHits, opHitsLeft, opHitsRght,
       opHitsTimeHist, opHitsTimeHistLeft, opHitsTimeHistRght)) return {};

  std::vector<FlashPredict::SimpleFlash> simpleFlashes;
  if(!findSimpleFlashes(simpleFlashes, opHits, kBothOpHs, opHitsTimeHist))
    return {};
  if(opHitsLeft.size() > 0 && opHitsTimeHistLeft->GetEntries() > 0)
    findSimpleFlashes(simpleFlashes, opHitsLeft, kLeftOpHs, opHitsTimeHistLeft);
  if(opHitsRght.size() > 0 && opHitsTimeHistRght->GetEntries() > 0)
    findSimpleFlashes(simpleFlashes, opHitsRght, kRghtOpHs, opHitsTimeHistRght);
  return simpleFlashes;
}


//ICARUS overload
std::vector<FlashPredict::SimpleFlash> FlashPredict::makeSimpleFlashes(
  std::vector<recob::OpHit>& opHits) const
{
  std::unique_ptr<TH1D> opHitsTimeHist = std::make_unique<TH1D>(
    "opHitsTimeHist", "ophittime", fTimeBins, fBeamWindowStart, fBeamWindowEnd);
  opHitsTimeHist->SetOption("HIST");
  opHitsTimeHist->SetDirectory(0);//turn off ROOT's object ownership
  if(!createOpHitsTimeHist(opHits, opHitsTimeHist)) return {};

  std::vector<FlashPredict::SimpleFlash> simpleFlashes;
  if(!findSimpleFlashes(simpleFlashes, opHits, kBothOpHs, opHitsTimeHist))
    return {};
  return simpleFlashes;
}


//SBND overload
bool FlashPredict::createOpHitsTimeHist(
  const std::vector<recob::OpHit>& opHits,
  std::vector<recob::OpHit>& opHitsLeft,
  std::vector<recob::OpHit>& opHitsRght,
  std::unique_ptr<TH1D>& opHitsTimeHist,
  std::unique_ptr<TH1D>& opHitsTimeHistLeft,
  std::unique_ptr<TH1D>& opHitsTimeHistRght) const
{
  for(const auto& oph : opHits) {
    auto ch = oph.OpChannel();
    if(!fUseUncoatedPMT &&
       fPDMapAlgPtr->isPDType(ch, "pmt_uncoated")) continue;
    opHitsTimeHist->Fill(oph.PeakTime(), oph.PE());
    if(sbndPDinTPC(ch) == 0){
      // opHitsLeft.push_back(oph);
      opHitsLeft.emplace_back(oph);
      opHitsTimeHistLeft->Fill(oph.PeakTime(), oph.PE());
    }
    else{
      // opHitsRght.push_back(oph);
      opHitsRght.emplace_back(oph);
      opHitsTimeHistRght->Fill(oph.PeakTime(), oph.PE());
    }
  }
  if(opHitsTimeHist->GetEntries() <= 0 ||
     opHitsTimeHist->Integral() < fMinFlashPE) return false;
  if(opHitsTimeHistLeft->GetEntries() <= 0 ||
     opHitsTimeHistLeft->Integral() < fMinFlashPE)
    opHitsTimeHistLeft->Reset();
  if(opHitsTimeHistRght->GetEntries() <= 0 ||
     opHitsTimeHistRght->Integral() < fMinFlashPE)
    opHitsTimeHistRght->Reset();
  return true;
}


//ICARUS overload
bool FlashPredict::createOpHitsTimeHist(
  const std::vector<recob::OpHit>& opHits,
  std::unique_ptr<TH1D>& opHitsTimeHist) const
{
  for(auto const& oph : opHits) {
    auto ch = oph.OpChannel();
    auto opDetXYZ = fGeometry->OpDetGeoFromOpChannel(ch).GetCenter();
    if(!fGeoCryo->ContainsPosition(opDetXYZ)) continue;
    opHitsTimeHist->Fill(oph.PeakTime(), oph.PE());
  }
  if(opHitsTimeHist->GetEntries() <= 0 ||
     opHitsTimeHist->Integral() < fMinFlashPE) return false;
  return true;
}


bool FlashPredict::findSimpleFlashes(
  std::vector<FlashPredict::SimpleFlash>& simpleFlashes,
  std::vector<recob::OpHit>& opHits,
  const unsigned volumeId,
  std::unique_ptr<TH1D>& opHitsTimeHist) const
{
  OpHitIt opH_beg = opHits.begin();
  for(unsigned flashId=0; flashId<fMaxFlashes; ++flashId){
    int ibin = opHitsTimeHist->GetMaximumBin();
    double maxpeak_time = opHitsTimeHist->GetBinCenter(ibin);
    double lowedge  = maxpeak_time + fFlashStart;
    double highedge = maxpeak_time + fFlashEnd;
    mf::LogDebug("FlashPredict")
      << "light window " << lowedge << " " << highedge << std::endl;
    int lowedge_bin = opHitsTimeHist->FindBin(lowedge);
    int highedge_bin = opHitsTimeHist->FindBin(highedge);
    // check if flash has enough PEs, return if is the first one
    if (opHitsTimeHist->Integral(lowedge_bin, highedge_bin) < fMinFlashPE){
      if(flashId == 0) return false;
      break;
    }
    // clear this peak to enforce non-overlapping flashes
    for(int i=lowedge_bin; i<highedge_bin; ++i){
      opHitsTimeHist->SetBinContent(i, 0.);
    }
    auto peakInsideEdges =
      [&lowedge, &highedge](const recob::OpHit& oph)-> bool
      { return ((lowedge <= oph.PeakTime()) && (oph.PeakTime() <= highedge)); };
    // partition container to move the hits of the flash
    // the iterators point to the boundaries of the partition
    OpHitIt opH_end = std::partition(opH_beg, opHits.end(),
                                     peakInsideEdges);
    simpleFlashes.emplace_back
      (SimpleFlash(flashId, flashId + volumeId,
                   opH_beg, opH_end, maxpeak_time));
    opH_beg = opH_end;
  }
  return true;
}


inline std::string FlashPredict::detectorName(const std::string detName) const
{
  if(detName.find("sbnd")   != std::string::npos) return "SBND";
  if(detName.find("icarus") != std::string::npos) return "ICARUS";
  return "";
}


bool FlashPredict::isPDInCryo(const int pdChannel) const
{
  if(fSBND) return true;
  else { // fICARUS
    // BUG: I believe this function is not working, every now and then
    // I get ophits from the other cryo
    auto& p = fGeometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
    // if the channel is in the Cryostat is relevant
    return fGeoCryo->ContainsPosition(p);
  }
  return false;
}


// TODO: find better, less hacky solution
unsigned FlashPredict::icarusPDinTPC(const int pdChannel) const
{
  auto p = fGeometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  if(fCryostat == 0) p.SetX((p.X() + 222.)/2. - 222.);//OpDets are outside the TPCs
  if(fCryostat == 1) p.SetX((p.X() - 222.)/2. + 222.);//OpDets are outside the TPCs
  return (fGeometry->PositionToTPCID(p)).TPC;
}


bool FlashPredict::isSBNDPDRelevant(const int pdChannel,
                                    const std::set<unsigned>& tpcWithHits) const
{
  // if there's hits on all TPCs all channels are relevant
  if(tpcWithHits.size() == fNTPC) return true;
  unsigned pdTPC = sbndPDinTPC(pdChannel);
  for(auto itpc: tpcWithHits) if(itpc == pdTPC) return true;
  return false;
}


// TODO: find better, less hacky solution
unsigned FlashPredict::sbndPDinTPC(const int pdChannel) const
{
  auto p = fGeometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  p.SetX(p.X()/2.);//OpDets are outside the TPCs
  return (fGeometry->PositionToTPCID(p)).TPC;
}


double FlashPredict::flashXGl(const double hypo_x,
                              const double flash_x) const
{
  double min = 10000;
  unsigned wId = 0;
  auto wIt = fWiresX_gl.begin();
  auto w = wIt;
  for(size_t i=0; i<fWiresX_gl.size(); i++){
    if(((*w<0) == (flash_x<0)) && std::abs(*w - flash_x) < min) {
      wIt = w;
      wId = i;
      min = std::abs(*w - flash_x);
    }
    w++;;
  }
  return (wId % 2) ? (*wIt - hypo_x) : (*wIt + hypo_x);
}


double FlashPredict::driftDistance(const double x) const
{
  auto wit = fWiresX_gl.begin();
  for(size_t i=0; i<fWiresX_gl.size()-1; i++){
    double wxl = *wit;
    wit++;
    double wxh = *wit;
    if(wxl < x && x<= wxh){
      double mid = (wxl + wxh)/2.;
      return (x<=mid) ? std::abs(x-wxl) : std::abs(x-wxh);
    }
  }
  return -10.;
}


unsigned FlashPredict::driftVolume(const double x) const
{
  auto wit = fWiresX_gl.begin();
  for(size_t i=0; i<fWiresX_gl.size()-1; i++){
    double wxl = *wit;
    wit++;
    double wxh = *wit;
    if(wxl < x && x<= wxh){
      double mid = (wxl + wxh)/2.;
      return (x<=mid) ? i*fCryostat : (i*fCryostat)+1;
    }
  }
  return 100;
}


// // TODO: no hardcoding
// // TODO: collapse with the next
// bool FlashPredict::isPDInCryoTPC(double pd_x, size_t itpc)
// {
//   // check whether this optical detector views the light inside this tpc.
//   std::ostringstream lostPDMessage;
//   lostPDMessage << "\nThere's a " << fDetector << " photo detector that belongs nowhere. \n"
//                 << "icryo: " << fCryostat << "\n"
//                 << "itpc:  " << itpc <<  "\n"
//                 << "pd_x:  " << pd_x <<  std::endl;
//   if (fICARUS) {
//     if (fCryostat == 0) {
//       if (itpc == 0 && -400. < pd_x && pd_x < -300. ) return true;
//       else if (itpc == 1 && -100. < pd_x && pd_x < 0.) return true;
//       else {std::cout << lostPDMessage.str(); return false;}
//     }
//     else if (fCryostat == 1) {
//       if (itpc == 0 && 0. < pd_x && pd_x < 100.) return true;
//       else if (itpc == 1 && 300. < pd_x && pd_x < 400.) return true;
//       else {std::cout << lostPDMessage.str(); return false;}
//     }
//   }
//   else if (fSBND) {
//     if (itpc == 0 && -213. < pd_x && pd_x < 0.) return true;
//     else if (itpc == 1 && 0. < pd_x && pd_x < 213.) return true;
//     else {std::cout << lostPDMessage.str(); return false;}
//   }
//   return false;
// }

// bool FlashPredict::isPDInCryoTPC(int pdChannel, size_t itpc)
// {
//   // check whether this optical detector views the light inside this tpc.
//   auto p = fGeometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
//   return isPDInCryoTPC(p.X(), itpc);
// }

// // TODO: no hardcoding
// // TODO: collapse with the previous
// // TODO: figure out what to do with the charge that falls into the crevices
// bool FlashPredict::isChargeInCryoTPC(double qp_x, int icryo, int itpc)
// {
//   std::ostringstream lostChargeMessage;
//   lostChargeMessage << "\nThere's " << fDetector << " charge that belongs nowhere. \n"
//                     << "icryo: " << fCryostat << "\n"
//                     << "itpc: "  << itpc << "\n"
//                     << "qp_x: " << qp_x << std::endl;

//   if (fICARUS) {
//     if (icryo == 0) {
//       if (itpc == 0 && -368.49 <= qp_x && qp_x <= -220.29 ) return true;
//       else if (itpc == 1 && -220.14 <= qp_x && qp_x <= -71.94) return true;
//       // else {std::cout << lostChargeMessage.str(); return false;}
//     }
//     else if (icryo == 1) {
//       if (itpc == 0 && 71.94 <= qp_x && qp_x <= 220.14) return true;
//       else if (itpc == 1 && 220.29 <= qp_x && qp_x <= 368.49) return true;
//       // else {std::cout << lostChargeMessage.str(); return false;}
//     }
//   }
//   else if (fSBND) {
//     if ((itpc == 0 && qp_x < 0) || (itpc == 1 && qp_x > 0) ) return true;
//     else {
//       return false;
//     }
//     //    else {std::cout << lostChargeMessage.str(); return false;}
//   }
//   return false;
// }


template <typename Stream>
void FlashPredict::printBookKeeping(Stream&& out)
{
  std::ostringstream m;
  m << "Bookkeeping\n";
  m << "----------------------------------------\n"
    << "Job Tally\n"
    << "\tEvents:       \t  " << bk.events << "\n";
  if(bk.nopfpneutrino) {
    m << "\tNo PFP Neutrino:  \t -" << bk.nopfpneutrino
      << ", scored as: " << kNoPFPInEvt << "\n";
  }
  // if(bk.noslice) {
  //   m << "\tNo Slice:         \t -" << bk.noslice
  //     << ", scored as: " << kNoSlcInEvt << "\n";
  // }
  if(bk.nonvalidophit) {
    m << "\tNon Valid OpHits: \t -" << bk.nonvalidophit
      << ", scored as: " << kNoOpHInEvt << "\n";
  }
  if(bk.nullophittime) {
    m << "\tNo OpHits in-time:\t -" << bk.nullophittime
      << ", scored as: " << kNoOpHInEvt << "\n";
  }
  m << "\t----------------------\n";
  if(bk.job_bookkeeping != bk.events_processed)
    m << "\tJob Bookkeeping:  \t" << bk.job_bookkeeping << " ERROR!\n";
  m << "\tEvents Processed: \t" << bk.events_processed << "\n"
    << "----------------------------------------\n"
    << "PFP Tally\n"
    << "\tPFP to Score: \t  " << bk.pfp_to_score << "\n";
  if(bk.no_oph_hits) {
    m << "\tHits w/o OpHits:\t -" << bk.no_oph_hits
      << ", scored as: " << kQNoOpHScr << "\n";
  }
  if(bk.no_charge) {
    m << "\tNo Charge:      \t -" << bk.no_charge
      << ", scored as: " << kNoChrgScr << "\n";
  }
  if(bk.no_flash_pe) {
    m << "\tNo VUV PE:      \t -" << bk.no_flash_pe
      << ", scored as: " << k0VUVPEScr << " ERROR!\n";
  }
  m << "\t----------------------\n";
  if(bk.pfp_bookkeeping != bk.scored_pfp)
    m << "\tPFP Bookkeeping: \t" << bk.pfp_bookkeeping << " ERROR!\n";
  m << "\tScored PFP:      \t" << bk.scored_pfp << "\n"
    << "----------------------------------------";
  out << m.str();
}


void FlashPredict::updateBookKeeping()
{
  // account for the reasons that an event could lack
  bk.job_bookkeeping = bk.events
    - bk.nopfpneutrino //- bk.noslice
    - bk.nonvalidophit - bk.nullophittime;

  // account for the reasons that a particle might lack
  bk.pfp_bookkeeping = bk.pfp_to_score
    - bk.no_oph_hits - bk.no_charge
    - bk.no_flash_pe;

  if(bk.events_processed != bk.job_bookkeeping ||
     bk.scored_pfp != bk.pfp_bookkeeping)
    printBookKeeping(mf::LogWarning("FlashPredict"));
}


template <typename Stream>
void FlashPredict::printMetrics(const std::string metric,
                                const int pdgc,
                                const std::set<unsigned>& tpcWithHits,
                                const double term,
                                Stream&& out) const
{
  int isl = int(fNBins * (_charge_x / fDriftDistance));
  std::string tpcs;
  for(auto itpc: tpcWithHits) tpcs += std::to_string(itpc) + ' ';
  out
    << "Big term " << metric << ":\t" << term << "\n"
    << std::left << std::setw(12) << std::setfill(' ')
    << "isl:        \t" << isl << "\n"
    << "pfp.PdgCode:\t" << pdgc << "\n"
    << "tpcWithHits:\t" << tpcs << "\n"
    // << "_slices:    \t" << std::setw(8) << _slices     << "\n"
    << "_flash_time:\t" << std::setw(8) << _flash_time << "\n"
    << "_charge_q:  \t" << std::setw(8) << _charge_q   << "\n"
    << "_flash_pe:  \t" << std::setw(8) << _flash_pe   << ",\t"
    << "_flash_unpe:\t" << std::setw(8) << _flash_unpe << "\n"
    << "_hypo_x:    \t" << std::setw(8) << _hypo_x     << ",\t"
    << "_charge_x:  \t" << std::setw(8) << _charge_x   << "\n"
    << "_flash_x:   \t" << std::setw(8) << _flash_x    << ",\t"
    << "_charge_x_gl\t" << std::setw(8) << _charge_x_gl<< "\n"
    << "_flash_y:   \t" << std::setw(8) << _flash_y    << ",\t"
    << "_charge_y:  \t" << std::setw(8) << _charge_y   << "\n"
    << "_flash_z:   \t" << std::setw(8) << _flash_z    << ",\t"
    << "_charge_z:  \t" << std::setw(8) << _charge_z   << "\n"
    << "_flash_r:   \t" << std::setw(8) << _flash_r    << ",\t"
    << "_flash_ratio\t" << std::setw(8) << _flash_ratio<< "\n";
}


void FlashPredict::beginJob()
{
  bk = BookKeeping();
}

void FlashPredict::endJob()
{
  printBookKeeping(mf::LogWarning("FlashPredict"));
}

DEFINE_ART_MODULE(FlashPredict)
