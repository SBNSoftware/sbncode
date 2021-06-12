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
  , fForceConcurrence(p.get<bool>("ForceConcurrence", false))
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

  // TODO: check params are sane:
  if(fFlashStart > 0. || fFlashEnd < 0.){
    throw cet::exception("FlashPredict");
  }

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
    if(fForceConcurrence || fUseUncoatedPMT || fUseARAPUCAS) {
      throw cet::exception("FlashPredict")
        << "ForceConcurrence: " << std::boolalpha << fForceConcurrence << ",\n"
        << "UseUncoatedPMT:   " << std::boolalpha << fUseUncoatedPMT << ",\n"
        << "UseARAPUCAS:      " << std::boolalpha << fUseARAPUCAS << "\n"
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
  _flash_rr      = -9999.;
  _flash_ratio   = -9999.;
  _hypo_x        = -9999.;
  _hypo_x_err    = -9999.;
  _hypo_x_rr     = -9999.;
  _hypo_x_ratio  = -9999.;
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
      sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
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
  //     sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
  //                             Flash(kNoScrPE), Score(kNoSlcInEvt)));
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
      sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
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


  std::vector<recob::OpHit> opHitsRght, opHitsLeft;
  const std::vector<SimpleFlash> simpleFlashes = (fSBND) ?
    makeSimpleFlashes(opHits, opHitsRght, opHitsLeft) : makeSimpleFlashes(opHits);
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
      sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
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

  // std::map<unsigned, FlashMetrics> flashMetricsMap;
  for(auto& chargeDigest : chargeDigestMap) {
    //const size_t pId = chargeDigest.second.pId;
    const int pfpPDGC = chargeDigest.second.pfpPDGC;
    const auto& pfp_ptr = chargeDigest.second.pfp_ptr;
    const auto& qClusters = chargeDigest.second.qClusters;
    const auto& tpcWithHits = chargeDigest.second.tpcWithHits;

    unsigned hitsInVolume = 0;
    if(fSBND){
      if(tpcWithHits.size() == 2){
        hitsInVolume = kActivityInBoth;
      }
      else{
        if(*tpcWithHits.begin() == kRght) hitsInVolume = kActivityInRght;
        else if(*tpcWithHits.begin() == kLeft) hitsInVolume = kActivityInLeft;
        else {
          mf::LogError("FlashPredict")
            << "ERROR!!! tpcWithHits.size() " << tpcWithHits.size();
        }
      }
    }

    ChargeMetrics charge = computeChargeMetrics(qClusters);
    if(!charge.metric_ok){
      mf::LogWarning("FlashPredict") << "Clusters with No Charge. Skipping...";
      bk.no_charge++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                              Flash(kNoScrPE), Score(kNoChrgScr)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      continue;
    }

    FlashMetrics flash = {};
    Score score = {std::numeric_limits<double>::max()};
    bool hits_ophits_concurrence = false;
    for(auto& simpleFlash : simpleFlashes) {
      unsigned ophsInVolume = simpleFlash.ophsInVolume;
      // TODO still need to check sanity with this, got to think if it
      // could be useful for ICARUS
      if(fSBND){
        if(hitsInVolume != ophsInVolume){
          if(fForceConcurrence) continue;
          else if((hitsInVolume < kActivityInBoth) &&
                  (ophsInVolume < kActivityInBoth)) {
            continue;
          }
        }
      }
      hits_ophits_concurrence = true;

      // TODO: check if metrics already computed, if so skip computing
      // and grab them. Depends on hitsInVolume, ophsInVolume and the
      // current flash number ... perhaps?
      // flashUId = something
      // if(auto it=flashMetricsMap.find(flashUId); it!=flashMetricsMap.end()){
      //   mf::LogDebug("FlashPredict") << "Reusing metrics previously computed, "
      //                                << "for flashUId " << flashUId;
      //   auto flash_tmp = it->second;
      //   flashMetricsMap[flashUId] = flash_tmp;
      // }

      FlashMetrics flash_tmp = computeFlashMetrics(simpleFlash);

      Score score_tmp = computeScore(charge, flash_tmp, tpcWithHits, pfpPDGC);
      if(score_tmp.total > 0. && score_tmp.total < score.total
         && flash_tmp.metric_ok){
        score = score_tmp;
        flash = flash_tmp;
      }
    } // for simpleFlashes
    // TODO: make sure this is ok, again got to think about ICARUS
    if(!hits_ophits_concurrence) {
      std::string extra_message = (fForceConcurrence) ? "" :
        "\nConsider setting ForceConcurrence to false to lower requirements";
      mf::LogWarning("FlashPredict")
        << "No OpHits where there's charge. Skipping..." << extra_message;
      bk.no_oph_hits++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                              Flash(kNoScrPE), Score(kQNoOpHScr)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      continue;
    }

    else if(!flash.metric_ok){
      printMetrics("ERROR", pfpPDGC, tpcWithHits, 0, mf::LogError("FlashPredict"));
      bk.no_flash_pe++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                              Flash(kNoScrPE), Score(k0VUVPEScr)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      continue;
    }

    if(score.total > 0. &&
       score.total < std::numeric_limits<double>::max()){
      if(fMakeTree) {
        updateChargeMetrics(charge);
        updateFlashMetrics(flash);
        updateScore(score);
        _flashmatch_nuslice_tree->Fill();
      }
      bk.scored_pfp++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      Charge c{charge.q, TVector3(charge.x_gl, charge.y, charge.z)};
      Flash f{flash.pe, TVector3(flash.x_gl, flash.y, flash.z)};
      sFM_v->emplace_back(sFM(true, flash.time, c, f, score));
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
  _flashmatch_nuslice_tree->Branch("flash_rr", &_flash_rr, "flash_rr/D");
  _flashmatch_nuslice_tree->Branch("flash_pe", &_flash_pe, "flash_pe/D");
  _flashmatch_nuslice_tree->Branch("flash_unpe", &_flash_unpe, "flash_unpe/D");
  _flashmatch_nuslice_tree->Branch("flash_ratio", &_flash_ratio, "flash_ratio/D");
  _flashmatch_nuslice_tree->Branch("charge_x_gl", &_charge_x_gl, "charge_x_gl/D");
  _flashmatch_nuslice_tree->Branch("charge_x", &_charge_x, "charge_x/D");
  _flashmatch_nuslice_tree->Branch("charge_y", &_charge_y, "charge_y/D");
  _flashmatch_nuslice_tree->Branch("charge_z", &_charge_z, "charge_z/D");
  _flashmatch_nuslice_tree->Branch("charge_q", &_charge_q, "charge_q/D");
  _flashmatch_nuslice_tree->Branch("hypo_x", &_hypo_x, "hypo_x/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_err", &_hypo_x_err, "hypo_x_err/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_rr", &_hypo_x_rr, "hypo_x_rr/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_ratio", &_hypo_x_ratio, "hypo_x_ratio/D");
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
    for (int ib = 1; ib <= bins; ++ib) {
      double me = temphisto->GetBinContent(ib);
      rr_means.push_back(me);
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in rr" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      rr_spreads.push_back(tt);
    }
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
    for (int ib = 1; ib <= bins; ++ib) {
      double me = temphisto->GetBinContent(ib);
      pe_means.push_back(me);
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in pe" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      pe_spreads.push_back(tt);
    }
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


FlashPredict::ChargeMetrics FlashPredict::computeChargeMetrics(
  const flashmatch::QCluster_t& qClusters) const
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
    double charge_q = norm;
    double charge_x_gl = xave / norm;
    double charge_x = driftDistance(charge_x_gl);
    double charge_y = yave / norm;
    double charge_z = zave / norm;
    return {charge_x, charge_x_gl, charge_y, charge_z, charge_q, true};
  }
  return {};
}


FlashPredict::FlashMetrics FlashPredict::computeFlashMetrics(
  const SimpleFlash& simpleFlash) const
{
  const OpHitIt opH_beg = simpleFlash.opH_beg;
  const OpHitIt opH_end = simpleFlash.opH_end;

  double peSumMax_wallX = wallXWithMaxPE(opH_beg, opH_end);
 
  double sum = 0.;
  double sum_PE = 0.;
  double sum_PE2 = 0.;
  double sum_unPE = 0.;
  double sum_visARA_PE = 0.;
  double sum_PE2Y  = 0.; double sum_PE2Z  = 0.;
  double sum_PE2Y2 = 0.; double sum_PE2Z2 = 0.;

  for(auto oph=opH_beg; oph!=opH_end; ++oph){
    int opChannel = oph->OpChannel();
    auto& opDet = fGeometry->OpDetGeoFromOpChannel(opChannel);
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
      if(fUseOppVolMetric &&
         std::abs(peSumMax_wallX-opDetXYZ.X()) > 5.) sum_unPE += ophPE;
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
    double flash_pe    = sum_PE;
    double flash_unpe  = sum_unPE;
    double flash_ratio = fOpDetNormalizer * flash_unpe / flash_pe;
    if(fUseARAPUCAS) {
      flash_unpe  += sum_visARA_PE;
      flash_ratio = (fOpDetNormalizer * sum_unPE  + sum_visARA_PE ) / flash_pe;
    }
    double flash_y  = sum_PE2Y / sum_PE2;
    double flash_z  = sum_PE2Z / sum_PE2;
    double flash_rr = std::sqrt(
      std::abs(sum_PE2Y2 + sum_PE2Z2 + sum_PE2 * (flash_y * flash_y + flash_z * flash_z)
               - 2.0 * (flash_y * sum_PE2Y + flash_z * sum_PE2Z) ) / sum_PE2);
    auto [hypo_x, hypo_x_err, hypo_x_rr, hypo_x_ratio] =
      hypoFlashX_fits(flash_rr, flash_ratio);
    // TODO: using _hypo_x make further corrections to _flash_time to
    // account for light transport time and/or rising edge
    // double flash_time = timeCorrections(simpleFlash.maxpeak_time, hypo_x);
    double flash_time = simpleFlash.maxpeak_time;
    double flash_x = peSumMax_wallX;
    double flash_x_gl = flashXGl(hypo_x, flash_x);

    return {flash_x, flash_x_gl, flash_y, flash_z, flash_rr, flash_pe, flash_unpe,
      flash_ratio, flash_time, hypo_x, hypo_x_err, hypo_x_rr, hypo_x_ratio, true};
  }
  else {
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
    return {};
  }
}


FlashPredict::Score FlashPredict::computeScore(
  const ChargeMetrics& charge,
  const FlashMetrics& flash,
  const std::set<unsigned>& tpcWithHits,
  const int pdgc) const
{
  double score = 0.;
  unsigned tcount = 0;
  int xbin = static_cast<int>(fNBins * (charge.x / fDriftDistance));
  auto out = mf::LogWarning("FlashPredict");

  double scr_y = scoreTerm(flash.y, charge.y, dy_means[xbin], dy_spreads[xbin]);
  if (scr_y > fTermThreshold) printMetrics("Y", pdgc, tpcWithHits, scr_y, out);
  score += scr_y;
  tcount++;
  double scr_z = scoreTerm(flash.z, charge.z, dz_means[xbin], dz_spreads[xbin]);
  if (scr_z > fTermThreshold) printMetrics("Z", pdgc, tpcWithHits, scr_z, out);
  score += scr_z;
  tcount++;
  double scr_rr = scoreTerm(flash.rr, rr_means[xbin], rr_spreads[xbin]);
  if (scr_rr > fTermThreshold) printMetrics("RR", pdgc, tpcWithHits, scr_rr, out);
  score += scr_rr;
  tcount++;
  double scr_ratio = 0.;
  if (fUseUncoatedPMT || fUseOppVolMetric) {
    scr_ratio = scoreTerm(flash.ratio, pe_means[xbin], pe_spreads[xbin]);
    if (scr_ratio > fTermThreshold) printMetrics("RATIO", pdgc, tpcWithHits, scr_ratio, out);
    score += scr_ratio;
    tcount++;
  }
  mf::LogDebug("FlashPredict")
    << "score:\t" << score << "using " << tcount << " terms";
  return {score, scr_y, scr_z, scr_rr, scr_ratio};
}


std::tuple<double, double, double, double> FlashPredict::hypoFlashX_fits(
  double flash_rr, double flash_ratio) const
{
  std::vector<double> rrXs;
  double rr_hypoX, rr_hypoXWgt;
  for(const auto& rrF : fRRFits){
    if(rrF.min < flash_rr && flash_rr < rrF.max){
      try{
        rrXs.emplace_back(rrF.f->GetX(flash_rr, 0., fDriftDistance, kEps));
      }catch (...) {
        mf::LogWarning("FlashPredict")
          << "Couldn't find root with fRRFits.\n"
          << "min/_flash_rr/max:"
          << rrF.min << "/" << flash_rr << "/" << rrF.max;
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
    if(flash_rr < fRRFits[2].min){//between: [l, m)
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
        ratioXs.emplace_back(ratioF.f->GetX(flash_ratio, 0., fDriftDistance, kEps));
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
    if(flash_ratio < fRatioFits[2].min){//between: [l, m)
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
void FlashPredict::updateChargeMetrics(const ChargeMetrics& chargeMetrics)
{
  const auto& c = chargeMetrics;
  _charge_x = c.x; _charge_x_gl = c.x_gl; _charge_y = c.y;
  _charge_z = c.z; _charge_q = c.q;
}


inline
void FlashPredict::updateFlashMetrics(const FlashMetrics& flashMetrics)
{
  const auto& f = flashMetrics;
  _flash_x = f.x; _flash_x_gl = f.x_gl; _flash_y = f.y;
  _flash_z = f.z; _flash_rr = f.rr; _flash_pe = f.pe; _flash_unpe = f.unpe;
  _flash_ratio = f.ratio; _flash_time = f.time; _hypo_x = f.hypo;
  _hypo_x_err = f.hypo_err; _hypo_x_rr  = f.hypo_rr;
  _hypo_x_ratio = f.hypo_ratio;
}


inline
void FlashPredict::updateScore(const Score& score)
{
  _score = score.total, _scr_y = score.y, _scr_z = score.z,
    _scr_rr = score.rr, _scr_ratio = score.ratio;
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
  std::vector<recob::OpHit>& opHitsRght,
  std::vector<recob::OpHit>& opHitsLeft) const
{
  std::unique_ptr<TH1D> opHitsTimeHist = std::make_unique<TH1D>(
    "opHitsTimeHist", "ophittime", fTimeBins, fBeamWindowStart, fBeamWindowEnd);
  opHitsTimeHist->SetOption("HIST");
  opHitsTimeHist->SetDirectory(0);//turn off ROOT's object ownership
  std::unique_ptr<TH1D> opHitsTimeHistRght = std::make_unique<TH1D>(
    "opHitsTimeHistRght", "ophittimer", fTimeBins, fBeamWindowStart, fBeamWindowEnd);
  opHitsTimeHistRght->SetOption("HIST");
  opHitsTimeHistRght->SetDirectory(0);//turn off ROOT's object ownership
  std::unique_ptr<TH1D> opHitsTimeHistLeft = std::make_unique<TH1D>(
    "opHitsTimeHistLeft", "ophittimel", fTimeBins, fBeamWindowStart, fBeamWindowEnd);
  opHitsTimeHistLeft->SetOption("HIST");
  opHitsTimeHistLeft->SetDirectory(0);//turn off ROOT's object ownership
  if(!createOpHitsTimeHist(
       opHits, opHitsRght, opHitsLeft,
       opHitsTimeHist, opHitsTimeHistRght, opHitsTimeHistLeft)) return {};

  bool oph_in_rght = false, oph_in_left = false;
  std::vector<FlashPredict::SimpleFlash> simpleFlashes;
  if(opHitsRght.size() > 0 && opHitsTimeHistRght->GetEntries() > 0){
    oph_in_rght = findSimpleFlashes(simpleFlashes, opHitsRght,
                                    kActivityInRght, opHitsTimeHistRght);
  }
  if(opHitsLeft.size() > 0 && opHitsTimeHistLeft->GetEntries() > 0){
    oph_in_left = findSimpleFlashes(simpleFlashes, opHitsLeft,
                                    kActivityInLeft, opHitsTimeHistLeft);
  }
  if(oph_in_rght == true && oph_in_left == true ){
    findSimpleFlashes(simpleFlashes, opHits, kActivityInBoth, opHitsTimeHist);
  }
  else if(oph_in_rght == false && oph_in_left == false){
    return {};
  }
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
  if(!findSimpleFlashes(simpleFlashes, opHits, kActivityInBoth, opHitsTimeHist))
    return {};
  return simpleFlashes;
}


//SBND overload
bool FlashPredict::createOpHitsTimeHist(
  const std::vector<recob::OpHit>& opHits,
  std::vector<recob::OpHit>& opHitsRght,
  std::vector<recob::OpHit>& opHitsLeft,
  std::unique_ptr<TH1D>& opHitsTimeHist,
  std::unique_ptr<TH1D>& opHitsTimeHistRght,
  std::unique_ptr<TH1D>& opHitsTimeHistLeft) const
{
  for(const auto& oph : opHits) {
    auto ch = oph.OpChannel();
    if(!fUseUncoatedPMT &&
       fPDMapAlgPtr->isPDType(ch, "pmt_uncoated")) continue;
    opHitsTimeHist->Fill(oph.PeakTime(), oph.PE());
    if(sbndPDinTPC(ch) == kRght){
      opHitsRght.emplace_back(oph);
      opHitsTimeHistRght->Fill(oph.PeakTime(), oph.PE());
    }
    else{// sbndPDinTPC(ch) == kLeft
      opHitsLeft.emplace_back(oph);
      opHitsTimeHistLeft->Fill(oph.PeakTime(), oph.PE());
    }
  }
  if(opHitsTimeHist->GetEntries() <= 0 ||
     opHitsTimeHist->Integral() < fMinFlashPE) return false;
  if(opHitsTimeHistRght->GetEntries() <= 0 ||
     opHitsTimeHistRght->Integral() < fMinFlashPE)
    opHitsTimeHistRght->Reset();
  if(opHitsTimeHistLeft->GetEntries() <= 0 ||
     opHitsTimeHistLeft->Integral() < fMinFlashPE)
    opHitsTimeHistLeft->Reset();
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
  const unsigned ophsInVolume,
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
      (SimpleFlash(flashId, ophsInVolume,
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


double FlashPredict::wallXWithMaxPE(const OpHitIt opH_beg,
                                    const OpHitIt opH_end) const
{
  std::map<double, double> opdetX_PE {{-99999., 0.}};
  for(auto oph=opH_beg; oph!=opH_end; ++oph){
    double ophPE = oph->PE();
    double ophPE2 = ophPE*ophPE;
    double opdetX = fGeometry->OpDetGeoFromOpChannel(
      oph->OpChannel()).GetCenter().X();
    bool stored = false;
    for(auto& m : opdetX_PE){
      if(std::abs(m.first - opdetX) < 5.) {
        m.second += ophPE2;
        stored = true;
        break;
      }
    }
    if(!stored) opdetX_PE[opdetX] = ophPE2;
  }
  auto maxIt = std::max_element(
    opdetX_PE.begin(), opdetX_PE.end(),
    [] (const auto& a, const auto& b) ->bool{return a.second < b.second;});
  return maxIt->first;
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
  int xbin = static_cast<int>(fNBins * (_charge_x / fDriftDistance));
  std::string tpcs;
  for(auto itpc: tpcWithHits) tpcs += std::to_string(itpc) + ' ';
  out
    << "Big term " << metric << ":\t" << term << "\n"
    << std::left << std::setw(12) << std::setfill(' ')
    << "xbin:        \t" << xbin << "\n"
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
    << "_flash_rr:  \t" << std::setw(8) << _flash_rr   << ",\t"
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
