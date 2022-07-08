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
  , fBeamSpillTimeStart(p.get<double>("BeamSpillTimeStart")) //us
  , fBeamSpillTimeEnd(p.get<double>("BeamSpillTimeEnd"))// us
  , fFlashFindingTimeStart(p.get<double>("FlashFindingTimeStart")) //us
  , fFlashFindingTimeEnd(p.get<double>("FlashFindingTimeEnd"))// us
  , fFlashStart(p.get<double>("FlashStart")) // in us w.r.t. flash time
  , fFlashEnd(p.get<double>("FlashEnd"))  // in us w.r.t flash time
  , fTimeBins(unsigned(1/fTickPeriod * (fFlashFindingTimeEnd - fFlashFindingTimeStart)))
  , fSelectNeutrino(p.get<bool>("SelectNeutrino", true)) // only attempt to match potential neutrino slices
  , fOnlyCollectionWires(p.get<bool>("OnlyCollectionWires", true))
  , fForceConcurrence(p.get<bool>("ForceConcurrence", false)) // require light and charge to coincide, different requirements for SBND and ICARUS
  , fCorrectDriftDistance(p.get<bool>("CorrectDriftDistance", false)) // require light and charge to coincide, different requirements for SBND and ICARUS
  , fUseUncoatedPMT(p.get<bool>("UseUncoatedPMT", false))
  , fUseOppVolMetric(p.get<bool>("UseOppVolMetric", false))
  , fUseARAPUCAS(p.get<bool>("UseARAPUCAS", false))
  , fStoreTrueNus(p.get<bool>("StoreTrueNus", false))
  , fStoreCheatMCT0(p.get<bool>("StoreCheatMCT0", false))
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
  , fOpHitTime(p.get<std::string>("OpHitTime", "RiseTime"))
  , fUseOpHitRiseTime("RiseTime" == fOpHitTime)
  , fUseOpHitPeakTime("PeakTime" == fOpHitTime)
  , fUseOpHitStartTime("StartTime" == fOpHitTime)
  , fMaxFlashes(p.get<unsigned>("MaxFlashes", 1))
  , fMinOpHPE(p.get<double>("MinOpHPE", 0.0))
  , fMinFlashPE(p.get<double>("MinFlashPE", 0.0))
  , fDetector(detectorName(fGeometry->DetectorName()))
  , fSBND((fDetector == "SBND") ? true : false )
  , fICARUS((fDetector == "ICARUS") ? true : false )
  , fPDMapAlgPtr(art::make_tool<opdet::PDMapAlg>(p.get<fhicl::ParameterSet>("PDMapAlg")))
  , fNTPC(fGeometry->NTPC())
  , fCryostat(p.get<int>("Cryostat", 0)) //set =0 or =1 for ICARUS to match reco chain selection
  , fGeoCryo(std::make_unique<geo::CryostatGeo>(fGeometry->Cryostat(fCryostat)))
  , fWiresX_gl(wiresXGl())
  , fDriftDistance(driftDistance())
  , fXBins(p.get<double>("XBins"))
  , fXBinWidth(fDriftDistance/fXBins)// cm
  , fRR_TF1_fit(p.get<std::string>("rr_TF1_fit", "pol3"))
  , fRatio_TF1_fit(p.get<std::string>("ratio_TF1_fit", "pol3"))
  , fYBins(p.get<unsigned>("YBins", 0.)) // roughly match the rows of opdets
  , fZBins(p.get<unsigned>("ZBins", 0.)) // roughly match the columns of opdets
  , fYLow(p.get<double>("YLow", 0.)) // lowest opdet position in cm
  , fYHigh(p.get<double>("YHigh", 0.)) // highest opdet position in cm
  , fZLow(p.get<double>("ZLow", 0.)) // most upstream opdet position in cm
  , fZHigh(p.get<double>("ZHigh", 0.)) // most downstream opdet position in cm
  , fPolCoefsY(p.get<std::vector<double>>("PolCoefsY", std::vector<double>{3, 0.})) // correcting coefficients for flash Y
  , fPolCoefsZ(p.get<std::vector<double>>("PolCoefsZ", std::vector<double>{3, 0.})) // correcting coefficients for flash Z
  , fSkewLimitY(p.get<double>("SkewLimitY", 10.)) //
  , fSkewLimitZ(p.get<double>("SkewLimitZ", 10.)) //
  , fOpDetNormalizer((fSBND) ? 4 : 1)
  , fTermThreshold(p.get<double>("ThresholdTerm", 30.))
{
  produces< std::vector<sbn::SimpleFlashMatch> >();
  produces< art::Assns <recob::PFParticle, sbn::SimpleFlashMatch> >();

  // TODO: check that
  // fBeamSpillTime(Start/End) and fFlashFindingTime(Start/End) are sane

  // TODO: check that all params are sane
  if(fFlashStart > 0. || fFlashEnd < 0.){
    throw cet::exception("FlashPredict")
      << "fFlashStart has to be non-positive, "
      << "and fFlashEnd has to be non-negative.";
  }

  if(std::abs(p.get<double>("DriftDistance")-driftDistance()) > 0.001){
    mf::LogError("FlashPredict")
      << "Provided driftDistance: " << p.get<double>("DriftDistance")
      << " is different than expected: " << driftDistance();
  }

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

  if (((fOpHitTime != "RiseTime" && fOpHitTime != "PeakTime" &&
        fOpHitTime != "StartTime")) ||
      (fUseOpHitRiseTime + fUseOpHitPeakTime + fUseOpHitStartTime  > 1) ||
      (fUseOpHitRiseTime + fUseOpHitPeakTime + fUseOpHitStartTime <= 0))
    throw cet::exception("FlashPredict")
      << "Supported OpHit times are 'RiseTime', 'PeakTime' or 'StartTime', "
      << "selected OpHitTime parameter is: " << fOpHitTime << std::endl;

  if (fMakeTree) initTree();

  loadMetrics();

  consumes<std::vector<recob::PFParticle>>(fPandoraProducer);
  consumes<std::vector<recob::Slice>>(fPandoraProducer);
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
  if(fMakeTree){
    _evt = evt.event();
    _sub = evt.subRun();
    _run = evt.run();
    _slices        = 0;
    _true_nus      = -9999;
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
    _mcT0          = -9999.;
  }
  bk.events++;

  if(fMakeTree && fStoreTrueNus){
    _true_nus = trueNus(evt);
  }

  // grab PFParticles in event
  const auto pfps_h =
    evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraProducer);
  if (pfps_h->size() == 0) {
      mf::LogWarning("FlashPredict")
        << "No recob:PFParticle on event. Skipping...";
      bk.noslice++;
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

  if(fMakeTree){
    auto const& slice_h = evt.getValidHandle<std::vector<recob::Slice>>(fPandoraProducer);
    _slices = slice_h.product()->size();
    if (_slices == 0) {
      mf::LogWarning("FlashPredict")
        << "No recob:Slice on event. Skipping...";
      bk.noslice++;
      updateBookKeeping();
      for(size_t pId=0; pId<pfps_h->size(); pId++) {
        if(!pfps_h->at(pId).IsPrimary()) continue;
        const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
        sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                                Flash(kNoScrPE), Score(kNoSlcInEvt)));
        util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      }
      evt.put(std::move(sFM_v));
      evt.put(std::move(pfp_sFM_assn_v));
      return;
    }
  }

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
  copyOpHitsInFlashFindingWindow(opHits, ophits_h);

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
      copyOpHitsInFlashFindingWindow(opHitsARA, ophits_ara_h);
      opHits.insert(opHits.end(),
                    opHitsARA.begin(), opHitsARA.end());
    }
  }


  std::vector<recob::OpHit> opHitsRght, opHitsLeft;
  const std::vector<SimpleFlash> simpleFlashes = (fSBND) ?
    makeSimpleFlashes(opHits, opHitsRght, opHitsLeft) : makeSimpleFlashes(opHits);
  if(simpleFlashes.empty() ||
     (simpleFlashes[0].maxpeak_time<fBeamSpillTimeStart || fBeamSpillTimeEnd<simpleFlashes[0].maxpeak_time)){
    mf::LogWarning("FlashPredict")
      << "\nNo OpHits in beam window [" << fBeamSpillTimeStart << ", " << fBeamSpillTimeEnd << "], "
      << "\nor the integral is less or equal to " << fMinFlashPE << " or 0."
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

  ChargeDigestMap chargeDigestMap = makeChargeDigest(evt, pfps_h);

  std::map<unsigned, FlashMetrics> flashMetricsMap;
  for(auto& chargeDigestPair : chargeDigestMap) {
    const auto& chargeDigest = chargeDigestPair.second;
    const auto& pfp_ptr = chargeDigest.pfp_ptr;
    const unsigned hitsInVolume = chargeDigest.hitsInVolume;
    bk.pfp_to_score++;
    if(chargeDigestPair.first < 0.){
      mf::LogDebug("FlashPredict") << "Not a nu candidate slice. Skipping...";
      bk.no_nu_candidate++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                              Flash(kNoScrPE), Score(kNotANuScr)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      continue;
    }

    ChargeMetrics charge = computeChargeMetrics(chargeDigest);
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
      if(hitsInVolume != ophsInVolume){
        if(fSBND){
          if(fForceConcurrence) continue;
          else if((hitsInVolume < kActivityInBoth) &&
                  (ophsInVolume < kActivityInBoth)) {
            continue;
          }
        }
        else if(fICARUS){
          if((hitsInVolume < kActivityInBoth) &&
             (ophsInVolume < kActivityInBoth)) {
            continue;
          }
          else if(fForceConcurrence && hitsInVolume == kActivityInBoth) continue;
        }
      }
      hits_ophits_concurrence = true;

      unsigned flashUId = simpleFlash.ophsInVolume * 10 + simpleFlash.flashId;
      bool mets_in_map = flashMetricsMap.find(flashUId) != flashMetricsMap.end();
      FlashMetrics flash_tmp = (mets_in_map) ?
        flashMetricsMap[flashUId] : computeFlashMetrics(simpleFlash);
      if(mets_in_map){
        mf::LogDebug("FlashPredict")
          << "Reusing metrics previously computed, for flashUId " << flashUId;
      }
      else
        flashMetricsMap[flashUId] = flash_tmp;

      Score score_tmp = computeScore(charge, flash_tmp);
      if(score_tmp.total > 0. && score_tmp.total < score.total
         && flash_tmp.metric_ok){

        // TODO: REMOVE IF BLOCK BELOW BEFORE MERGING
        if (score.total < std::numeric_limits<double>::max()) {
          // // TODO:
          // // HACK: to avoid matching late light
          // double flashes_time_diff = flash_tmp.time - flash.time;
          // double f0_light_charge_ratio = std::log(flash.pe)/std::log(charge.q);
          // double f1_light_charge_ratio = std::log(flash_tmp.pe)/std::log(charge.q);
          // if ((0.05 < flashes_time_diff && flashes_time_diff < 1.6) &&
          //     (0.55<f0_light_charge_ratio && f1_light_charge_ratio<0.55))
          //   continue;
          // // TODO:
          // // HACK: to avoid matching late light
          std::cout << "\n----Updating the score and matched flash"<< "\n";
          std::cout << "_slices:\t" << _slices << "\n";
          std::cout << "_mcT0:\t" << _mcT0 << "\n";
          std::cout << "_petoq:\t" << _petoq << "\n";
          std::cout << "_true_nus:\t" << _true_nus << "\n";
          std::cout << "charge:\n";
          std::cout << charge.dumpMetrics() << "\n";
          std::cout << "score.total:    \t" << score.total << "  \t" << score_tmp.total << "\n";
          std::cout << "  score.y:      \t" << score.y << "  \t" << score_tmp.y << "\n";
          std::cout << "  score.z:      \t" << score.z << "  \t" << score_tmp.z << "\n";
          std::cout << "  score.rr:     \t" << score.rr << "  \t" << score_tmp.rr << "\n";
          std::cout << "  score.ratio:  \t" << score.ratio << "  \t" << score_tmp.ratio << "\n";
          std::cout << "  score.slope:  \t" << score.slope << "  \t" << score_tmp.slope << "\n";
          std::cout << "  score.petoq:  \t" << score.petoq << "  \t" << score_tmp.petoq << "\n";
          std::cout << "----from " << "\n";
          std::cout << "flash:\n"
                    << "  id:      " << flash.id          << "\t" << flash_tmp.id          << "\n"
                    << "  activity:" << flash.activity    << "\t" << flash_tmp.activity    << "\n"
                    << "  x:       " << flash.x           << "\t" << flash_tmp.x           << "\n"
                    << "  h_x:     " << flash.h_x         << "\t" << flash_tmp.h_x         << "\n"
                    << "  x_gl:    " << flash.x_gl        << "\t" << flash_tmp.x_gl        << "\n"
                    << "  y:       " << flash.y           << "\t" << flash_tmp.y           << "\n"
                    << "  z:       " << flash.z           << "\t" << flash_tmp.z           << "\n"
                    << "  yb:      " << flash.yb          << "\t" << flash_tmp.yb          << "\n"
                    << "  zb:      " << flash.zb          << "\t" << flash_tmp.zb          << "\n"
                    << "  xw:      " << flash.xw          << "\t" << flash_tmp.xw          << "\n"
                    << "  yw:      " << flash.yw          << "\t" << flash_tmp.yw          << "\n"
                    << "  zw:      " << flash.zw          << "\t" << flash_tmp.zw          << "\n"
                    << "  rr:      " << flash.rr          << "\t" << flash_tmp.rr          << "\n"
                    << "  pe:      " << flash.pe          << "\t" << flash_tmp.pe          << "\n"
                    << "  unpe:    " << flash.unpe        << "\t" << flash_tmp.unpe        << "\n"
                    << "  ratio:   " << flash.ratio       << "\t" << flash_tmp.ratio       << "\n"
                    << "  slope:   " << flash.slope       << "\t" << flash_tmp.slope       << "\n"
                    << "  time:    " << flash.time        << "\t" << flash_tmp.time        << "\n"
                    << "  h_xerr:  " << flash.h_xerr      << "\t" << flash_tmp.h_xerr      << "\n"
                    << "  h_xrr:   " << flash.h_xrr       << "\t" << flash_tmp.h_xrr       << "\n"
                    << "  h_xratio:" << flash.h_xratio    << "\t" << flash_tmp.h_xratio    << "\n"
                    << "  y_skew:  " << flash.y_skew      << "\t" << flash_tmp.y_skew      << "\n"
                    << "  z_skew:  " << flash.z_skew      << "\t" << flash_tmp.z_skew      << "\n"
                    << "  metric_ok: "<< std::boolalpha << flash.metric_ok << "\t" << flash_tmp.metric_ok  << "\n";
        } // if (score.total < std::numeric_limits<double>::max())
        // TODO: REMOVE ABOVE IF BEFORE MERGING
        score = score_tmp;
        flash = flash_tmp;
        // TODO: update x position if new match is made
        // charge.x = driftCorrection(charge.xb, flash.time);
        // charge.x_gl  = xGlCorrection(charge.x_gl, charge.xb, flash.time);
      }
    } // for simpleFlashes
    if(!hits_ophits_concurrence) {
      std::string extra_message = (!fForceConcurrence) ? "" :
        "\nConsider setting ForceConcurrence to false to lower requirements";
      mf::LogInfo("FlashPredict")
        << "No OpHits where there's charge. Skipping..." << extra_message;
      bk.no_oph_hits++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      sFM_v->emplace_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                              Flash(kNoScrPE), Score(kQNoOpHScr)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      continue;
    }
    else if(!flash.metric_ok){
      printMetrics("ERROR", charge, flash, 0, mf::LogError("FlashPredict"));
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
        _mcT0 = chargeDigest.mcT0;
        _petoq = PEToQ(flash.pe, charge.q);
        updateChargeMetrics(charge);
        updateFlashMetrics(flash);
        updateScore(score);
        _flashmatch_nuslice_tree->Fill();
      }
      bk.scored_pfp++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      Charge c{charge.q, TVector3(charge.x_gl, charge.y, charge.z),
               TVector3(charge.x_glw, charge.yw, charge.zw)};
      Flash f{flash.pe, TVector3(flash.x_gl, flash.y, flash.z),
              TVector3(flash.xw, flash.yw, flash.zw)};
      sFM_v->emplace_back(sFM(true, flash.time, c, f, score));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
    }
    else{
      mf::LogError("FlashPredict")
        << "ERROR: score <= 0. Dumping info.\n"
        << "score:       " << score.total << "\n"
        << "score.y:     " << score.y << "\t"
        << "score.z:     " << score.z << "\n"
        << "score.rr:    " << score.rr << "\t"
        << "score.ratio: " << score.ratio << "\n"
        << "score.slope: " << score.slope << "\t"
        << "score.petoq: " << score.petoq << "\n";
      printMetrics("ERROR", charge, flash, 0, mf::LogError("FlashPredict"));
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
  _flashmatch_nuslice_tree->Branch("slices", &_slices, "slices/I");
  _flashmatch_nuslice_tree->Branch("flash_id", &_flash_id, "flash_id/I");
  _flashmatch_nuslice_tree->Branch("flash_activity", &_flash_activity, "flash_activity/I");
  _flashmatch_nuslice_tree->Branch("flash_x", &_flash_x, "flash_x/D");
  _flashmatch_nuslice_tree->Branch("flash_yb", &_flash_yb, "flash_yb/D");
  _flashmatch_nuslice_tree->Branch("flash_zb", &_flash_zb, "flash_zb/D");
  _flashmatch_nuslice_tree->Branch("flash_x_gl", &_flash_x_gl, "flash_x_gl/D");
  _flashmatch_nuslice_tree->Branch("flash_y", &_flash_y, "flash_y/D");
  _flashmatch_nuslice_tree->Branch("flash_z", &_flash_z, "flash_z/D");
  _flashmatch_nuslice_tree->Branch("flash_xw", &_flash_xw, "flash_xw/D");
  _flashmatch_nuslice_tree->Branch("flash_yw", &_flash_yw, "flash_yw/D");
  _flashmatch_nuslice_tree->Branch("flash_zw", &_flash_zw, "flash_zw/D");
  _flashmatch_nuslice_tree->Branch("flash_rr", &_flash_rr, "flash_rr/D");
  _flashmatch_nuslice_tree->Branch("flash_ratio", &_flash_ratio, "flash_ratio/D");
  _flashmatch_nuslice_tree->Branch("flash_slope", &_flash_slope, "flash_slope/D");
  _flashmatch_nuslice_tree->Branch("flash_pe", &_flash_pe, "flash_pe/D");
  _flashmatch_nuslice_tree->Branch("flash_unpe", &_flash_unpe, "flash_unpe/D");
  _flashmatch_nuslice_tree->Branch("flash_time", &_flash_time, "flash_time/D");
  _flashmatch_nuslice_tree->Branch("hypo_x", &_hypo_x, "hypo_x/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_err", &_hypo_x_err, "hypo_x_err/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_rr", &_hypo_x_rr, "hypo_x_rr/D");
  _flashmatch_nuslice_tree->Branch("hypo_x_ratio", &_hypo_x_ratio, "hypo_x_ratio/D");
  _flashmatch_nuslice_tree->Branch("y_skew", &_y_skew, "y_skew/D");
  _flashmatch_nuslice_tree->Branch("z_skew", &_z_skew, "z_skew/D");
  _flashmatch_nuslice_tree->Branch("y_kurt", &_y_kurt, "y_kurt/D");
  _flashmatch_nuslice_tree->Branch("z_kurt", &_z_kurt, "z_kurt/D");
  _flashmatch_nuslice_tree->Branch("charge_id", &_charge_id, "charge_id/I");
  _flashmatch_nuslice_tree->Branch("charge_activity", &_charge_activity, "charge_activity/I");
  _flashmatch_nuslice_tree->Branch("charge_x", &_charge_x, "charge_x/D");
  _flashmatch_nuslice_tree->Branch("charge_x_gl", &_charge_x_gl, "charge_x_gl/D");
  _flashmatch_nuslice_tree->Branch("charge_y", &_charge_y, "charge_y/D");
  _flashmatch_nuslice_tree->Branch("charge_z", &_charge_z, "charge_z/D");
  _flashmatch_nuslice_tree->Branch("charge_x_glw", &_charge_x_glw, "charge_xglw/D");
  _flashmatch_nuslice_tree->Branch("charge_yw", &_charge_yw, "charge_yw/D");
  _flashmatch_nuslice_tree->Branch("charge_zw", &_charge_zw, "charge_zw/D");
  _flashmatch_nuslice_tree->Branch("charge_slope", &_charge_slope, "charge_slope/D");
  _flashmatch_nuslice_tree->Branch("charge_q", &_charge_q, "charge_q/D");
  _flashmatch_nuslice_tree->Branch("petoq", &_petoq, "petoq/D");
  _flashmatch_nuslice_tree->Branch("score", &_score, "score/D");
  _flashmatch_nuslice_tree->Branch("scr_y", &_scr_y, "scr_y/D");
  _flashmatch_nuslice_tree->Branch("scr_z", &_scr_z, "scr_z/D");
  _flashmatch_nuslice_tree->Branch("scr_rr", &_scr_rr, "scr_rr/D");
  _flashmatch_nuslice_tree->Branch("scr_ratio", &_scr_ratio, "scr_ratio/D");
  _flashmatch_nuslice_tree->Branch("scr_slope", &_scr_slope, "scr_slope/D");
  _flashmatch_nuslice_tree->Branch("scr_petoq", &_scr_petoq, "scr_petoq/D");
  _flashmatch_nuslice_tree->Branch("true_nus", &_true_nus, "true_nus/I");
  _flashmatch_nuslice_tree->Branch("mcT0", &_mcT0, "mcT0/D");
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
     !metricsInFile->Contains("ratio_h1") ||
     !metricsInFile->Contains("rr_fit_l") ||
     !metricsInFile->Contains("rr_fit_m") ||
     !metricsInFile->Contains("rr_fit_h") ||
     !metricsInFile->Contains("ratio_fit_l") ||
     !metricsInFile->Contains("ratio_fit_m") ||
     !metricsInFile->Contains("ratio_fit_h"))
  {
    mf::LogError("FlashPredict")
      << "The metrics file '" << fname << "' lacks at least one metric.";
    throw cet::exception("FlashPredict")
      << "The metrics file '" << fname << "'lacks at least one metric.";
  }
  //
  TH1 *temphisto = (TH1*)infile->Get("dy_h1");
  int bins = temphisto->GetNbinsX();
  if (bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for dy " << bins << " bins " << std::endl;
    bins = 1;
    fdYMeans.push_back(0);
    fdYSpreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= bins; ++ib) {
      fdYMeans.push_back(temphisto->GetBinContent(ib));
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in dy" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      fdYSpreads.push_back(tt);
    }
  }
  //
  temphisto = (TH1*)infile->Get("dz_h1");
  bins = temphisto->GetNbinsX();
  if (bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for dz " << bins << " bins " << std::endl;
    bins = 1;
    fdZMeans.push_back(0);
    fdZSpreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= bins; ++ib) {
      fdZMeans.push_back(temphisto->GetBinContent(ib));
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in dz" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      fdZSpreads.push_back(tt);
    }
  }
  //
  temphisto = (TH1*)infile->Get("rr_h1");
  bins = temphisto->GetNbinsX();
  if (bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for rr " << bins << " bins " << std::endl;
    bins = 1;
    fRRMeans.push_back(0);
    fRRSpreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= bins; ++ib) {
      double me = temphisto->GetBinContent(ib);
      fRRMeans.push_back(me);
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in rr" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      fRRSpreads.push_back(tt);
    }
    unsigned s = 0;
    for(auto& rrF : fRRFits){
      std::string nold = "rr_fit_" + kSuffixes[s];
      std::string nnew = "rrFit_" + kSuffixes[s];
      TF1* tempF1 = (TF1*)infile->Get(nold.c_str());
      auto params = tempF1->GetParameters();
      rrF.f = std::make_unique<TF1>(nnew.c_str(), fRR_TF1_fit.c_str(),
                                    0., fDriftDistance);
      rrF.f->SetParameters(params);
      rrF.min = rrF.f->GetMinimum(0., fDriftDistance, kEps);
      rrF.max = rrF.f->GetMaximum(0., fDriftDistance, kEps);
      s++;
    }
  }
  //
  temphisto = (TH1*)infile->Get("ratio_h1");
  bins = temphisto->GetNbinsX();
  if (bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for ratio " << bins << " bins " << std::endl;
    bins = 1;
    fRatioMeans.push_back(0);
    fRatioSpreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= bins; ++ib) {
      double me = temphisto->GetBinContent(ib);
      fRatioMeans.push_back(me);
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in ratio" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      fRatioSpreads.push_back(tt);
    }
    unsigned s = 0;
    for(auto& ratioF : fRatioFits){
      std::string nold = "ratio_fit_" + kSuffixes[s];
      std::string nnew = "ratioFit_" + kSuffixes[s];
      TF1* tempF1 = (TF1*)infile->Get(nold.c_str());
      auto params = tempF1->GetParameters();
      ratioF.f = std::make_unique<TF1>(nnew.c_str(), fRatio_TF1_fit.c_str(),
                                       0., fDriftDistance);
      ratioF.f->SetParameters(params);
      ratioF.min = ratioF.f->GetMinimum(0., fDriftDistance, kEps);
      ratioF.max = ratioF.f->GetMaximum(0., fDriftDistance, kEps);
      s++;
    }
  }
  // // TODO: enable the block below once some data is available
  // //
  // temphisto = (TH1*)infile->Get("slope_h1");
  // bins = temphisto->GetNbinsX();
  // if (bins <= 0 || fNoAvailableMetrics) {
  //   std::cout << " problem with input histos for slope " << bins << " bins " << std::endl;
  //   bins = 1;
  //   fSlopeMeans.push_back(0);
  //   fSlopeSpreads.push_back(0.001);
  // }
  // else {
  //   for (int ib = 1; ib <= bins; ++ib) {
  //     fSlopeMeans.push_back(temphisto->GetBinContent(ib));
  //     double tt = temphisto->GetBinError(ib);
  //     if (tt <= 0) {
  //       std::cout << "zero value for bin spread in slope" << std::endl;
  //       std::cout << "ib:\t" << ib << "\n";
  //       std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
  //       std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
  //       tt = 100.;
  //     }
  //     fSlopeSpreads.push_back(tt);
  //   }
  // }
  // //
  // temphisto = (TH1*)infile->Get("petoq_h1");
  // bins = temphisto->GetNbinsX();
  // if (bins <= 0 || fNoAvailableMetrics) {
  //   std::cout << " problem with input histos for petoq " << bins << " bins " << std::endl;
  //   bins = 1;
  //   fPEToQMeans.push_back(0);
  //   fPEToQSpreads.push_back(0.001);
  // }
  // else {
  //   for (int ib = 1; ib <= bins; ++ib) {
  //     fPEToQMeans.push_back(temphisto->GetBinContent(ib));
  //     double tt = temphisto->GetBinError(ib);
  //     if (tt <= 0) {
  //       std::cout << "zero value for bin spread in PEToQ" << std::endl;
  //       std::cout << "ib:\t" << ib << "\n";
  //       std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
  //       std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
  //       tt = 100.;
  //     }
  //     fPEToQSpreads.push_back(tt);
  //   }
  // }
  // // TODO: enable the block above

  TH2* tmp1_h2 = (TH2*)infile->Get("rr_h2");
  fRRH2 = (TH2D*)tmp1_h2->Clone("fRRH2");
  TH2* tmp2_h2 = (TH2*)infile->Get("ratio_h2");
  fRatioH2 = (TH2D*)tmp2_h2->Clone("fRatioH2");

  infile->Close();
  delete infile;
  mf::LogInfo("FlashPredict") << "Finish loading metrics";
}


double FlashPredict::cheatMCT0(const std::vector<art::Ptr<recob::Hit>>& hits,
                               const std::vector<art::Ptr<simb::MCParticle>>& mcParticles)
{
  auto clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  int pidMaxEnergy =
    TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, hits, true);
  for(auto& mcp: mcParticles){
    if(mcp->TrackId() == pidMaxEnergy){
      return mcp->Position().T();
    }
  }
  return -9999999.;
}


FlashPredict::ChargeMetrics FlashPredict::computeChargeMetrics(
  const ChargeDigest& chargeDigest) const
{
 // TODO: WIP. this function needs to be cleanud up before merging
  const unsigned pId = chargeDigest.pId;
  const int pfpPDGC = chargeDigest.pfpPDGC;
  const unsigned hitsInVolume = chargeDigest.hitsInVolume;
  const auto& qClusters = chargeDigest.qClusters;
  using Q_t = const flashmatch::QPoint_t&;
  double charge_q = std::accumulate(qClusters.begin(), qClusters.end(), 0.,
                                    [](double x, Q_t q) {return x+q.q;});
  if (charge_q <= 0.) return {};

  ChargeMetrics charge;
  charge.id = pId; charge.activity = hitsInVolume; charge.pdgc = pfpPDGC;
  charge.q = charge_q; charge.metric_ok = true;

  // double xmin = 9999; double ymin = 99999; double zmin = 999999;
  // double xmax = -9999; double ymax = -99999; double zmax = -999999;
  // double outlier_fraction = 0.025;
  // for (auto& qp : qClusters) {
  //   if(qp.x < xmin) xmin = qp.x;
  //   if(qp.y < ymin) ymin = qp.y;
  //   if(qp.z < zmin) zmin = qp.z;
  //   if(qp.x > xmax) xmax = qp.x;
  //   if(qp.y > ymax) ymax = qp.y;
  //   if(qp.z > zmax) zmax = qp.z;
  // }
  charge.x_gl =
    (std::accumulate(qClusters.begin(), qClusters.end(), 0.,
                     [](double x, Q_t q) {return x+q.q*q.x;}))/charge_q;
  charge.x = foldXGl(charge.x_gl);
  // std::vector<double> x_diff(qClusters.size());
  // std::transform(qClusters.begin(), qClusters.end(), x_diff.begin(),
  //                [charge](Q_t q) { return q.x - charge.x_gl; });
  // double x_sqsum = std::inner_product(x_diff.begin(), x_diff.end(), x_diff.begin(), 0.0);
  // charge.x_glw = std::sqrt(x_sqsum / (qClusters.size()-1));

  // {
  //   auto qCSorted = qClusters;
  //   std::sort(qCSorted.begin(), qCSorted.end(),
  //             [] (Q_t q1, Q_t q2){ return (q1.x < q2.x); });
  //   double qSum = 0.;
  //   for (size_t i=0; i<qCSorted.size(); ++i) {
  //     qSum += qCSorted[i].q;
  //     xmin = qCSorted[i].x;
  //     if (qSum >= charge_q*outlier_fraction) break;
  //   }
  //   qSum = 0.;
  //   for (size_t i=qCSorted.size()-1; i>=0; --i) {
  //     qSum += qCSorted[i].q;
  //     xmax = qCSorted[i].x;
  //     if (qSum >= charge_q*outlier_fraction) break;
  //   }
  // }
  // charge.x_glw = (xmax-xmin);
  // // charge.x_glw = (x_stdev<1.) ? 1. : x_stdev;

  charge.y =
    (std::accumulate(qClusters.begin(), qClusters.end(), 0.,
                     [](double y, Q_t q) {return y+q.q*q.y;}))/charge_q;
  // std::vector<double> y_diff(qClusters.size());
  // std::transform(qClusters.begin(), qClusters.end(), y_diff.begin(),
  //                [charge](Q_t q) { return q.y - charge.y; });
  // double y_sqsum = std::inner_product(y_diff.begin(), y_diff.end(), y_diff.begin(), 0.0);
  // charge.yw = std::sqrt(y_sqsum / (qClusters.size()-1));

  // {
  //   auto qCSorted = qClusters;
  //   std::sort(qCSorted.begin(), qCSorted.end(),
  //             [] (Q_t q1, Q_t q2){ return (q1.y < q2.y); });
  //   double qSum = 0.;
  //   for (size_t i=0; i<qCSorted.size(); ++i) {
  //     qSum += qCSorted[i].q;
  //     ymin = qCSorted[i].y;
  //     if (qSum >= charge_q*outlier_fraction) break;
  //   }
  //   qSum = 0.;
  //   for (size_t i=qCSorted.size()-1; i>=0; --i) {
  //     qSum += qCSorted[i].q;
  //     ymax = qCSorted[i].y;
  //     if (qSum >= charge_q*outlier_fraction) break;
  //   }
  // }
  // charge.yw = (ymax-ymin);
  // // charge.yw = (y_stdev<1.) ? 1. : y_stdev;

  charge.z =
    (std::accumulate(qClusters.begin(), qClusters.end(), 0.,
                     [](double z, Q_t q) {return z+q.q*q.z;}))/charge_q;
  // std::vector<double> z_diff(qClusters.size());
  // std::transform(qClusters.begin(), qClusters.end(), z_diff.begin(),
  //                [charge](Q_t q) { return q.z - charge.z; });
  // double z_sqsum = std::inner_product(z_diff.begin(), z_diff.end(), z_diff.begin(), 0.0);
  // charge.zw = std::sqrt(z_sqsum / (qClusters.size()-1));
  // {
  //   auto qCSorted = qClusters;
  //   std::sort(qCSorted.begin(), qCSorted.end(),
  //             [] (Q_t q1, Q_t q2){ return (q1.z < q2.z); });
  //   double qSum = 0.;
  //   for (size_t i=0; i<qCSorted.size(); ++i) {
  //     qSum += qCSorted[i].q;
  //     zmin = qCSorted[i].z;
  //     if (qSum >= charge_q*outlier_fraction) break;
  //   }
  //   qSum = 0.;
  //   for (size_t i=qCSorted.size()-1; i>=0; --i) {
  //     qSum += qCSorted[i].q;
  //     zmax = qCSorted[i].z;
  //     if (qSum >= charge_q*outlier_fraction) break;
  //   }
  // }
  // charge.zw = (zmax-zmin);
  // // charge.zw = (z_stdev<1.) ? 1. : z_stdev;

  // charge.slope = charge.zw / charge.yw;
  // charge.slope = charge.zw / (std::sqrt(charge.x_glw*charge.x_glw+charge.yw*charge.yw));
  //
  // double min_axis = std::min(charge.yw, charge.zw);
  // double max_axis = std::max(charge.yw, charge.zw);
  // double dir = (charge.zw>=charge.yw) ? 1. : -1.;
  // charge.slope = dir * std::sqrt(1 - ((min_axis*min_axis)/(max_axis*max_axis)));
  //
  //   charge.slope = charge.zw / csize;
  //

  std::unique_ptr<TH1F> qX = std::make_unique<TH1F>("qX", "", 80, -200., 200.);
  std::unique_ptr<TH1F> qY = std::make_unique<TH1F>("qY", "", 80, -200., 200.);
  std::unique_ptr<TH1F> qZ = std::make_unique<TH1F>("qZ", "", 100, 0., 500.);
  for (size_t i=0; i<qClusters.size(); ++i) {
    // double q2 = qClusters[i].q * qClusters[i].q;
    double q2 = qClusters[i].q;
    qX->Fill(qClusters[i].x, q2);
    qY->Fill(qClusters[i].y, q2);
    qZ->Fill(qClusters[i].z, q2);
  }
  charge.x_glw = qX->GetStdDev();
  charge.yw = qY->GetStdDev();
  charge.zw = qZ->GetStdDev();


  double csize = std::sqrt(charge.x_glw*charge.x_glw + charge.yw*charge.yw + charge.zw*charge.zw);
  // double c_xy_size = std::sqrt(charge.x_glw*charge.x_glw + charge.yw*charge.yw);
  // double csize = (charge.x_glw + charge.yw + charge.zw);
  // double csize = std::sqrt(charge.yw*charge.yw + charge.zw*charge.zw);
  double cx_fac = charge.x_glw / csize;
  double cy_fac = charge.yw / csize;
  double cz_fac = charge.zw / csize;
  // double cxy_fac = c_xy_size / csize;
  // double cy_fac = std::sqrt(charge.x_glw*charge.x_glw + charge.yw*charge.yw) / csize;
  // double cz_fac = std::sqrt(0.25*charge.x_glw*charge.x_glw + charge.zw*charge.zw) / csize;
  double min_axis = std::min(cy_fac, cz_fac);
  double max_axis = std::max(cy_fac, cz_fac);
  // double dir = ((1.2*cz_fac>=cx_fac) || (1.2*cz_fac>=cy_fac)) ? 1. : -1.;
  // double dir = ( (cz_fac>=cy_fac)) ? (max_axis-min_axis)/max_axis : (-max_axis+min_axis)/max_axis;
  double dir = ((cz_fac>=cy_fac)) ? 1. : -1.;
  // double dir = 1.;
  // if      ((cz_fac >= cy_fac) && (cz_fac >= cx_fac)) dir = 1.;
  // else if ((cz_fac >= cy_fac) && (cz_fac <  cx_fac)) dir = 0.5;
  // // else if ((cz_fac <  cy_fac) && (cz_fac <  cx_fac)) dir = -1.;
  // else dir = -1.;

  // charge.slope = dir * std::sqrt(1 - ((min_axis*min_axis)/(max_axis*max_axis)));
  charge.slope = dir * std::sqrt((max_axis*max_axis) - (min_axis*min_axis));
  // charge.slope = dir * (max_axis-min_axis)/max_axis;
  // charge.slope = dir * max_axis / min_axis;

  charge.x_glw = cx_fac;
  charge.yw = cy_fac;
  charge.zw = cz_fac;

  return charge;
}


FlashPredict::FlashMetrics FlashPredict::computeFlashMetrics(
  const SimpleFlash& simpleFlash) const
{
  // TODO: WIP. this function needs to be cleanud up before merging
  const OpHitIt opH_beg = simpleFlash.opH_beg;
  const OpHitIt opH_end = simpleFlash.opH_end;

  std::unique_ptr<TH1F> ophY = std::make_unique<TH1F>("ophY", "", fYBins, fYLow, fYHigh);
  std::unique_ptr<TH1F> ophZ = std::make_unique<TH1F>("ophZ", "", fZBins, fZLow, fZHigh);
  std::unique_ptr<TH1F> oph2Y = std::make_unique<TH1F>("oph2Y", "", fYBins, fYLow, fYHigh);
  std::unique_ptr<TH1F> oph2Z = std::make_unique<TH1F>("oph2Z", "", fZBins, fZLow, fZHigh);

  // std::list<double> opdetY;
  // std::list<double> opdetZ;

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

    ophY->Fill(opDetXYZ.Y(), ophPE);// TODO: ophPE or ophPE2
    ophZ->Fill(opDetXYZ.Z(), ophPE);// TODO: ophPE or ophPE2
    // if (!is_pmt_vis || !is_ara_vis){
      oph2Y->Fill(opDetXYZ.Y(), ophPE2);// TODO: ophPE or ophPE2
      oph2Z->Fill(opDetXYZ.Z(), ophPE2);// TODO: ophPE or ophPE2
    // }
    // opdetY.push_back(opDetXYZ.Y());
    // opdetZ.push_back(opDetXYZ.Z());

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

  // opdetY.sort();
  // opdetY.unique([](double l, double r) { return std::abs(l - r) < 0.00001;});
  // opdetZ.sort();
  // opdetZ.unique([](double l, double r) { return std::abs(l - r) < 0.00001;});
  // std::cout << "opdetY:\t";
  // for(auto Y: opdetY){
  //   std::cout << Y << "   ";
  // }
  // std::cout << std::endl;
  // std::cout << "opdetZ:\t";
  // for(auto Z: opdetZ){
  //   std::cout << Z << "   ";
  // }
  // std::cout << std::endl;

  if (sum_PE > 0.) {
    FlashMetrics flash;
    flash.metric_ok = true;
    flash.id = simpleFlash.flashId;
    flash.activity = simpleFlash.ophsInVolume;
    flash.pe    = sum_PE;
    flash.unpe  = sum_unPE;
    flash.y_skew = ophY->GetSkewness();
    flash.z_skew = ophZ->GetSkewness();
    flash.y_kurt = ophY->GetKurtosis();
    flash.z_kurt = ophZ->GetKurtosis();
    flash.xw = fractTimeWithFractionOfLight(simpleFlash, flash.pe, 0.8); // TODO: unharcode
    // TODO: low values of flash.xw <0.5 are indicative of good
    // mcT0-flash_time matching, so akin to matching to prompt light
    // Note that around the middle of the detector (~(+-100, 0, 250)
    // cm for SBND) the values of flash.xw are slightly larger, this
    // is natural and has to do with the way light disperses

    // TODO: decide which is the most useful
    // flash.yw = ophY->GetStdDev() / oph2Y->GetStdDev();
    // flash.zw = ophZ->GetStdDev() / oph2Z->GetStdDev();
    flash.yw = oph2Y->GetStdDev();
    flash.zw = oph2Z->GetStdDev();

    flash.ratio = fOpDetNormalizer * flash.unpe / flash.pe;
    if(fSBND){
      flash.ratio = (fOpDetNormalizer * flash.unpe) / (flash.pe - flash.unpe);
      if(fUseARAPUCAS) {
        // TODO: come up with a sensible way to mix PMT and ARAPUCA
        // flash ratios a potential better approach is to take the
        // average of the two separate ratios
        flash.unpe += sum_visARA_PE;
        flash.ratio = (fOpDetNormalizer * sum_unPE  + sum_visARA_PE ) / flash.pe;
      }
    }
    flash.yb  = sum_PE2Y / sum_PE2;
    flash.zb  = sum_PE2Z / sum_PE2;
    flash.rr = std::sqrt(
      std::abs(sum_PE2Y2 + sum_PE2Z2 + sum_PE2 * (flash.yb * flash.yb + flash.zb * flash.zb)
               - 2.0 * (flash.yb * sum_PE2Y + flash.zb * sum_PE2Z) ) / sum_PE2);

    //  TODO: corrections for short distance depositions that occur in
    //  front of problematic places.
    // SUPER EXPERIMENTAL.
    if(fSBND && false &&
       flash.xw < 0.236 &&
       // ((flash.yw<40) && (flash.zw<50)) &&  NOT GOOD:  widths grow when in front of uncoated
       ((flash.y_kurt>-1.) && (flash.z_kurt>0.))){// proxy for short distance
      // TODO these proxies using kurtosis have naturally low kurt values around the middles (y ~ 0, z ~ 250), that is the region we want to amend
      std::array<double, 3> hard_rr_y = {-135., 0., 135.}; // uncoated pmts
      std::array<double, 2> medm_rr_y = {-67.5, 67.5}; // gaps
      std::array<double, 3> hard_rr_z = {124.51, 254.7, 384.88}; // gaps
      std::array<double, 4> medm_rr_z = {57.87, 191.16, 318.24, 451.52}; // uncoated pmts

      std::array<double, 3> hard_ratio_y = hard_rr_y; // uncoated pmts
      // std::array<double, 3> medm_ratio_y = nahhh;
      std::array<double, 3> hard_ratio_z = hard_rr_z;
      std::array<double, 2> medm_ratio_z = {191.16, 318.24}; // uncoated pmts

      // double y_vicinity = 10000.; double z_vicinity = 10000.;
      // double hard_vicinity = 5.;
      // double hard_correction = 0.7;
      // double medm_vicinity = 10.;
      // double medm_correction = 0.82;

      for (double rr_y : hard_rr_y) {
        double vicinity = std::abs(flash.yb - rr_y);
        if(vicinity <  5.) {flash.rr = flash.rr * 0.70; break;}
        if(vicinity < 10.) {flash.rr = flash.rr * 0.80; break;}
        // if(vicinity < 15.) {flash.rr = flash.rr * 0.90; break;}
      }
      for (double rr_y : medm_rr_y) {
        double vicinity = std::abs(flash.yb - rr_y);
        if(vicinity <  5.) {flash.rr = flash.rr * 0.8; break;}
        if(vicinity < 10.) {flash.rr = flash.rr * 0.9; break;}
      }
      for (double rr_z : hard_rr_z) {
        double vicinity = std::abs(flash.zb - rr_z);
        if(vicinity <  5.) {flash.rr = flash.rr * 0.7; break;}
        if(vicinity < 10.) {flash.rr = flash.rr * 0.8; break;}
        // if(vicinity < 17.) {flash.rr = flash.rr * 0.9; break;}
      }
      for (double rr_z : medm_rr_z) {
        double vicinity = std::abs(flash.zb - rr_z);
        if(vicinity <  5.) {flash.rr = flash.rr * 0.8; break;}
        if(vicinity < 10.) {flash.rr = flash.rr * 0.9; break;}
      }

      for (double ratio_y : hard_ratio_y) {
        double vicinity = std::abs(flash.yb - ratio_y);
        if(vicinity <  5.) {flash.ratio = flash.ratio * 0.70; break;}
        if(vicinity < 10.) {flash.ratio = flash.ratio * 0.80; break;}
        // if(vicinity < 15.) {flash.ratio = flash.ratio * 0.90; break;}
      }
      // for (double ratio_y : medm_ratio_y) {
      //   double vicinity = std:abs(flash.yb - ratio_y);
      //   if(vicinity <  5.) {flash.ratio = flash.ratio * 0.85; break;}
      //   if(vicinity < 10.) {flash.ratio = flash.ratio * 0.95; break;}
      // }
      for (double ratio_z : hard_ratio_z) {
        double vicinity = std::abs(flash.zb - ratio_z);
        if(vicinity <  7.) {flash.ratio = flash.ratio * 0.82; break;}
        if(vicinity < 12.) {flash.ratio = flash.ratio * 0.89; break;}
        if(vicinity < 17.) {flash.ratio = flash.ratio * 0.96; break;}
      }
      for (double ratio_z : medm_ratio_z) {
        double vicinity = std::abs(flash.zb - ratio_z);
        if(vicinity <  7.) {flash.ratio = flash.ratio * 0.93; break;}
        if(vicinity < 12.) {flash.ratio = flash.ratio * 0.99; break;}
      }
    }
    //  TODO: corrections for short distance depositions that occur in
    //  front of problematic places


    std::tie(flash.h_x, flash.h_xerr, flash.h_xrr, flash.h_xratio) =
      hypoFlashX_H2(flash.rr, flash.ratio);

    // TODO: using _hypo_x make further corrections to _flash_time to
    // account for light transport time and/or rising edge
    // double flash_time = timeCorrections(simpleFlash.maxpeak_time, hypo_x);
    flash.time = simpleFlash.maxpeak_time;
    flash.x = peSumMax_wallX;
    flash.x_gl = flashXGl(flash.h_x, flash.x);

    // // TODO: maybe need to update the range before getting the skew to
    // // get a better representation? to only include bins that have
    // // hits
    // int firstNonEmptyYBin = 1;
    // for (size_t i=1; i<fYBins; ++i) {
    //   if (ophY->GetBinContent(i) > 0.) {
    //     firstNonEmptyYBin = i;
    //     break;
    //   }
    // }
    // int lastNonEmptyYBin = fYBins;
    // for (size_t i=fYBins; i>1; --i) {
    //   if (ophY->GetBinContent(i) > 0.) {
    //     lastNonEmptyYBin = i;
    //     break;
    //   }
    // }
    // int firstNonEmptyZBin = 1;
    // for (size_t i=1; i<fZBins; ++i) {
    //   if (ophZ->GetBinContent(i) > 0.) {
    //     firstNonEmptyZBin = i;
    //     break;
    //   }
    // }
    // int lastNonEmptyZBin = fZBins;
    // for (size_t i=fZBins; i>1; --i) {
    //   if (ophZ->GetBinContent(i) > 0.) {
    //     lastNonEmptyZBin = i;
    //     break;
    //   }
    // }
    // ophY->GetXaxis()->SetRange(firstNonEmptyYBin, lastNonEmptyYBin);
    // ophZ->GetXaxis()->SetRange(firstNonEmptyZBin, lastNonEmptyZBin);

    // ophY->ResetStats();
    // ophZ->ResetStats();

    // slope option 1
    // flash.slope = flash.zw / flash.yw;
    //
    // slope option 2
    // double min_axis = std::min(flash.yw, flash.zw);
    // double max_axis = std::max(flash.yw, flash.zw);
    // double dir = (flash.zw>=flash.yw) ? 1. : -1.;
    // flash.slope = dir * std::sqrt(1 - ((min_axis*min_axis)/(max_axis*max_axis)));
    //
    // another option
    double fsize = std::sqrt(flash.yw*flash.yw + flash.zw*flash.zw);
    //double fsize = (flash.yw + flash.zw);
    double fy_fac = flash.yw / fsize;
    double fz_fac = flash.zw / fsize;
    double min_axis = std::min(fy_fac, fz_fac);
    double max_axis = std::max(fy_fac, fz_fac);
    double dir = (fz_fac>=fy_fac) ? 1. : -1.;
    // double dir = ( (fz_fac>=fy_fac)) ? (max_axis-min_axis)/max_axis : (-max_axis+min_axis)/max_axis;
    // double dir = 1.;
    // flash.slope = dir * std::sqrt(1 - ((min_axis*min_axis)/(max_axis*max_axis)));
    flash.slope = dir * std::sqrt((max_axis*max_axis) - (min_axis*min_axis));
    // flash.slope = std::sqrt((max_axis*max_axis) - (min_axis*min_axis));
    // flash.slope = dir * (max_axis-min_axis)/max_axis;
    // flash.slope = fz_fac / fy_fac;

    // flash.yw = fy_fac;
    // flash.zw = fz_fac;

    flash.y = flash.yb - polinomialCorrection(flash.y_skew, flash.h_x,
                                              fPolCoefsY, fSkewLimitY);
    flash.z = flash.zb - polinomialCorrection(flash.z_skew, flash.h_x,
                                              fPolCoefsZ, fSkewLimitZ);

    return flash;
  }
  else {
    std::string channels;
    for(auto oph=opH_beg; oph!=opH_end; ++oph) channels += std::to_string(oph->OpChannel()) + ' ';
    mf::LogError("FlashPredict")
      << "Really odd that I landed here, this shouldn't had happen.\n"
      << "sum:          \t" << sum << "\n"
      << "sum_PE:       \t" << sum_PE << "\n"
      << "sum_unPE:     \t" << sum_unPE << "\n"
      << "opHits size:  \t" << std::distance(opH_beg, opH_end) << "\n"
      << "channels:     \t" << channels << std::endl;
    return {};
  }
}


FlashPredict::Score FlashPredict::computeScore(
  const ChargeMetrics& charge,
  const FlashMetrics& flash) const
{
  // TODO: WIP. this function needs to be cleanud up before merging
  Score score{0.};
  unsigned tcount = 0;
  double charge_x = (fCorrectDriftDistance) ?
    driftCorrection(charge.x, flash.time) : charge.x;
  int xbin = static_cast<int>(fXBins * (charge_x / fDriftDistance));
  if (charge_x < 0.) xbin = 0;
  else if (charge_x > fDriftDistance) xbin = fXBins - 1;

  score.y = scoreTerm(flash.y, charge.y, fdYMeans[xbin], fdYSpreads[xbin]);
  if(score.y > fTermThreshold) printMetrics("Y", charge, flash, score.y,
                                            mf::LogWarning("FlashPredict"));
  score.total += score.y;
  tcount++;
  score.z = scoreTerm(flash.z, charge.z, fdZMeans[xbin], fdZSpreads[xbin]);
  if(score.z > fTermThreshold) printMetrics("Z", charge, flash, score.z,
                                            mf::LogWarning("FlashPredict"));
  score.total += score.z;
  tcount++;
  score.rr = scoreTerm(flash.rr, fRRMeans[xbin], fRRSpreads[xbin]);
  if(score.rr > fTermThreshold) printMetrics("RR", charge, flash, score.rr,
                                             mf::LogWarning("FlashPredict"));
  score.total += score.rr;
  tcount++;
  if(fUseUncoatedPMT || fUseOppVolMetric) {
    score.ratio = scoreTerm(flash.ratio, fRatioMeans[xbin], fRatioSpreads[xbin]);
    if(fICARUS && !std::isnan(flash.h_x)){
      // TODO HACK to penalise matches with flash and charge on opposite volumes
      double charge_x_gl = (fCorrectDriftDistance) ?
        xGlCorrection(charge.x_gl, charge.x, flash.time) : charge.x_gl;
      double x_gl_diff = std::abs(flash.x_gl-charge_x_gl);
      double x_diff = std::abs(flash.h_x-charge_x);
      double cathode_tolerance = 30.;
      if(x_gl_diff > x_diff + cathode_tolerance) { // ok if close to the cathode
        double penalization = scoreTerm((flash.pe-flash.unpe)/flash.pe,
                                        fRatioMeans[xbin], fRatioSpreads[xbin]);
        score.ratio += penalization;
        mf::LogWarning("FlashPredict")
          << "HACK: Penalizing match with flash and charge in opposite volumes."
          << "\nflash.x_gl: " << flash.x_gl << " charge.x_gl: " << charge.x_gl
          << "\nX distance between them: " << x_gl_diff
          << "\nscore.ratio: " << score.ratio
          << "\nscore penalization: " << penalization;
      }
    }
    if(score.ratio > fTermThreshold) printMetrics("RATIO", charge, flash, score.ratio,
                                                  mf::LogWarning("FlashPredict"));
    score.total += score.ratio;
    tcount++;
  }
  // TODO: add these two new scores
  // score.slope = scoreTerm(flash.slope, charge.slope, fSlopeMeans[xbin], fSlopeSpreads[xbin]);
  score.slope = scoreTerm(flash.slope, charge.slope, -.18, .30);
  if(score.slope > fTermThreshold) printMetrics("SLOPE", charge, flash, score.slope,
                                                mf::LogWarning("FlashPredict"));
  // score.total += score.slope;
  // tcount++;
  // score.petoq = scoreTerm(std::log(flash.pe)/std::log(charge.q), fPEToQMeans[xbin], fPEToQSpreads[xbin]);
  score.petoq = (charge_x<70.) ? scoreTerm(PEToQ(flash.pe, charge.q), 0.71, 0.054) : scoreTerm(PEToQ(flash.pe, charge.q), 0.65, 0.046);
  if(score.petoq > fTermThreshold) printMetrics("LIGHT/CHARGE", charge, flash, score.petoq,
                                                mf::LogWarning("FlashPredict"));
  score.total += score.petoq;
  tcount++;
  mf::LogDebug("FlashPredict")
    << "score:       " << score.total << "using " << tcount << " terms\n"
    << "score.y:     " << score.y << "\t"
    << "score.z:     " << score.z << "\n"
    << "score.rr:    " << score.rr << "\t"
    << "score.ratio: " << score.ratio << "\n"
    << "score.slope: " << score.slope << "\t"
    << "score.petoq: " << score.petoq << "\n";
  return score;
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


std::tuple<double, double, double, double> FlashPredict::hypoFlashX_H2(
  double flash_rr, double flash_ratio) const
{
  auto[rr_hypoX, rr_hypoXRMS] = xEstimateAndRMS(flash_rr, fRRH2);
  auto[ratio_hypoX, ratio_hypoXRMS] = xEstimateAndRMS(flash_ratio, fRatioH2);

  double drr2 = rr_hypoXRMS*rr_hypoXRMS;
  double dratio2 = ratio_hypoXRMS*ratio_hypoXRMS;
  double rr_hypoXWgt = 1./drr2;
  double ratio_hypoXWgt = 1./dratio2;
  double sum_weights = rr_hypoXWgt + ratio_hypoXWgt;
  double hypo_x =
    (rr_hypoX*rr_hypoXWgt + ratio_hypoX*ratio_hypoXWgt) / sum_weights;
  if (std::abs(rr_hypoX - ratio_hypoX)<= 2.*std::sqrt(drr2+dratio2)){
    // consistent estimates, resulting error is smaller
    double hypo_x_err = std::sqrt(1./sum_weights);
    return {hypo_x, hypo_x_err, rr_hypoX, ratio_hypoX};
  }
  else {
    // inconsistent estimates, resulting error is larger
    double hypo_x_err = std::sqrt(drr2 + dratio2);
    return {hypo_x, hypo_x_err, rr_hypoX, ratio_hypoX};
  }
}


std::tuple<double, double> FlashPredict::xEstimateAndRMS(
  double metric_value, const TH2D* metric_h2) const
{
  int bin = metric_h2->GetYaxis()->FindBin(metric_value);
  int bins = metric_h2->GetNbinsY();
  int bin_buff = 0;
  // TODO: figure out a better method to make the estimates near the edges
  while(0 < bin-bin_buff || bin+bin_buff <= bins){
    int low_bin = (0 < bin-bin_buff) ? bin-bin_buff : 0;
    int high_bin = (bin+bin_buff <= bins) ? bin+bin_buff : -1;
    auto metric_px = std::make_unique<TH1D>(
      *(metric_h2->ProjectionX("metric_px", low_bin, high_bin)));
    if(metric_px->GetEntries() > kMinEntriesInProjection){
      double metric_hypoX = metric_px->GetRandom();
      // metric_hypoX = metric_px->GetMean(); // TODO: which one is more justified?
      double metric_rmsX = metric_px->GetRMS();
      if(metric_rmsX < fXBinWidth){//something went wrong // TODO: better test to see if things are OK
        mf::LogDebug("FlashPredict")
          << "metric_h2 projected on metric_value: "<< metric_value
          << ", bin: " << bin
          << ", bin_buff: " << bin_buff
          << "; has " << metric_px->GetEntries() << " entries."
          << "\nmetric_hypoX: " << metric_hypoX
          << ",  metric_rmsX: " << metric_rmsX;
        return {-1., 0.}; // no estimate
      }
      return {metric_hypoX, metric_rmsX};
    }
    bin_buff += 1;
  }
  return {-1., fDriftDistance}; // no estimate
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

  std::vector<art::Ptr<simb::MCParticle>> mcParticles;
  if(fStoreCheatMCT0){
    lar_pandora::LArPandoraHelper::CollectMCParticles(evt, "largeant",
                                                      mcParticles);
    // TODO: print particles identities, helpful for developing
  }

  std::unordered_map<size_t, size_t> pfpMap;
  for(size_t pId=0; pId<pfps_h->size(); pId++) {
    pfpMap[pfps_h->at(pId).Self()] = pId;
  }

  // TODO: this block is meant to only load the charge related objects
  // of the slices that are potentially neutrinos, this to improve
  // performance. Loading these objects is currently what that takes
  // the most time. Unfortunately, the objects are not very robust and
  // frequently there's out of bounds requests.
  //
  // std::vector<art::Ptr<recob::PFParticle>> particles_in_slices;
  // for(size_t pId=0; pId<pfps_h->size(); pId++) {
  //   if(!pfps_h->at(pId).IsPrimary()) continue;
  //   const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
  //   unsigned pfpPDGC = std::abs(pfp_ptr->PdgCode());
  //   if(fSelectNeutrino &&
  //      (pfpPDGC != 12) && (pfpPDGC != 14) && (pfpPDGC != 16)) continue;
  //   std::vector<art::Ptr<recob::PFParticle>> particles;
  //   addDaughters(pfpMap, pfp_ptr, pfps_h, particles);
  //   particles_in_slices.insert(particles_in_slices.end(),
  //                              particles.begin(), particles.end());
  // }
  // const art::FindManyP<recob::SpacePoint>
  //   pfp_spacepoints_assns(particles_in_slices, evt, fPandoraProducer);
  // std::vector<art::Ptr<recob::SpacePoint>> spacepoints_in_slices;
  // // for(auto& p: particles_in_slices) {
  // for(size_t p=0; p<particles_in_slices.size(); ++p) {
  //   auto& sps = pfp_spacepoints_assns.at(p);
  //   spacepoints_in_slices.insert(spacepoints_in_slices.end(), sps.begin(), sps.end());
  // }
  // const art::FindManyP<recob::Hit>
  //   spacepoint_hits_assns(spacepoints_in_slices, evt, fSpacePointProducer);
  // TODO: this block is meant to only load the charge related objects

  ChargeDigestMap chargeDigestMap;
  // Loop over pandora pfp particles
  for(size_t pId=0; pId<pfps_h->size(); pId++) {
    if(!pfps_h->at(pId).IsPrimary()) continue;
    const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
    unsigned pfpPDGC = std::abs(pfp_ptr->PdgCode());
    if(fSelectNeutrino &&
       (pfpPDGC != 12) && (pfpPDGC != 14) && (pfpPDGC != 16)){
      chargeDigestMap[-10.-pId] =
        ChargeDigest(pId, pfpPDGC, pfp_ptr, flashmatch::QCluster_t{}, 0, -9999.);
      continue;
    }
    std::set<unsigned> tpcWithHits;

    std::vector<art::Ptr<recob::PFParticle>> particles_in_slice;
    addDaughters(pfpMap, pfp_ptr, pfps_h, particles_in_slice);

    double sliceQ = 0.;
    std::vector<art::Ptr<recob::Hit>> hits_in_slice;
    flashmatch::QCluster_t particlesClusters;
    for(const auto& particle: particles_in_slice) {
      const auto particle_key = particle.key();
      const auto& particle_spacepoints = pfp_spacepoints_assns.at(particle_key);
      double particleQ = 0.;
      flashmatch::QCluster_t spsClusters;
      for(const auto& spacepoint : particle_spacepoints) {
        const auto spacepoint_key = spacepoint.key();
        const auto& hits = spacepoint_hits_assns.at(spacepoint_key);
        if(fStoreCheatMCT0){
          hits_in_slice.insert(hits_in_slice.end(),
                               hits.begin(), hits.end());
        }
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
    double mcT0 = (fStoreCheatMCT0) ? cheatMCT0(hits_in_slice, mcParticles)/1000.
      : -9999.;
    unsigned hitsInVolume = 0;
    bool in_right = false, in_left = false;
    for(unsigned t : tpcWithHits){
      if(t/fTPCPerDriftVolume == kRght) in_right = true;
      else if(t/fTPCPerDriftVolume == kLeft) in_left = true;
    }
    if(in_right && in_left) hitsInVolume = kActivityInBoth;
    else if(in_right && !in_left) hitsInVolume = kActivityInRght;
    else if(!in_right && in_left) hitsInVolume = kActivityInLeft;
    else {
      mf::LogError("FlashPredict")
        << "ERROR!!! tpcWithHits.size() " << tpcWithHits.size();
    }
    chargeDigestMap[sliceQ] =
      ChargeDigest(pId, pfpPDGC, pfp_ptr, particlesClusters, hitsInVolume, mcT0);
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


unsigned FlashPredict::trueNus(art::Event& evt) const
{
  art::Handle<std::vector<simb::MCTruth> > mctruthList_h;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel("generator", mctruthList_h))
    art::fill_ptr_vector(mclist, mctruthList_h);
  unsigned true_nus = 0;
  for(auto const& mc: mclist){
    if(mc->Origin() == simb::kBeamNeutrino) ++true_nus;
  }
  return true_nus;
}


inline
void FlashPredict::updateChargeMetrics(const ChargeMetrics& chargeMetrics)
{
  const auto& c = chargeMetrics;
  _charge_id = c.id; _charge_activity = c.activity; _charge_pdgc = c.pdgc;
  _charge_x = c.x; _charge_x_gl = c.x_gl; _charge_y = c.y; _charge_z = c.z;
  _charge_x_glw = c.x_glw; _charge_yw = c.yw; _charge_zw = c.zw;
  _charge_slope = c.slope; _charge_q = c.q;
}


inline
void FlashPredict::updateFlashMetrics(const FlashMetrics& flashMetrics)
{
  const auto& f = flashMetrics;
  _flash_id = f.id; _flash_activity = f.activity;
  _flash_x = f.x; _flash_yb = f.yb; _flash_zb = f.zb;
  _flash_x_gl = f.x_gl; _flash_y = f.y; _flash_z = f.z;
  _flash_xw = f.xw; _flash_yw = f.yw; _flash_zw = f.zw;
  _flash_rr = f.rr; _flash_ratio = f.ratio; _flash_slope = f.slope;
  _flash_pe = f.pe; _flash_unpe = f.unpe; _flash_time = f.time;
  _hypo_x = f.h_x; _hypo_x_err = f.h_xerr; _hypo_x_rr = f.h_xrr;
  _hypo_x_ratio = f.h_xratio;
  _y_skew = f.y_skew; _z_skew = f.z_skew; _y_kurt = f.y_kurt; _z_kurt = f.z_kurt;
}


inline
void FlashPredict::updateScore(const Score& score)
{
  _score = score.total, _scr_y = score.y, _scr_z = score.z,
    _scr_rr = score.rr, _scr_ratio = score.ratio,
    _scr_slope = score.slope, _scr_petoq = score.petoq;
}


inline
double FlashPredict::scoreTerm(const double m, const double n,
                               const double mean, const double spread) const
{
  return std::abs((m - n) - mean) / spread;
}


inline
double FlashPredict::scoreTerm(const double m,
                               const double mean, const double spread) const
{
  return std::abs(m - mean) / spread;
}


inline
double FlashPredict::PEToQ(const double pe, const double q) const
{
  return std::log(pe)/std::log(q);
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


void FlashPredict::copyOpHitsInFlashFindingWindow(
  std::vector<recob::OpHit>& opHits,
  const art::Handle<std::vector<recob::OpHit>>& ophits_h) const
{
  double s = fFlashFindingTimeStart;
  double e = fFlashFindingTimeEnd;
  double m = fMinOpHPE;
  // copy ophits that are inside the time window and with PEs
  auto opHitInWindow =
    [s, e, m, this](const recob::OpHit& oph)-> bool
      {return ((opHitTime(oph) > s) &&
               (opHitTime(oph) < e) &&
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
    "opHitsTimeHist", "ophittime", fTimeBins, fFlashFindingTimeStart, fFlashFindingTimeEnd);
  opHitsTimeHist->SetOption("HIST");
  opHitsTimeHist->SetDirectory(0);//turn off ROOT's object ownership
  std::unique_ptr<TH1D> opHitsTimeHistRght = std::make_unique<TH1D>(
    "opHitsTimeHistRght", "ophittimer", fTimeBins, fFlashFindingTimeStart, fFlashFindingTimeEnd);
  opHitsTimeHistRght->SetOption("HIST");
  opHitsTimeHistRght->SetDirectory(0);//turn off ROOT's object ownership
  std::unique_ptr<TH1D> opHitsTimeHistLeft = std::make_unique<TH1D>(
    "opHitsTimeHistLeft", "ophittimel", fTimeBins, fFlashFindingTimeStart, fFlashFindingTimeEnd);
  opHitsTimeHistLeft->SetOption("HIST");
  opHitsTimeHistLeft->SetDirectory(0);//turn off ROOT's object ownership
  if(!createOpHitsTimeHist(
       opHits, opHitsRght, opHitsLeft,
       opHitsTimeHist, opHitsTimeHistRght, opHitsTimeHistLeft)) return {};

  bool oph_in_rght = false, oph_in_left = false;
  std::vector<FlashPredict::SimpleFlash> simpleFlashes;
  if(opHitsRght.size() > 0 && opHitsTimeHistRght->GetEntries() > 0){
    oph_in_rght = true;
    findSimpleFlashes(simpleFlashes, opHitsRght,
                      kActivityInRght, opHitsTimeHistRght);
  }
  if(opHitsLeft.size() > 0 && opHitsTimeHistLeft->GetEntries() > 0){
    oph_in_left = true;
    findSimpleFlashes(simpleFlashes, opHitsLeft,
                      kActivityInLeft, opHitsTimeHistLeft);
  }
  if(oph_in_rght && oph_in_left ){
    findSimpleFlashes(simpleFlashes, opHits, kActivityInBoth, opHitsTimeHist);
  }
  else if(!oph_in_rght && !oph_in_left){
    return {};
  }
  return simpleFlashes;
}


//ICARUS overload
std::vector<FlashPredict::SimpleFlash> FlashPredict::makeSimpleFlashes(
  std::vector<recob::OpHit>& opHits) const
{
  std::unique_ptr<TH1D> opHitsTimeHist = std::make_unique<TH1D>(
    "opHitsTimeHist", "ophittime", fTimeBins, fFlashFindingTimeStart, fFlashFindingTimeEnd);
  opHitsTimeHist->SetOption("HIST");
  opHitsTimeHist->SetDirectory(0);//turn off ROOT's object ownership
  unsigned ophsInVolume = createOpHitsTimeHist(opHits, opHitsTimeHist);
  if(ophsInVolume == 0) return {};

  std::vector<FlashPredict::SimpleFlash> simpleFlashes;
  if(!findSimpleFlashes(simpleFlashes, opHits, ophsInVolume, opHitsTimeHist))
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
    opHitsTimeHist->Fill(opHitTime(oph), oph.PE());
    if(sbndPDinTPC(ch) == kRght){
      opHitsRght.emplace_back(oph);
      opHitsTimeHistRght->Fill(opHitTime(oph), oph.PE());
    }
    else{// sbndPDinTPC(ch) == kLeft
      opHitsLeft.emplace_back(oph);
      opHitsTimeHistLeft->Fill(opHitTime(oph), oph.PE());
    }
  }
  if(opHitsTimeHist->GetEntries() <= 0 ||
     opHitsTimeHist->Integral() <= 0. ||
     opHitsTimeHist->Integral() <= fMinFlashPE) return false;
  if(opHitsTimeHistRght->GetEntries() <= 0 ||
     opHitsTimeHistRght->Integral() <= 0. ||
     opHitsTimeHistRght->Integral() <= fMinFlashPE)
    opHitsTimeHistRght->Reset();
  if(opHitsTimeHistLeft->GetEntries() <= 0 ||
     opHitsTimeHistLeft->Integral() <= 0. ||
     opHitsTimeHistLeft->Integral() <= fMinFlashPE)
    opHitsTimeHistLeft->Reset();
  return true;
}


//ICARUS overload
unsigned FlashPredict::createOpHitsTimeHist(
  const std::vector<recob::OpHit>& opHits,
  std::unique_ptr<TH1D>& opHitsTimeHist) const
{
  bool in_right = false, in_left = false;
  for(auto const& oph : opHits) {
    auto ch = oph.OpChannel();
    auto opDetXYZ = fGeometry->OpDetGeoFromOpChannel(ch).GetCenter();
    if(!fGeoCryo->ContainsPosition(opDetXYZ)) continue;
    opHitsTimeHist->Fill(opHitTime(oph), oph.PE());
    unsigned t = icarusPDinTPC(ch);
    if(t/fTPCPerDriftVolume == kRght) in_right = true;
    else if(t/fTPCPerDriftVolume == kLeft) in_left = true;
  }
  if(opHitsTimeHist->GetEntries() <= 0 ||
     opHitsTimeHist->Integral() <= 0. ||
     opHitsTimeHist->Integral() <= fMinFlashPE) return 0;
  if(in_right && in_left) return kActivityInBoth;
  else if(in_right && !in_left) return kActivityInRght;
  else if(!in_right && in_left) return kActivityInLeft;
  else return 0;
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
    if (flashId == 0 || flashId == 1) { // First flashes have to be within the beam spill
      int beam_start_bin = opHitsTimeHist->FindBin(fBeamSpillTimeStart);
      int beam_end_bin = opHitsTimeHist->FindBin(fBeamSpillTimeEnd);
      opHitsTimeHist->GetXaxis()->SetRange(beam_start_bin, beam_end_bin);
      int ibin_beam = opHitsTimeHist->GetMaximumBin();
      maxpeak_time = opHitsTimeHist->GetBinCenter(ibin_beam);
      opHitsTimeHist->GetXaxis()->SetRange(0, 0); // reset range
    }
    double lowedge  = maxpeak_time + fFlashStart;
    double highedge = maxpeak_time + fFlashEnd;
    mf::LogDebug("FlashPredict")
      << "light window " << lowedge << " " << highedge << std::endl;
    int lowedge_bin = opHitsTimeHist->FindBin(lowedge);
    int highedge_bin = opHitsTimeHist->FindBin(highedge);
    // check if flash has enough PEs, skip if is the first flash
    if (opHitsTimeHist->Integral(lowedge_bin, highedge_bin) <= fMinFlashPE ||
        opHitsTimeHist->Integral(lowedge_bin, highedge_bin) <= 0. ){
      if(flashId == 0 || flashId == 1) continue;
      break;
    }
    // clear this peak to enforce non-overlapping flashes
    for(int i=lowedge_bin; i<highedge_bin; ++i){
      opHitsTimeHist->SetBinContent(i, 0.);
    }
    auto peakInsideEdges =
      [lowedge, highedge, this](const recob::OpHit& oph)-> bool
        { return ((lowedge <= opHitTime(oph)) && (opHitTime(oph) <= highedge)); };
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
unsigned FlashPredict::sbndPDinTPC(const int pdChannel) const
{
  auto p = fGeometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  p.SetX(p.X()/2.);//OpDets are outside the TPCs
  return (fGeometry->PositionToTPCID(p)).TPC;
}


// TODO: find better, less hacky solution
unsigned FlashPredict::icarusPDinTPC(const int pdChannel) const
{
  auto p = fGeometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  if(fCryostat == 0) p.SetX((p.X() + 222.)/2. - 222.);//OpDets are outside the TPCs
  if(fCryostat == 1) p.SetX((p.X() - 222.)/2. + 222.);//OpDets are outside the TPCs
  return (fGeometry->PositionToTPCID(p)).TPC;
}


// bool FlashPredict::isSBNDPDRelevant(const int pdChannel,
//                                     const std::set<unsigned>& tpcWithHits) const
// {
//   // if there's hits on all TPCs all channels are relevant
//   if(tpcWithHits.size() == fNTPC) return true;
//   unsigned pdTPC = sbndPDinTPC(pdChannel);
//   for(auto itpc: tpcWithHits) if(itpc == pdTPC) return true;
//   return false;
// }


double FlashPredict::opHitTime(const recob::OpHit& oph) const
{
  // TODO: maybe add some guards to prevent breaking backwards
  // compatibility
  if      (fUseOpHitRiseTime)  return oph.RiseTime();
  else if (fUseOpHitPeakTime)  return oph.PeakTime();
  else if (fUseOpHitStartTime) return oph.StartTime();
  return oph.PeakTime();
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


double FlashPredict::fractTimeWithFractionOfLight(
  const SimpleFlash& simpleFlash, double sum_pe, double fraction_pe) const
{
  std::vector<recob::OpHit> timeSortedOpH(simpleFlash.opH_beg, simpleFlash.opH_end);
  std::sort(timeSortedOpH.begin(), timeSortedOpH.end(),
            [this] (const recob::OpHit& oph1, const recob::OpHit& oph2)
              { return (opHitTime(oph1) < opHitTime(oph2)); });
  const double threshold_fraction = fraction_pe * sum_pe;
  const double flash_start_time = opHitTime(timeSortedOpH[0]);
  const double flash_end_time = opHitTime(timeSortedOpH[timeSortedOpH.size()-1]);
  double threshold_time = 0.;
  double running_sum = 0.;
  for (auto& oph : timeSortedOpH) {
    running_sum += oph.PE();
    if (running_sum >= threshold_fraction){
      threshold_time = opHitTime(oph);
      break;
    }
  }
  return (threshold_time-flash_start_time) / (flash_end_time-flash_start_time);
}


inline
double FlashPredict::polinomialCorrection(
  const double skew, const double hypo_x,
  const std::vector<double>& polCoefs,
  const double skew_limit) const
{
  if (std::abs(skew)>skew_limit || std::isnan(hypo_x)) return 0.;
  // TODO: maybe some other condition to prevent corrections?
  double correction = 0.;
  double exponent = 1.;
  for (double coef : polCoefs) {
    correction += coef * exponent;
    exponent *= hypo_x;
  }
  return correction * skew;
}


std::list<double> FlashPredict::wiresXGl() const
{
  std::list<double> wiresX_gl;
  for (size_t t = 0; t < fNTPC; t++) {
    const geo::TPCGeo& tpcg = fGeoCryo->TPC(t);
    wiresX_gl.push_back(tpcg.LastPlane().GetCenter().X());
  }
  wiresX_gl.unique([](double l, double r) { return std::abs(l - r) < 0.00001;});
  return wiresX_gl;
}


double FlashPredict::driftDistance() const
{
  auto wit = fWiresX_gl.begin();
  auto wite = fWiresX_gl.begin();
  wite++;
  return std::abs(*wite - *wit)/2.;
}


inline
double FlashPredict::driftCorrection(const double c_xb, const double f_time) const
{
  // TODO: best would be to use some SCE service function to
  // adequately take into account SCE
  return c_xb - (0.16  * f_time); // TODO: use drift speed from larsoft
}


inline
double FlashPredict::xGlCorrection(const double c_x_glb, const double c_xb,
                                   const double f_time) const
{
  double c_x = driftCorrection(c_xb, f_time);
  if (std::abs(foldXGl(driftCorrection(c_x_glb, f_time))-c_x)<0.001)
    return driftCorrection(c_x_glb, f_time);
  else
    return driftCorrection(c_x_glb, -1.*f_time);
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


double FlashPredict::foldXGl(const double x_gl) const
{
  auto wit = fWiresX_gl.begin();
  for(size_t i=0; i<fWiresX_gl.size()-1; i++){
    double wxl = *wit;
    wit++;
    double wxh = *wit;
    if(wxl < x_gl && x_gl<= wxh){
      double mid = (wxl + wxh)/2.;
      return (x_gl<=mid) ? std::abs(x_gl-wxl) : std::abs(x_gl-wxh);
    }
  }
  mf::LogWarning("FlashPredict")
    << "Global charge center X position is out of bounds: " << x_gl;
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
  if(bk.noslice) {
    m << "\tNo Slice:         \t -" << bk.noslice
      << ", scored as: " << kNoSlcInEvt << "\n";
  }
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
  if(bk.no_nu_candidate) {
    m << "\tNo Nu Candidate:\t -" << bk.no_nu_candidate
      << ", scored as: " << kNotANuScr << "\n";
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
    - bk.nopfpneutrino - bk.noslice
    - bk.nonvalidophit - bk.nullophittime;

  // account for the reasons that a particle might lack
  bk.pfp_bookkeeping = bk.pfp_to_score
    - bk.no_oph_hits - bk.no_charge
    - bk.no_flash_pe - bk.no_nu_candidate;

  if(1-std::abs(bk.events_processed-bk.job_bookkeeping) == 0 ||
     1-std::abs(bk.scored_pfp-bk.pfp_bookkeeping) == 0){
    printBookKeeping(mf::LogWarning("FlashPredict"));
  }
}


template <typename Stream>
void FlashPredict::printMetrics(const std::string metric,
                                const ChargeMetrics& charge,
                                const FlashMetrics& flash,
                                const double term,
                                Stream&& out) const
{
  int xbin = static_cast<int>(fXBins * (_charge_x / fDriftDistance));
  out
    << "Big term " << metric << ":\t" << term << "\n"
    << std::left << std::setw(12) << std::setfill(' ')
    << "xbin:        \t" << xbin << "\n"
    << "_slices:    \t" << std::setw(8) << _slices   << "\n"
    << "_true_nus:  \t" << std::setw(8) << _true_nus << "\n"
    << "_mcT0:      \t" << std::setw(8) << _mcT0     << "\n"
    << "_petoq:     \t" << std::setw(8) << _petoq    << "\n"
    << "charge metrics:\n" << charge.dumpMetrics()   << "\n"
    << "flash metrics:\n"  << flash.dumpMetrics()    << "\n";
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
