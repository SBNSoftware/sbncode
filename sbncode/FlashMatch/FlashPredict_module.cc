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
  , fOnlyPrimaries(p.get<bool>("OnlyPrimaries", true))
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
  , fMinSliceQ(p.get<double>("MinSliceQ", 0.0))
  , fQScale(p.get<double>("QScale", 1.0))
  , fMinOpHPE(p.get<double>("MinOpHPE", 0.0))
  , fMinFlashPE(p.get<double>("MinFlashPE", 0.0))
  , fPEScale(p.get<double>("PEScale", 1.0))
  , fCryostat(p.get<int>("Cryostat", 0)) //set =0 ot =1 for ICARUS to match reco chain selection
  , fNBins(p.get<int>("n_bins"))
  , fDriftDistance(p.get<double>("DriftDistance"))// rounded up for binning
  , fOpDetNormalizer(p.get<unsigned>("OpDetNormalizer", 1))
  , fTermThreshold(p.get<double>("ThresholdTerm", 30.))
{
  produces< std::vector<sbn::SimpleFlashMatch> >();
  produces< art::Assns <recob::PFParticle, sbn::SimpleFlashMatch> >();
  // fFlashProducer         = p.get<art::InputTag>("FlashProducer");

  fPDMapAlgPtr = art::make_tool<opdet::PDMapAlg>(p.get<fhicl::ParameterSet>("PDMapAlg"));
  fGeoCryo = std::make_unique<geo::CryostatGeo>(geometry->Cryostat(fCryostat));
  fNTPC = geometry->NTPC();
  for (size_t t = 0; t < fNTPC; t++) {
    const geo::TPCGeo& tpcg = fGeoCryo->TPC(t);
    fWiresX_gl.push_back(tpcg.LastPlane().GetCenter().X());
  }
  fWiresX_gl.unique([](double l, double r) { return std::abs(l - r) < 0.00001;});

  fDetector = geometry->DetectorName();
  if(fDetector.find("sbnd") != std::string::npos) {
    fDetector = "SBND";
    fSBND = true;
    fICARUS = false;
    fTPCPerDriftVolume = 1;
    fDriftVolumes = fNTPC/fTPCPerDriftVolume;
  }
  else if (fDetector.find("icarus") != std::string::npos) {
    fDetector = "ICARUS";
    fSBND = false;
    fICARUS = true;
    fTPCPerDriftVolume = 2;
    fDriftVolumes = fNTPC/fTPCPerDriftVolume;
  }
  else {
      throw cet::exception("FlashPredict")
        << "Detector: " << fDetector
        << ", not supported. Stopping.\n";
  }

  // TODO no point on having fCryostat as parameter, user whatever comes from geometry
  if (fSBND && fCryostat == 1) {
    throw cet::exception("FlashPredict")
      << "SBND has only one cryostat. \n"
      << "Check Detector and Cryostat parameter." << std::endl;
  }
  else if (fICARUS && fCryostat > 1) {
    throw cet::exception("FlashPredict")
      << "ICARUS has only two cryostats. \n"
      << "Check Detector and Cryostat parameter." << std::endl;
  }

  fOpHitsTimeHist = std::make_unique<TH1D>("fOpHitsTimeHist", "ophittime", fTimeBins,
                                           fBeamWindowStart, fBeamWindowEnd); // in us
  fOpHitsTimeHist->SetOption("HIST");
  fOpHitsTimeHist->SetDirectory(0);//turn off ROOT's object ownership

  if(!fOnlyPrimaries){
    mf::LogWarning("FlashPredict")
      << "The fcl option OnlyPrimaries is useless, sorry. I'll remove it soon.";
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

  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfps_h, evt,
                                                          fPandoraProducer);

  auto const& spacepoint_h =
    evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointProducer);
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h,
                                                   evt, fSpacePointProducer);

  // grab tracks associated with PFParticles
  // auto const& track_h = evt.getValidHandle<std::vector<recob::Track> >(fTrackProducer);
  // art::FindManyP<recob::Track> pfp_track_assn_v(track_h, evt, fTrackProducer);

  // grab calorimetry info for tracks
  // auto const& calo_h =
  //   evt.getValidHandle<std::vector<anab::Calorimetry> >(fCaloProducer);
  // art::FindManyP<anab::Calorimetry>  track_calo_assn_v(calo_h,
  //                                                      evt, fCaloProducer);

  // load OpHits previously created
  art::Handle<std::vector<recob::OpHit>> ophit_h;
  evt.getByLabel(fOpHitProducer, ophit_h);
  if(!ophit_h.isValid()) {
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

  std::vector<recob::OpHit> opHits(ophit_h->size());
  copyOpHitsInBeamWindow(opHits, ophit_h);

  if(fUseARAPUCAS && !fOpHitARAProducer.empty()){
    art::Handle<std::vector<recob::OpHit>> ophitara_h;
    evt.getByLabel(fOpHitARAProducer, ophitara_h);
    if(!ophitara_h.isValid()) {
      mf::LogWarning("FlashPredict")
        << "Non valid ophits from ARAPUCAS"
        << "\nfUseARAPUCAS: " << std::boolalpha << fUseARAPUCAS
        << "\nfOpHitARAProducer: " << fOpHitARAProducer;
    }
    else{
      std::vector<recob::OpHit> opHitsARA(ophitara_h->size());
      copyOpHitsInBeamWindow(opHitsARA, ophitara_h);
      opHits.insert(opHits.end(),
                    opHitsARA.begin(), opHitsARA.end());
    }
  }

  fOpHitsTimeHist->Reset();
  fPeakCounter = 0;
  if(!getOpHitsInFlash(opHits)){
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
  // TODO: for SBND optically isolated volumes create 2 more opHits:
  // opHitsOnlyOnTheLeft, opHitsOnlyOnTheRight... or some better name/object

  std::set<unsigned> tpcWithOpH;
  if(fSBND) {// no point for ICARUS
    for(auto oph=fOpH_beg; oph!=fOpH_end; ++oph){
      tpcWithOpH.insert(sbndPDinTPC(oph->OpChannel()));
      if(tpcWithOpH.size() == fNTPC) break;
    }
  }

  std::map<size_t, size_t> pfpMap;
  for (size_t pId=0; pId<pfps_h->size(); pId++) {
    pfpMap[pfps_h->at(pId).Self()] = pId;
  }

  std::map<double, ChargeDigest, std::greater<double>> chargeDigestMap;
  // Loop over pandora pfp particles
  for(size_t pId=0; pId<pfps_h->size(); pId++) {
    if(!pfps_h->at(pId).IsPrimary()) continue;
    const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pId);
    unsigned pfpPDGC = std::abs(pfp_ptr->PdgCode());
    if(fSelectNeutrino &&
        (pfpPDGC != 12) && (pfpPDGC != 14) && (pfpPDGC != 16) ) continue;
    bk.pfp_to_score++;
    flashmatch::QCluster_t qClusters;
    std::set<unsigned> tpcWithHits;

    {//TODO: pack this into a function
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    AddDaughters(pfpMap, pfp_ptr, pfps_h, pfp_ptr_v);

    double chargeToNPhotons = lar_pandora::LArPandoraHelper::IsTrack(pfp_ptr) ?
      fChargeToNPhotonsTrack : fChargeToNPhotonsShower;
    double totalCharge = 0;
    //  loop over all mothers and daughters, fill qCluster
    for (auto& pfp_md: pfp_ptr_v) {
      auto key = pfp_md.key();
      /*
        if ( fUseCalo && lar_pandora::LArPandoraHelper::IsTrack(pfp_ptr)) {
        // grab tracks associated with pfp particle
        auto const& track_ptr_v = pfp_track_assn_v.at(key);
        for (size_t tr=0; tr < track_ptr_v.size(); tr++) {
        auto mytrack = track_ptr_v[tr];
        auto const& trackkey = mytrack.key();
        // grab calo objects associated with tracks
        const std::vector< art::Ptr<anab::Calorimetry> > calo_ptr_v = track_calo_assn_v.at( trackkey );
        for (size_t ca=0;  ca <  calo_ptr_v.size(); ca++) {
        auto mycalo = calo_ptr_v.at( ca );
        int npts = mycalo->dEdx().size();
        for (int ip=0;ip<npts;++ip) {
        Point_t pxyz=mycalo->fXYZ[ip];
        double ds = mycalo->fTrkPitch[ip];
        double dQdx = mycalo->fdQdx[ip];
        double dEdx = mycalo->fdEdx[ip];
        double alpha = dQdx/dEdx;
        double charge = (1-alpha)*dEdx*ds;
        // hardcode for now for SBND
        double xpos = 0.0;
        if (pxyz[0]<0) xpos = fabs(pxyz[0]+200.0);
        else xpos = fabs(200.0-pxyz[0]);
        qClusterInTPC[tpcindex].emplace_back(xpos, position[1], position[2], charge);
        }
        }
        }
        }
        else { // this is a shower
      */

      auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      for (auto& SP : spacepoint_ptr_v) {
        auto const& spkey = SP.key();
        const auto& hit_ptr_v = spacepoint_hit_assn_v.at(spkey);
        for (auto& hit : hit_ptr_v) {
          // TODO: Add hits from induction wires too.
          // Only use hits from the collection plane
          geo::WireID wid = hit->WireID();
          if (geometry->SignalType(wid) != geo::kCollection) continue;
          const auto& pos(SP->XYZ());
          auto itpc = wid.TPC;
          tpcWithHits.insert(itpc);
          const auto charge(hit->Integral());
          if (charge < fMinHitQ) continue;
          totalCharge += charge;
          qClusters.emplace_back(pos[0], pos[1], pos[2], charge * chargeToNPhotons);
        } // for all hits associated to this spacepoint
      } // for all spacepoints
      //      }  // if track or shower
    } // for all pfp pointers
    if (totalCharge < fMinSliceQ) continue;
    chargeDigestMap[totalCharge] = ChargeDigest(pId, pfpPDGC, pfp_ptr,
                                                qClusters, tpcWithHits);
    }//TODO: pack this into a function
  } // over all PFParticles

  for(auto& chargeDigest : chargeDigestMap) {
    //const size_t pId = chargeDigest.second.pId;
    const int pfpPDGC = chargeDigest.second.pfpPDGC;
    const auto& pfp_ptr = chargeDigest.second.pfp_ptr;
    const auto& qClusters = chargeDigest.second.qClusters;
    const auto& tpcWithHits = chargeDigest.second.tpcWithHits;

    if(fSBND){// because SBND has an opaque cathode
      std::set<unsigned> tpcWithHitsOpH;
      std::set_intersection(tpcWithHits.begin(), tpcWithHits.end(),
                            tpcWithOpH.begin(), tpcWithOpH.end(),
                            std::inserter(tpcWithHitsOpH, tpcWithHitsOpH.begin()));
      if (tpcWithHitsOpH.size() == 0) {
        mf::LogWarning("FlashPredict") << "No OpHits where there's charge. Skipping...";
        bk.no_oph_hits++;
        mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
        sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                             Flash(kNoScrPE), Score(kQNoOpHScr)));
        util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
        continue;
      }
    }

    if(!computeChargeMetrics(qClusters)){
      mf::LogWarning("FlashPredict") << "Clusters with No Charge. Skipping...";
      bk.no_charge++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                           Flash(kNoScrPE), Score(kNoChrgScr)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      continue;
    }

    if(!computeFlashMetrics(tpcWithHits)){
      printMetrics("ERROR", pfpPDGC, tpcWithHits, 0, mf::LogError("FlashPredict"));
      bk.no_flash_pe++;
      mf::LogDebug("FlashPredict") << "Creating sFM and PFP-sFM association";
      sFM_v->push_back(sFM(kNoScr, kNoScrTime, Charge(kNoScrQ),
                           Flash(kNoScrPE), Score(k0VUVPEScr)));
      util::CreateAssn(*this, evt, *sFM_v, pfp_ptr, *pfp_sFM_assn_v);
      continue;
    }

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
    for(auto& rrF : rrFits){
      std::string nold = "rr_fit_" + suffixes[s];
      std::string nnew = "rrFit_" + suffixes[s];
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
    for(auto& peF : peFits){
      std::string nold = "pe_fit_" + suffixes[s];
      std::string nnew = "peFit_" + suffixes[s];
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
  _charge_q = 0.;
  for (auto& qp : qClusters) {
    xave += fQScale * qp.q * qp.x;
    yave += fQScale * qp.q * qp.y;
    zave += fQScale * qp.q * qp.z;
    norm += fQScale * qp.q;
    _charge_q += qp.q;
  }
  if (norm > 0) {
    _charge_x_gl = xave / norm;
    _charge_x = driftDistance(_charge_x_gl);
    _charge_y = yave / norm;
    _charge_z = zave / norm;
    return true;
  }
  return false;
}


bool FlashPredict::computeFlashMetrics(const std::set<unsigned>& tpcWithHits)
{
  // NOTE: _flash_x holds the X coordinate of the opdet that registered
  //       the most PEs in the flash.
  // TODO: change _flash_x to hold the X coordinate of the wall that
  //       registered the most PEs overall
  auto compareOpHits = [] (const recob::OpHit& oph1, const recob::OpHit& oph2)->bool
    { return oph1.PE() < oph2.PE(); };
  auto opHMax = std::max_element(fOpH_beg, fOpH_end, compareOpHits);
  _flash_x = geometry->OpDetGeoFromOpChannel(opHMax->OpChannel()).GetCenter().X();

  double sum = 0.;
  double sum_PE = 0.;
  double sum_PE2 = 0.;
  double sum_unPE = 0.;
  double sum_visARA_PE = 0.;
  double sum_PE2Y  = 0.; double sum_PE2Z  = 0.;
  double sum_PE2Y2 = 0.; double sum_PE2Z2 = 0.;

  for(auto oph=fOpH_beg; oph!=fOpH_end; ++oph){
    if(fSBND && !isSBNDPDRelevant(oph->OpChannel(), tpcWithHits)) continue;
    auto opDet = geometry->OpDetGeoFromOpChannel(oph->OpChannel());
    auto opDetXYZ = opDet.GetCenter();

    std::string op_type = "pmt"; // the label ICARUS has
    if(fSBND){
      op_type = fPDMapAlgPtr->pdType(oph->OpChannel());
      if(!fUseUncoatedPMT && op_type == "pmt_uncoated") continue;
    }

    double ophPE2 = oph->PE() * oph->PE();
    sum       += 1.0;
    sum_PE    += oph->PE();
    sum_PE2   += ophPE2;
    sum_PE2Y  += ophPE2 * opDetXYZ.Y();
    sum_PE2Z  += ophPE2 * opDetXYZ.Z();
    sum_PE2Y2 += ophPE2 * opDetXYZ.Y() * opDetXYZ.Y();
    sum_PE2Z2 += ophPE2 * opDetXYZ.Z() * opDetXYZ.Z();

    if(fICARUS){
      // TODO: change this if block
      if(fUseOppVolMetric){
        unsigned pdVolume = icarusPDinTPC(oph->OpChannel())/fTPCPerDriftVolume;
        geo::Point_t q(_charge_x_gl, _charge_y, _charge_z);
        unsigned qVolume = driftVolume(_charge_x_gl);
        if(qVolume < fDriftVolumes){
          if (pdVolume != qVolume){
            sum_unPE += oph->PE();
          }
        }
      }
      // // TODO: for this block instead, needs testing!
      // if(fUseOppVolMetric && _flash_x != opDetXYZ.X()){
      //   sum_unPE += oph->PE();
      // }
    }
    else if(fSBND){
      if(fUseUncoatedPMT && op_type == "pmt_uncoated") {
        sum_unPE += oph->PE();
      }
      else if(fUseARAPUCAS &&
              (op_type == "arapuca_vis" || op_type == "xarapuca_vis")) {
        sum_visARA_PE += oph->PE();
      }
    }
  } // for opHits

  if (sum_PE > 0) {
    _flash_pe    = sum_PE   * fPEScale;
    _flash_unpe  = sum_unPE * fPEScale;
    _flash_ratio = fOpDetNormalizer * _flash_unpe / _flash_pe;
    if(fUseARAPUCAS) {
      _flash_unpe  += sum_visARA_PE * fPEScale;
      _flash_ratio = (fOpDetNormalizer * sum_unPE  + sum_visARA_PE )* fPEScale / _flash_pe;
    }
    _flash_y  = sum_PE2Y / sum_PE2;
    _flash_z  = sum_PE2Z / sum_PE2;
    _flash_r = std::sqrt(
      std::abs(sum_PE2Y2 + sum_PE2Z2 + sum_PE2 * (_flash_y * _flash_y + _flash_z * _flash_z)
       - 2.0 * (_flash_y * sum_PE2Y + _flash_z * sum_PE2Z) ) / sum_PE2);
    // _hypo_x = hypoFlashX_splines();
    _hypo_x = hypoFlashX_fits();
    _flash_x_gl = flashXGl(_hypo_x, _flash_x);
    return true;
  }
  else {
    std::string channels;
    for(auto oph=fOpH_beg; oph!=fOpH_end; ++oph) channels += std::to_string(oph->OpChannel()) + ' ';
    std::string tpcs;
    for(auto itpc: tpcWithHits) tpcs += std::to_string(itpc) + ' ';
    mf::LogError("FlashPredict")
      << "Really odd that I landed here, this shouldn't had happen.\n"
      << "sum:          \t" << sum << "\n"
      << "sum_PE:       \t" << sum_PE << "\n"
      << "sum_unPE:     \t" << sum_unPE << "\n"
      << "tpcWithHits:  \t" << tpcs << "\n"
      << "opHits size:  \t" << std::distance(fOpH_beg, fOpH_end) << "\n"
      << "channels:     \t" << channels << std::endl;
    _flash_y = 0;
    _flash_z = 0;
    _flash_r = 0;
    _flash_pe = 0;
    _flash_unpe = 0;
    _flash_ratio = 0;
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
  if ((fSBND && fUseUncoatedPMT) ||
      (fICARUS && fUseOppVolMetric)) {
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


double FlashPredict::hypoFlashX_fits()
{
  std::vector<double> rrXs;
  double rr_hypoX, rr_hypoXWgt;
  for(const auto& rrF : rrFits){
    if(rrF.min < _flash_r && _flash_r < rrF.max){
      try{
        rrXs.push_back(rrF.f->GetX(_flash_r, 0., fDriftDistance, kEps));
      }catch (...) {
        mf::LogWarning("FlashPredict")
          << "Couldn't find root with rrFits.\n"
          << "min/_flash_ratio/max:"
          << rrF.min << "/" << _flash_r << "/" << rrF.max;
      }
    }
  }
  if(rrXs.size() > 1){
    rr_hypoX = (rrXs[0] + rrXs[1])/2.;
    rr_hypoXWgt =  1./std::abs(rrXs[0] - rrXs[1]);
  }
  else if(rrXs.size() == 0){
    rr_hypoX = 0.;
    rr_hypoXWgt = 0.;
  }
  else{//(rrXs.size() == 1)
    if(_flash_r < rrFits[2].min){
      rr_hypoX =  rrXs[0]/2.;
      rr_hypoXWgt = 1./std::abs(rrXs[0]);
    }
    else{
      rr_hypoX = (rrXs[0] + fDriftDistance)/2.;
      rr_hypoXWgt = 1./std::abs(rrXs[0] - fDriftDistance);
    }
  }

  _hypo_x_rr = rr_hypoX;
  if(!fUseUncoatedPMT && !fUseOppVolMetric)
    return rr_hypoX;

  std::vector<double> peXs;
  double pe_hypoX, pe_hypoXWgt;
  for(const auto& peF : peFits){
    if(peF.min < _flash_ratio && _flash_ratio < peF.max){
      try{
        peXs.push_back(peF.f->GetX(_flash_ratio, 0., fDriftDistance, kEps));
      }catch (...) {
        mf::LogWarning("FlashPredict")
          << "Couldn't find root with peFits.\n"
          << "min/_flash_ratio/max:"
          << peF.min << "/" << _flash_ratio << "/" << peF.max;
      }
    }
  }
  if(peXs.size() > 1){
    pe_hypoX = (peXs[0] + peXs[1])/2.;
    pe_hypoXWgt =  1./std::abs(peXs[0] - peXs[1]);
  }
  else if(peXs.size() == 0){
    pe_hypoX = 0.;
    pe_hypoXWgt = 0.;
  }
  else{// peXs.size() == 1
    if(_flash_ratio < peFits[2].min){
      pe_hypoX =   peXs[0]/2.;
      pe_hypoXWgt =  1./std::abs(peXs[0]);
    }
    else{
      pe_hypoX = (peXs[0] + fDriftDistance)/2.;
      pe_hypoXWgt = 1./std::abs(peXs[0] - fDriftDistance);
    }
  }
  _hypo_x_ratio = pe_hypoX;

  return (rr_hypoX*rr_hypoXWgt + pe_hypoX*pe_hypoXWgt) / (rr_hypoXWgt + pe_hypoXWgt);
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


void FlashPredict::AddDaughters(
  const std::map<size_t, size_t>& pfpMap,
  const art::Ptr<recob::PFParticle>& pfp_ptr,
  const art::ValidHandle<std::vector<recob::PFParticle>>& pfps_h,
  std::vector<art::Ptr<recob::PFParticle>>& pfp_v) const
{
  auto daughters = pfp_ptr->Daughters();
  pfp_v.push_back(pfp_ptr);
  for(auto const& daughterid : daughters) {
    if (pfpMap.find(daughterid) == pfpMap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error" << std::endl;
      continue;
    }
    const art::Ptr<recob::PFParticle> pfp_ptr(pfps_h, pfpMap.at(daughterid) );
    AddDaughters(pfpMap, pfp_ptr, pfps_h, pfp_v);
  } // for all daughters
  return;
} // void FlashPredict::AddDaughters


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
  const art::Handle<std::vector<recob::OpHit>>& ophit_h) const
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
  auto it = std::copy_if(ophit_h->begin(), ophit_h->end(), opHits.begin(),
                         opHitInWindow);
  opHits.resize(std::distance(opHits.begin(), it));
}


bool FlashPredict::getOpHitsInFlash(std::vector<recob::OpHit>& opHits)
{
  if(fOpHitsTimeHist->GetEntries() == 0){// only create it once per event
    if(!createOpHitsTimeHist(opHits)) return false;
  }
  if(!findMaxPeak(opHits)) return false;
  return true;
}


bool FlashPredict::createOpHitsTimeHist(
  const std::vector<recob::OpHit>& opHits) const
{
  for(auto const& oph : opHits) {
    auto ch = oph.OpChannel();
    auto opDetXYZ = geometry->OpDetGeoFromOpChannel(ch).GetCenter();
    if (fICARUS &&
        !fGeoCryo->ContainsPosition(opDetXYZ)) continue;
    if(fSBND && !fUseUncoatedPMT &&
       !fPDMapAlgPtr->isPDType(oph.OpChannel(), "pmt_coated")) continue;
    fOpHitsTimeHist->Fill(oph.PeakTime(), fPEScale * oph.PE());
  }
  if (fOpHitsTimeHist->GetEntries() <= 0 ||
      fOpHitsTimeHist->Integral() < fMinFlashPE) return false;
  return true;
}


bool FlashPredict::findMaxPeak(std::vector<recob::OpHit>& opHits)
{
  int ibin = fOpHitsTimeHist->GetMaximumBin();
  double maxpeak_time = fOpHitsTimeHist->GetBinCenter(ibin);
  _flash_time = maxpeak_time; // in us
  double lowedge  = _flash_time + fFlashStart;
  double highedge = _flash_time + fFlashEnd;
  mf::LogDebug("FlashPredict") << "light window " << lowedge << " " << highedge << std::endl;
  int lowedge_bin = fOpHitsTimeHist->FindBin(lowedge);
  int highedge_bin = fOpHitsTimeHist->FindBin(highedge);
  if (fOpHitsTimeHist->Integral(lowedge_bin, highedge_bin) < fMinFlashPE){
    return false;
  }
  for(int i=lowedge_bin; i<highedge_bin; ++i){// clear this peak
    fOpHitsTimeHist->SetBinContent(i, 0.);
  }
  auto peakInsideEdges =
    [&lowedge, &highedge](const recob::OpHit& oph)-> bool
      { return ((lowedge <= oph.PeakTime()) && (oph.PeakTime() <= highedge)); };
  // partition container to move the hits of the peak
  // the iterators point to the boundaries of the partition
  fOpH_beg = (fPeakCounter > 0) ? fOpH_end : opHits.begin();
  fOpH_end = std::partition(fOpH_beg, opHits.end(),
                            peakInsideEdges);
  fPeakCounter++;
  return true;
}


bool FlashPredict::isPDInCryo(const int pdChannel) const
{
  if(fSBND) return true;
  else if(fICARUS){
    // BUG: I believe this function is not working, every now and then
    // I get ophits from the other cryo
    auto& p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
    // if the channel is in the Cryostat is relevant
    return fGeoCryo->ContainsPosition(p);
  }
  return false;
}


// TODO: find better, less hacky solution
unsigned FlashPredict::icarusPDinTPC(const int pdChannel) const
{
  auto p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  if(fCryostat == 0) p.SetX((p.X() + 222.)/2. - 222.);//OpDets are outside the TPCs
  if(fCryostat == 1) p.SetX((p.X() - 222.)/2. + 222.);//OpDets are outside the TPCs
  return (geometry->PositionToTPCID(p)).TPC;
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
  auto p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  p.SetX(p.X()/2.);//OpDets are outside the TPCs
  return (geometry->PositionToTPCID(p)).TPC;
}


double FlashPredict::flashXGl(const double hypo_x,
                              const double flash_x) const
{
  double min = 10000;
  unsigned wId = 0;
  auto wIt = fWiresX_gl.begin();
  auto w = wIt;
  for(size_t i=0; i<fWiresX_gl.size(); i++){
    if(((*w<0) == (flash_x<0))&& std::abs(*w - flash_x) < min) {
      wIt = w;
      wId = i;
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
//   auto p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
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
