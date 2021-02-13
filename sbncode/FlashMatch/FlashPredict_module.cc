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
  , fLightWindowStart(p.get<double>("LightWindowStart")) // in us w.r.t. flash time
  , fLightWindowEnd(p.get<double>("LightWindowEnd"))  // in us w.r.t flash time
  , fTimeBins(unsigned(1/fTickPeriod * (fBeamWindowEnd - fBeamWindowStart)))
  , fSelectNeutrino(p.get<bool>("SelectNeutrino", true))
  , fUseUncoatedPMT(p.get<bool>("UseUncoatedPMT", false))
  , fUseOppVolMetric(p.get<bool>("UseOppVolMetric", false))
  , fUseARAPUCAS(p.get<bool>("UseARAPUCAS", false))
  // , fUseCalo(p.get<bool>("UseCalo", false))
  , fInputFilename(p.get<std::string>("InputFileName")) // root file with score metrics
  , fNoAvailableMetrics(p.get<bool>("NoAvailableMetrics", false))
  , fMakeTree(p.get<bool>("MakeTree", false))
  , fMinFlashPE(p.get<double>("MinFlashPE", 0.0))
  , fMinOpHPE(p.get<double>("MinOpHPE", 0.0))
  , fPEscale(p.get<double>("PEscale", 1.0))
  , fChargeToNPhotonsShower(p.get<double>("ChargeToNPhotonsShower", 1.0))  // ~40000/1600
  , fChargeToNPhotonsTrack(p.get<double>("ChargeToNPhotonsTrack", 1.0))  // ~40000/1600
  , fCryostat(p.get<int>("Cryostat", 0)) //set =0 ot =1 for ICARUS to match reco chain selection
  , fNBins(p.get<int>("n_bins"))
  , fDriftDistance(p.get<double>("DriftDistance"))// rounded up for binning
  , fVUVToVIS(p.get<unsigned>("VUVToVIS", 4))
  , fTermThreshold(p.get<double>("ThresholdTerm", 30.))
{
  produces< std::vector< anab::T0 > >();
  produces< art::Assns <recob::PFParticle, anab::T0> >();
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

  if (fMakeTree) initTree();

  loadMetrics();

  consumes<std::vector<recob::PFParticle>>(fPandoraProducer);
  consumes<art::Assns<recob::SpacePoint, recob::PFParticle>>(fPandoraProducer);
  consumes<std::vector<recob::SpacePoint>>(fSpacePointProducer);
  consumes<art::Assns<recob::Hit, recob::SpacePoint>>(fSpacePointProducer);
  consumes<std::vector<recob::OpHit>>(fOpHitProducer);
} // FlashPredict::FlashPredict(fhicl::ParameterSet const& p)


void FlashPredict::produce(art::Event & e)
{

  std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< art::Assns <recob::PFParticle, anab::T0> > pfp_t0_assn_v( new art::Assns<recob::PFParticle, anab::T0>  );

  // reset TTree variables
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  _flash_time    = -9999.;
  _flash_pe      = -9999.;
  _flash_unpe    = -9999.;
  _flash_r       = -9999.;
  _flash_ratio   = -9999.;
  _hypo_x        = -9999.;
  _score         = -9999.;
  bk.events++;

  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPandoraProducer);

  if (fSelectNeutrino &&
      !pfpNeutrinoOnEvent(pfp_h)) {
    mf::LogInfo("FlashPredict")
      << "No pfp neutrino on event. Skipping...";
    bk.nopfpneutrino++;
    updateBookKeeping();
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPandoraProducer);

  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointProducer);
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointProducer);

  // grab tracks associated with PFParticles
  // auto const& track_h = e.getValidHandle<std::vector<recob::Track> >(fTrackProducer);
  // art::FindManyP<recob::Track> pfp_track_assn_v(track_h, e, fTrackProducer);

  // grab calorimetry info for tracks
  // auto const& calo_h = e.getValidHandle<std::vector<anab::Calorimetry> >(fCaloProducer);
  // art::FindManyP<anab::Calorimetry>  track_calo_assn_v(calo_h, e, fCaloProducer);

  // load OpHits previously created
  art::Handle<std::vector<recob::OpHit>> ophit_h;
  e.getByLabel(fOpHitProducer, ophit_h);
  if(!ophit_h.isValid()) {
    mf::LogError("FlashPredict")
      << "No optical hits from producer module "
      << fOpHitProducer;
    bk.nonvalidophit++;
    updateBookKeeping();
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  std::vector<recob::OpHit> opHits(ophit_h->size());
  copyOpHitsInBeamWindow(opHits, ophit_h);

  if(fUseARAPUCAS && !fOpHitARAProducer.empty()){
    art::Handle<std::vector<recob::OpHit>> ophitara_h;
    e.getByLabel(fOpHitARAProducer, ophitara_h);
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
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  std::set<unsigned> tpcWithOpH;
  if(fSBND) {// no point for ICARUS
    for(auto oph=fOpH_beg; oph!=fOpH_end; ++oph){
      tpcWithOpH.insert(sbndPDinTPC(oph->OpChannel()));
      if(tpcWithOpH.size() == fNTPC) break;
    }
  }

  _pfpmap.clear();
  for (size_t p=0; p<pfp_h->size(); p++) _pfpmap[pfp_h->at(p).Self()] = p;

  // Loop over pandora pfp particles
  for (size_t p=0; p<pfp_h->size(); p++) {
    auto const& pfp = pfp_h->at(p);
    unsigned pfpPDGC = std::abs(pfp_h->at(p).PdgCode());
    if (!pfp.IsPrimary()) continue;
    if (fSelectNeutrino &&
        (pfpPDGC != 12) &&
        (pfpPDGC != 14) &&
        (pfpPDGC != 16) ) continue;
    bk.pfp_to_score++;
    flashmatch::QCluster_t qClusters;
    std::set<unsigned> tpcWithHits;

    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);

    double chargeToNPhotons = lar_pandora::LArPandoraHelper::IsTrack(pfp_ptr) ? fChargeToNPhotonsTrack : fChargeToNPhotonsShower;
    //  loop over all mothers and daughters, fill qCluster
    for (auto& pfp_md: pfp_ptr_v) {
      auto key = pfp_md.key();
      pfp_v.push_back(*pfp_md);

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
      std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
      for (auto& SP : spacepoint_ptr_v) {
        auto const& spkey = SP.key();
        const auto& this_hit_ptr_v = spacepoint_hit_assn_v.at(spkey);
        for (auto& hit : this_hit_ptr_v) {
          // TODO: Add hits from induction wires too.
          // Only use hits from the collection plane
          geo::WireID wid = hit->WireID();
          if (geometry->SignalType(wid) != geo::kCollection) continue;
          const auto& pos(SP->XYZ());
          auto itpc = wid.TPC;
          tpcWithHits.insert(itpc);
          const auto charge(hit->Integral());
          qClusters.emplace_back(pos[0], pos[1], pos[2], charge * chargeToNPhotons);
        } // for all hits associated to this spacepoint
      } // for all spacepoints
      //      }  // if track or shower
    } // for all pfp pointers

    if(fSBND){// because SBND has an opaque cathode
      std::set<unsigned> tpcWithHitsOpH;
      std::set_intersection(tpcWithHits.begin(), tpcWithHits.end(),
                            tpcWithOpH.begin(), tpcWithOpH.end(),
                            std::inserter(tpcWithHitsOpH, tpcWithHitsOpH.begin()));
      if (tpcWithHitsOpH.size() == 0) {
        mf::LogWarning("FlashPredict") << "No OpHits where there's charge. Skipping...";
        bk.no_oph_hits++;
        mf::LogDebug("FlashPredict") << "Creating T0 and PFP-T0 association";
        T0_v->push_back(anab::T0(-9999., 0, p, -9999., kQNoOpHScr));
        util::CreateAssn(*this, e, *T0_v, pfp_ptr, *pfp_t0_assn_v);
        continue;
      }
    }

    if(!computeChargeMetrics(qClusters)){
      mf::LogWarning("FlashPredict") << "Clusters with No Charge. Skipping...";
      bk.no_charge++;
      mf::LogDebug("FlashPredict") << "Creating T0 and PFP-T0 association";
      T0_v->push_back(anab::T0(-9999., 0, p, -9999., kNoChrgScr));
      util::CreateAssn(*this, e, *T0_v, pfp_ptr, *pfp_t0_assn_v);
      continue;
    }

    if(!computeFlashMetrics(tpcWithHits)){
      printMetrics("ERROR", pfp.PdgCode(), tpcWithHits, 0, mf::LogError("FlashPredict"));
      bk.no_flash_pe++;
      mf::LogDebug("FlashPredict") << "Creating T0 and PFP-T0 association";
      T0_v->push_back(anab::T0(-9999., 0, p, -9999., k0VUVPEScr));
      util::CreateAssn(*this, e, *T0_v, pfp_ptr, *pfp_t0_assn_v);
      continue;
    }

    _hypo_x = hypoFlashX_splines();

    if(computeScore(tpcWithHits, pfp.PdgCode())){
      if (fMakeTree) {_flashmatch_nuslice_tree->Fill();}
      bk.scored_pfp++;
      mf::LogDebug("FlashPredict") << "Creating T0 and PFP-T0 association";
      T0_v->push_back(anab::T0(_flash_time, icountPE, p, _hypo_x, _score));
      util::CreateAssn(*this, e, *T0_v, pfp_ptr, *pfp_t0_assn_v);
    }
  } // over all PFparticles

  bk.events_processed++;
  updateBookKeeping();

  e.put(std::move(T0_v));
  e.put(std::move(pfp_t0_assn_v));

}// end of producer module


void FlashPredict::initTree(void)
{
  art::ServiceHandle<art::TFileService> tfs;
  _flashmatch_nuslice_tree = tfs->make<TTree>("nuslicetree", "nu FlashPredict tree");
  _flashmatch_nuslice_tree->Branch("evt", &_evt, "evt/I");
  _flashmatch_nuslice_tree->Branch("run", &_run, "run/I");
  _flashmatch_nuslice_tree->Branch("sub", &_sub, "sub/I");
  _flashmatch_nuslice_tree->Branch("flash_time", &_flash_time, "flash_time/D");
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
  _flashmatch_nuslice_tree->Branch("score", &_score, "score/D");
  _flashmatch_nuslice_tree->Branch("scr_y", &_scr_y, "scr_y/D");
  _flashmatch_nuslice_tree->Branch("scr_z", &_scr_z, "scr_z/D");
  _flashmatch_nuslice_tree->Branch("scr_rr", &_scr_rr, "scr_rr/D");
  _flashmatch_nuslice_tree->Branch("scr_ratio", &_scr_ratio, "scr_ratio/D");
}


void FlashPredict::loadMetrics()
{
  // TODO: Set a better way to run with no metrics

  // TODO: fill histos with less repetition and range for loops
  // read histograms and fill vectors for match score calculation
  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fInputFilename, fname);
  mf::LogInfo("FlashPredict") << "Opening file with metrics: " << fname;
  TFile *infile = new TFile(fname.c_str(), "READ");
  if(!infile->IsOpen()) {
    throw cet::exception("FlashPredict")
      << "Could not find the light-charge match root file '"
      << fname << "'!\n";
  }
  auto metricsInFile = infile->GetListOfKeys();
  if(!metricsInFile->Contains("dy_h1") ||
     !metricsInFile->Contains("dz_h1") ||
     !metricsInFile->Contains("rr_h1") ||
     !metricsInFile->Contains("pe_h1"))
  {
    throw cet::exception("FlashPredict")
      << "The metrics file lacks at least one metric.";
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
    std::vector<double> x, yH, yL;
    for (int ib = 1; ib <= bins; ++ib) {
      double me = temphisto->GetBinContent(ib);
      rr_means.push_back(me);
      x.push_back(temphisto->GetBinCenter(ib));
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in rr" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      rr_spreads.push_back(tt);
      yL.push_back(me - tt);
      yH.push_back(me + tt);
    }
    rr_InvSpl = TSpline3("rr_InvSpl", rr_means.data(), x.data(), rr_means.size());
    rr_l_InvSpl = TSpline3("rr_l_InvSpl", yL.data(), x.data(), yL.size());
    rr_h_InvSpl = TSpline3("rr_h_InvSpl", yH.data(), x.data(),  yH.size());
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
    std::vector<double> x, yH, yL;
    for (int ib = 1; ib <= bins; ++ib) {
      double me = temphisto->GetBinContent(ib);
      pe_means.push_back(me);
      x.push_back(temphisto->GetBinCenter(ib));
      double tt = temphisto->GetBinError(ib);
      if (tt <= 0) {
        std::cout << "zero value for bin spread in pe" << std::endl;
        std::cout << "ib:\t" << ib << "\n";
        std::cout << "temphisto->GetBinContent(ib):\t" << temphisto->GetBinContent(ib) << "\n";
        std::cout << "temphisto->GetBinError(ib):\t" << temphisto->GetBinError(ib) << "\n";
        tt = 100.;
      }
      pe_spreads.push_back(tt);
      yL.push_back(me - tt);
      yH.push_back(me + tt);
    }
    pe_InvSpl = TSpline3("pe_InvSpl", pe_means.data(), x.data(), pe_means.size());
    pe_l_InvSpl = TSpline3("pe_l_InvSpl", yL.data(), x.data(), yL.size());
    pe_h_InvSpl = TSpline3("pe_h_InvSpl", yH.data(), x.data(), yH.size());
  }

  infile->Close();
  delete infile;
  mf::LogInfo("FlashPredict") << "Finish loading metrics";
}


bool FlashPredict::computeChargeMetrics(flashmatch::QCluster_t& qClusters)
{
  double xave = 0.; double yave = 0.;
  double zave = 0.; double norm = 0.;
  double scale = 0.001;
  _charge_q = 0.;
  for (auto& qp : qClusters) {
    xave += scale * qp.q * qp.x;
    yave += scale * qp.q * qp.y;
    zave += scale * qp.q * qp.z;
    norm += scale * qp.q;
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


bool FlashPredict::computeFlashMetrics(std::set<unsigned>& tpcWithHits)
{
  double sum = 0.;
  double sum_PE = 0.;
  double sum_PE2 = 0.;
  double sum_unPE = 0.;
  double sum_visARA_PE = 0.;
  double sum_PE2Y  = 0.; double sum_PE2Z  = 0.;
  double sum_PE2Y2 = 0.; double sum_PE2Z2 = 0.;

  double maxPE = -9999.;
  for(auto oph=fOpH_beg; oph!=fOpH_end; ++oph){
    if(!isPDRelevant(oph->OpChannel(), tpcWithHits)) continue;
    auto opDet = geometry->OpDetGeoFromOpChannel(oph->OpChannel());
    auto opDetXYZ = opDet.GetCenter();

    std::string op_type = "pmt"; // the label ICARUS has
    if(fSBND) op_type = fPDMapAlgPtr->pdType(oph->OpChannel());

    // _flash_x is the X coord of the opdet where most PE are deposited
    if(oph->PE()>maxPE){
      _flash_x = opDetXYZ.X();
      maxPE = oph->PE();
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
    _flash_pe    = sum_PE   * fPEscale;
    _flash_unpe  = sum_unPE * fPEscale;
    _flash_ratio = fVUVToVIS * _flash_unpe / _flash_pe;
    if(fUseARAPUCAS) {
      _flash_unpe  += sum_visARA_PE * fPEscale;
      _flash_ratio = (fVUVToVIS * sum_unPE  + sum_visARA_PE )* fPEscale / _flash_pe;
    }
    _flash_y  = sum_PE2Y / sum_PE2;
    _flash_z  = sum_PE2Z / sum_PE2;
    _flash_r = std::sqrt(
      std::abs(sum_PE2Y2 + sum_PE2Z2 + sum_PE2 * (_flash_y * _flash_y + _flash_z * _flash_z)
       - 2.0 * (_flash_y * sum_PE2Y + _flash_z * sum_PE2Z) ) / sum_PE2);
    icountPE = std::round(_flash_pe);
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


bool FlashPredict::computeScore(std::set<unsigned>& tpcWithHits, int pdgc)
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


double FlashPredict::hypoFlashX_splines()
{
  double hX, mX, lX;
  double rr_X = rr_InvSpl.Eval(_flash_r);
  double rr_hX = rr_h_InvSpl.Eval(_flash_r);
  double rr_lX = rr_l_InvSpl.Eval(_flash_r);
  if(0.<rr_X && rr_X<fDriftDistance) mX = rr_X;
  else mX = (_flash_r<rr_means[0]) ? 0. : fDriftDistance;
  hX = (0.<rr_hX && rr_hX<fDriftDistance) ? rr_hX : mX;
  lX = (0.<rr_lX && rr_lX<fDriftDistance) ? rr_lX : mX;
  double rr_hypoX =  (lX + hX)/2.;
  double rr_hypoXWgt =  2./std::abs(lX - hX);

  if(!fUseUncoatedPMT && !fUseOppVolMetric)
    return rr_hypoX;

  double pe_X = pe_InvSpl.Eval(_flash_ratio);
  double pe_hX = pe_h_InvSpl.Eval(_flash_ratio);
  double pe_lX = pe_l_InvSpl.Eval(_flash_ratio);
  if(0.<pe_X && pe_X<fDriftDistance) mX = pe_X;
  else mX = (_flash_ratio<pe_means[0]) ? 0. : fDriftDistance;
  hX = (0.<pe_hX && pe_hX<fDriftDistance) ? pe_hX : mX;
  lX = (0.<pe_lX && pe_lX<fDriftDistance) ? pe_lX : mX;
  double pe_hypoX =  (lX + hX)/2.;
  double pe_hypoXWgt =  2./std::abs(lX - hX);

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
  const art::Ptr<recob::PFParticle>& pfp_ptr,
  const art::ValidHandle<std::vector<recob::PFParticle>>& pfp_h,
  std::vector<art::Ptr<recob::PFParticle>>& pfp_v)
{
  auto daughters = pfp_ptr->Daughters();
  pfp_v.push_back(pfp_ptr);
  for(auto const& daughterid : daughters) {
    if (_pfpmap.find(daughterid) == _pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error" << std::endl;
      continue;
    }
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
    AddDaughters(pfp_ptr, pfp_h, pfp_v);
  } // for all daughters
  return;
} // void FlashPredict::AddDaughters


inline
double FlashPredict::scoreTerm(double m, double n,
                               double mean, double spread)
{
  return std::abs(std::abs(m - n) - mean) / spread;
}


inline
double FlashPredict::scoreTerm(double m,
                               double mean, double spread)
{
  return std::abs(m - mean) / spread;
}


bool FlashPredict::pfpNeutrinoOnEvent(const art::ValidHandle<std::vector<recob::PFParticle>>& pfp_h)
{
  for (auto const& p : (*pfp_h)) {
    unsigned pfpPDGC = std::abs(p.PdgCode());
    if ((pfpPDGC == 12) ||
        (pfpPDGC == 14) ||
        (pfpPDGC == 16)) {
      return true;
    }
  }
  return false;
}


void FlashPredict::copyOpHitsInBeamWindow(std::vector<recob::OpHit>& opHits,
                                          art::Handle<std::vector<recob::OpHit>>& ophit_h)
{
  double s = fBeamWindowStart;
  double e = fBeamWindowEnd;
  double m = fMinOpHPE;
  // copy ophits that are inside the time window and with PEs
  auto peakInWindow =
    [s, e, m](const recob::OpHit& oph)-> bool
      {return ((oph.PeakTime() > s) &&
               (oph.PeakTime() < e) &&
               (oph.PE() > m)); };
  auto it = std::copy_if(ophit_h->begin(), ophit_h->end(), opHits.begin(),
                         peakInWindow);
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


bool FlashPredict::createOpHitsTimeHist(std::vector<recob::OpHit>& opHits)
{
  for(auto const& oph : opHits) {
    auto ch = oph.OpChannel();
    auto opDetXYZ = geometry->OpDetGeoFromOpChannel(ch).GetCenter();
    if (fICARUS &&
        !fGeoCryo->ContainsPosition(opDetXYZ)) continue;
    if(fSBND && !fUseUncoatedPMT &&
       !fPDMapAlgPtr->isPDType(oph.OpChannel(), "pmt_coated")) continue;
    fOpHitsTimeHist->Fill(oph.PeakTime(), fPEscale * oph.PE());
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
  double lowedge  = _flash_time + fLightWindowStart;
  double highedge = _flash_time + fLightWindowEnd;
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
  fOpH_end = std::partition(opHits.begin(), opHits.end(),
                            peakInsideEdges);
  fPeakCounter++;
  return true;
}


bool FlashPredict::isPDRelevant(int pdChannel,
                                std::set<unsigned>& tpcWithHits)
{
  if (fICARUS) {
    // BUG: I believe this function is not working, every now and then
    // I get ophits from the other cryo
    auto& p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
    // if the channel is in the Cryostat is relevant
    return fGeoCryo->ContainsPosition(p);
  }
  else if (fSBND) {
    // if there's hits on all TPCs all channels are relevant
    if(tpcWithHits.size() == fNTPC) return true;
    unsigned pdTPC = sbndPDinTPC(pdChannel);
    for(auto itpc: tpcWithHits){
      if(itpc == pdTPC) return true;
    }
  }
  return false;
}


// TODO: find better, less hacky solution
unsigned FlashPredict::icarusPDinTPC(int pdChannel)
{
  auto p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  if(fCryostat == 0) p.SetX((p.X() + 222.)/2. - 222.);//OpDets are outside the TPCs
  if(fCryostat == 1) p.SetX((p.X() - 222.)/2. + 222.);//OpDets are outside the TPCs
  return (geometry->PositionToTPCID(p)).TPC;
}


// TODO: find better, less hacky solution
unsigned FlashPredict::sbndPDinTPC(int pdChannel)
{
  auto p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  p.SetX(p.X()/2.);//OpDets are outside the TPCs
  return (geometry->PositionToTPCID(p)).TPC;
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
  if(bk.nopfpneutrino) m << "\tNo PFP Neutrino:  \t -" << bk.nopfpneutrino << "\n";
  if(bk.nonvalidophit) m << "\tNon Valid OpHits: \t -" << bk.nonvalidophit << "\n";
  if(bk.nullophittime) m << "\tNo OpHits in-time:\t -" << bk.nullophittime << "\n";
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
    - bk.nopfpneutrino - bk.nonvalidophit
    - bk.nullophittime;

  // account for the reasons that a particle might lack
  bk.pfp_bookkeeping = bk.pfp_to_score
    - bk.no_oph_hits - bk.no_charge
    - bk.no_flash_pe;

  if(bk.events_processed != bk.job_bookkeeping ||
     bk.scored_pfp != bk.pfp_bookkeeping)
    printBookKeeping(mf::LogWarning("FlashPredict"));
}


template <typename Stream>
void FlashPredict::printMetrics(std::string metric, int pdgc,
                                std::set<unsigned>& tpcWithHits,
                                double term,
                                Stream&& out)
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
