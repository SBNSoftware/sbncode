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
  : EDProducer{p}  // ,
    // More initializers here.
{
  produces< std::vector< anab::T0 > >();
  produces< art::Assns <recob::PFParticle, anab::T0> >();
  // fFlashProducer         = p.get<art::InputTag>("FlashProducer");
  fOpHitProducer            = p.get<std::string>("OpHitProducer", "ophit");
  fPandoraProducer          = p.get<std::string>("PandoraProducer", "pandora");
  fTrackProducer            = p.get<std::string>("TrackProducer", "pandoraTrack");
  fCaloProducer             = p.get<std::string>("CaloProducer", "pandoraCalo");
  fSpacePointProducer       = p.get<std::string>("SpacePointProducer", "pandora");
  fInputFilename            = p.get<std::string>("InputFileName", "FlashMatch/fm_metrics_sbnd.root"); // root file with score metrics
  fBeamWindowStart          = p.get<double>("BeamWindowStart", 0.0);
  fBeamWindowEnd            = p.get<double>("BeamWindowEnd", 4000.0);  // in ns
  fMinFlashPE               = p.get<double>("MinFlashPE", 0.0);
  fChargeToNPhotonsShower   = p.get<double>("ChargeToNPhotonsShower", 1.0);  // ~40000/1600
  fChargeToNPhotonsTrack    = p.get<double>("ChargeToNPhotonsTrack", 1.0);   // ~40000/1600
  fNoAvailableMetrics       = p.get<bool>("NoAvailableMetrics", false);
  fMakeTree                 = p.get<bool>("MakeTree", true);
  fUseCalo                  = p.get<bool>("UseCalo", false);
  fSelectNeutrino           = p.get<bool>("SelectNeutrino", true);
  fUseUncoatedPMT           = p.get<bool>("UseUncoatedPMT", false);
  fLightWindowStart         = p.get<double>("LightWindowStart", -0.010);  // in us w.r.t. flash time
  fLightWindowEnd           = p.get<double>("LightWindowEnd", 0.090);  // in us w.r.t flash time
  fDriftDistance            = p.get<double>("DriftDistance", 200.);
  fCryostat                 = p.get<int>("Cryostat", 0); //set =0 ot =1 for ICARUS to match reco chain selection
  fPEscale                  = p.get<double>("PEscale", 1.0);
  fTermThreshold            = p.get<double>("ThresholdTerm", 30.);

  fPDMapAlgPtr = art::make_tool<opdet::PDMapAlg>(p.get<fhicl::ParameterSet>("PDMapAlg"));

  fDetector = geometry->DetectorName();
  if(fDetector.find("sbnd") != std::string::npos) fDetector = "SBND";
  else if (fDetector.find("icarus") != std::string::npos) fDetector = "ICARUS";
  else {
      throw cet::exception("FlashPredict") << "Detector: " << fDetector
                                           << ", not supported. Stopping.\n";
  }
  fNTPC = geometry->NTPC();
  if (fDetector == "SBND" && fCryostat == 1) {
    throw cet::exception("FlashPredictSBND") << "SBND has only one cryostat. \n"
                                             << "Check Detector and Cryostat parameter." << std::endl;
  }
  else if (fDetector == "ICARUS" && fCryostat > 1) {
    throw cet::exception("FlashPredictICARUS") << "ICARUS has only two cryostats. \n"
                                               << "Check Detector and Cryostat parameter." << std::endl;
  }

  art::ServiceHandle<art::TFileService> tfs;

  int time_bins = int(500 * (fBeamWindowEnd - fBeamWindowStart));
  ophittime = tfs->make<TH1D>("ophittime", "ophittime", time_bins, fBeamWindowStart, fBeamWindowEnd); // in us
  ophittime->SetOption("HIST");
  ophittime2 = tfs->make<TH1D>("ophittime2", "ophittime2", 5 * time_bins, -5.0, +10.0); // in us
  ophittime2->SetOption("HIST");

  if (fMakeTree) {
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
    // TODO: add charge_time?
    _flashmatch_nuslice_tree->Branch("charge_x", &_charge_x, "charge_x/D");
    _flashmatch_nuslice_tree->Branch("charge_y", &_charge_y, "charge_y/D");
    _flashmatch_nuslice_tree->Branch("charge_z", &_charge_z, "charge_z/D");
    _flashmatch_nuslice_tree->Branch("charge_q", &_charge_q, "charge_q/D");
    _flashmatch_nuslice_tree->Branch("score", &_score, "score/D");
  }

  // TODO: Set a better way to run with no metrics

  // TODO: fill histos with less repetition and range for loops
  // read histograms and fill vectors for match score calculation
  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fInputFilename, fname);
  mf::LogInfo("FlashPredict") << "Opening file with metrics: " << fname;
  TFile *infile = new TFile(fname.c_str(), "READ");
  if(!infile->IsOpen()) {
    throw cet::exception("FlashPredictSBND") << "Could not find the light-charge match root file '"
                                             << fname << "'!\n";
  }
  //
  TH1 *temphisto = (TH1*)infile->Get("dy_h1");
  n_bins = temphisto->GetNbinsX();
  if (n_bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for dy " << n_bins << " bins " << std::endl;
    n_bins = 1;
    dy_means.push_back(0);
    dy_spreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= n_bins; ++ib) {
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
  n_bins = temphisto->GetNbinsX();
  if (n_bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for dz " << n_bins << " bins " << std::endl;
    n_bins = 1;
    dz_means.push_back(0);
    dz_spreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= n_bins; ++ib) {
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
  n_bins = temphisto->GetNbinsX();
  if (n_bins <= 0 || fNoAvailableMetrics) {
    std::cout << " problem with input histos for rr " << n_bins << " bins " << std::endl;
    n_bins = 1;
    rr_means.push_back(0);
    rr_spreads.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= n_bins; ++ib) {
      rr_means.push_back(temphisto->GetBinContent(ib));
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
  }
  //
  if (fDetector == "SBND" ) {
    temphisto = (TH1*)infile->Get("pe_h1");
    n_bins = temphisto->GetNbinsX();
    if (n_bins <= 0 || fNoAvailableMetrics) {
      std::cout << " problem with input histos for pe " << n_bins << " bins " << std::endl;
      n_bins = 1;
      pe_means.push_back(0);
      pe_spreads.push_back(0.001);
    }
    else {
      for (int ib = 1; ib <= n_bins; ++ib) {
        pe_means.push_back(temphisto->GetBinContent(ib));
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
    }
  }
  else if (fDetector == "ICARUS" ) {
    n_bins = 1;
    pe_means.push_back(0);
    pe_spreads.push_back(0.001);
  }
  //
  infile->Close();

  event_counter = 0;
  fill_counter = 0;
  bookkeeping = 0;

  nopfpneutrino_counter = 0;
  nullophittime_counter = 0;
  nonvalidophit_counter = 0;
  no_ophit_counter = 0;
  no_charge_counter = 0;
  multiple_fill_counter = 0;
  // Call appropriate produces<>() functions here.

  // Call appropriate consumes<>() for any products to be retrieved by this module.
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
  _score         = -9999.;
  event_counter++;

  geo::CryostatGeo geo_cryo = geometry->Cryostat(fCryostat);

  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPandoraProducer);

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
  art::Handle<std::vector<recob::OpHit> > ophit_h;
  e.getByLabel(fOpHitProducer, ophit_h);
  if(!ophit_h.isValid()) {
    mf::LogError("FlashPredict") << "No optical hits from producer module "
                                 << fOpHitProducer;
    nonvalidophit_counter++;
    updateBookKeeping();
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }
  std::vector<recob::OpHit> const& OpHitCollection(*ophit_h);
  std::vector<recob::OpHit> OpHitSubset(OpHitCollection.size());

  // copy ophits that are inside the time window and with PEs
  auto it = std::copy_if(OpHitCollection.begin(), OpHitCollection.end(), OpHitSubset.begin(),
                         [this](const recob::OpHit& oph)-> bool
                           { return ((oph.PeakTime() > fBeamWindowStart) &&
                                     (oph.PeakTime() < fBeamWindowEnd)   &&
                                     (oph.PE() > 0)); });
  OpHitSubset.resize(std::distance(OpHitSubset.begin(), it));
  // TODO: release OpHitCollection memory now

  _pfpmap.clear();
  for (size_t p=0; p<pfp_h->size(); p++) _pfpmap[pfp_h->at(p).Self()] = p;

  // get flash time
  ophittime->Reset();
  ophittime2->Reset();
  for(auto const& oph : OpHitSubset) {
    double PMTxyz[3];
    geometry->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
    if (fDetector == "SBND" && fPDMapAlgPtr->isPDType(oph.OpChannel(), "pmt_uncoated"))
      ophittime2->Fill(oph.PeakTime(), fPEscale * oph.PE());
    if (fDetector == "SBND" && !fPDMapAlgPtr->isPDType(oph.OpChannel(), "pmt_coated")) continue; // use only coated PMTs for SBND for flash_time
    if (!geo_cryo.ContainsPosition(PMTxyz)) continue;   // use only PMTs in the specified cryostat for ICARUS
    //    std::cout << "op hit " << j << " channel " << oph.OpChannel() << " time " << oph.PeakTime() << " pe " << fPEscale*oph.PE() << std::endl;

    ophittime->Fill(oph.PeakTime(), fPEscale * oph.PE());
    // double thisPE = fPEscale*oph.PE();
    // if (thisPE>1) ophittime->Fill(oph.PeakTime(),thisPE);
  }

  if (ophittime->GetEntries() <= 0 || ophittime->Integral() < fMinFlashPE) {
    mf::LogWarning("FlashPredict") << "\nOpHitTime has no entries: " << ophittime->GetEntries()
                                   << "\nor the integral: " << ophittime->Integral()
                                   << " is less than " << fMinFlashPE
                                   << "\nSkipping...";
    nullophittime_counter++;
    updateBookKeeping();
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  auto ibin =  ophittime->GetMaximumBin();
  _flash_time = (ibin * 0.002) + fBeamWindowStart; // in us
  double lowedge = _flash_time + fLightWindowStart;
  double highedge = _flash_time + fLightWindowEnd;
  mf::LogDebug("FlashPredict") << "light window " << lowedge << " " << highedge << std::endl;

  // only use optical hits around the flash time
  OpHitSubset.erase(
    std::remove_if(OpHitSubset.begin(), OpHitSubset.end(),
                   [lowedge, highedge](const recob::OpHit& oph)-> bool
                     { return ((oph.PeakTime() < lowedge) || (oph.PeakTime() > highedge)); }),
    OpHitSubset.end());

  // check if the drift volumes have OpHits
  std::array<bool, fDriftVolumes> lightInTPC = {false};
  for (size_t t=0; t<fDriftVolumes; t++){
    // TODO: lambda function should use a geometry service or a pdMap function
    auto it = std::find_if (OpHitSubset.begin(), OpHitSubset.end(),
                            [t, this](const recob::OpHit& oph)-> bool
                              { return isPDInCryoTPC(oph.OpChannel(), fCryostat, t, fDetector); });
    lightInTPC[t] = (it != OpHitSubset.end());
  }
  if(std::none_of(lightInTPC.begin(), lightInTPC.end(), [](bool v) { return v; })){
    mf::LogWarning("FlashPredict") << "No OpHits on event. Skipping...";
    no_ophit_counter++;
    updateBookKeeping();
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  // TODO: nice to this check in a single block
  bool pfpneutrino = false;
  for (size_t p=0; p<pfp_h->size(); p++) {
    unsigned pfpPDGC = std::abs(pfp_h->at(p).PdgCode());
    if ((pfpPDGC == 12) ||
        (pfpPDGC == 14) ||
        (pfpPDGC == 16)) {
      pfpneutrino = true;
      break;
    }
  }
  if (!pfpneutrino) {
    mf::LogWarning("FlashPredict") << "No pfp neutrino on event. Skipping...";
    nopfpneutrino_counter++;
    updateBookKeeping();
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  unsigned filled_tree = 0;
  std::vector<bool> qInTPC(fNTPC, false);
  // Loop over pandora pfp particles
  for (unsigned int p=0; p<pfp_h->size(); p++) {
    auto const& pfp = pfp_h->at(p);
    unsigned pfpPDGC = std::abs(pfp_h->at(p).PdgCode());
    if (!pfp.IsPrimary()) continue;
    if (fSelectNeutrino &&
        (pfpPDGC != 12) &&
        (pfpPDGC != 14) &&
        (pfpPDGC != 16) ) continue;
    std::vector<flashmatch::QCluster_t> qClusterInTPC(fNTPC);

    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);

    //  loop over all mothers and daughters, fill qCluster
    for (size_t i=0; i<pfp_ptr_v.size(); i++) {
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
      auto key = pfp_ptr_v.at(i).key();
      recob::PFParticle pfp = *pfp_ptr_v.at(i);
      pfp_v.push_back(pfp);

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
      // TODO: refactor this loop over spacepoints so that it's not
      // necessary to query the wire position every time.
      // There's just two different X wire positions on any given
      // cryostat for the collection wires.
      for (auto& SP : spacepoint_ptr_v) {
        auto const& spkey = SP.key();
        const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
        for (auto& hit : this_hit_ptr_v) {
          // Only use hits from the collection plane
          geo::WireID wid = hit->WireID();
          if (geometry->SignalType(wid) != geo::kCollection) continue;
          // Add the charged point to the vector
          const auto &position(SP->XYZ());
          const auto tpcindex = wid.TPC;

          // TODO BUG!: SBND has 2 TPC, ICARUS 4 per cryo, the next functions breaks ICARUS
          // throw the charge coming from another TPC
          // if (!isChargeInCryoTPC(position[0], fCryostat, tpcindex, fDetector)) continue;

          const auto charge(hit->Integral());
          double Wxyz[3];
          geometry->WireIDToWireGeo(wid).GetCenter(Wxyz);
          double wires_distance_X = std::abs(position[0] - Wxyz[0]);
          qClusterInTPC[tpcindex].emplace_back(
            wires_distance_X, position[1], position[2],
            charge * (lar_pandora::LArPandoraHelper::IsTrack(pfp_ptr)
                      ? fChargeToNPhotonsTrack : fChargeToNPhotonsShower));
        } // for all hits associated to this spacepoint
      } // for all spacepoints
      //      }  // if track or shower
    } // for all pfp pointers

    std::vector<double> mscore(fNTPC, 0.);
    // double charge[nMaxTPCs] = {0.}; // TODO: Use this
    for (size_t itpc=0; itpc<fNTPC; ++itpc) {

      // TODO BUG!: SBND has 2 TPC, ICARUS 4 per cryo, the next functions breaks ICARUS
      // if (!lightInTPC[itpc]) {
      //   // mf::LogDebug("FlashPredict") << "No LIGHT in the " << itpc << " TPC, continue.";
      //   continue;
      // }

      double xave = 0.0; double yave = 0.0; double zave = 0.0; double norm = 0.0;
      _charge_q = 0;
      // TODO: use accumulators instead of this for loop
      for (auto& qp : qClusterInTPC[itpc]) {
        xave += 0.001 * qp.q * qp.x;
        yave += 0.001 * qp.q * qp.y;
        zave += 0.001 * qp.q * qp.z;
        norm += 0.001 * qp.q;
        _charge_q += qp.q;
      }
      if (norm <= 0) {
        qInTPC[itpc] = false;
        mf::LogDebug("FlashPredict") << "No CHARGE in the " << itpc << " TPC, continue.";
        continue;
      }
      qInTPC[itpc] = true;
      _charge_x = xave / norm;
      _charge_y = yave / norm;
      _charge_z = zave / norm;
      // charge[itpc] = _charge_q; //TODO: Use this

      computeFlashMetrics(itpc, OpHitSubset);

      // calculate match score here, put association on the event
      double slice = _charge_x;
      _score = 0.; int icount = 0;
      double term;
      std::ostringstream thresholdMessage;
      thresholdMessage << std::left << std::setw(12) << std::setfill(' ');
      thresholdMessage << "pfp.PdgCode:\t" << pfp.PdgCode() << "\n"
                       << "_run:       \t" << _run << "\n"
                       << "_sub:       \t" << _sub << "\n"
                       << "_evt:       \t" << _evt << "\n"
                       << "itpc:       \t" << itpc << "\n"
                       << "_flash_y:   \t" << std::setw(8) << _flash_y   << ",\t"
                       << "_charge_y:  \t" << std::setw(8) << _charge_y  << "\n"
                       << "_flash_z:   \t" << std::setw(8) << _flash_z   << ",\t"
                       << "_charge_z:  \t" << std::setw(8) << _charge_z  << "\n"
                       << "_flash_x:   \t" << std::setw(8) << _flash_x   << ",\t"
                       << "_charge_x:  \t" << std::setw(8) << _charge_x  << "\n"
                       << "_flash_pe:  \t" << std::setw(8) << _flash_pe  << ",\t"
                       << "_charge_q:  \t" << std::setw(8) << _charge_q  << "\n"
                       << "_flash_r:   \t" << std::setw(8) << _flash_r   << "\n"
                       << "_flash_time: \t" << std::setw(8) << _flash_time << "\n" << std::endl;
      int isl = int(n_bins * (slice / fDriftDistance));
      if (dy_spreads[isl] > 0) {
        term = std::abs(std::abs(_flash_y - _charge_y) - dy_means[isl]) / dy_spreads[isl];
        if (term > fTermThreshold) std::cout << "\nBig term Y:\t" << term << ",\tisl:\t" << isl << "\n" << thresholdMessage.str();
        _score += term;
      }
      icount++;
      isl = int(n_bins * (slice / fDriftDistance));
      if (dz_spreads[isl] > 0) {
        term = std::abs(std::abs(_flash_z - _charge_z) - dz_means[isl]) / dz_spreads[isl];
        if (term > fTermThreshold) std::cout << "\nBig term Z:\t" << term << ",\tisl:\t" << isl << "\n" << thresholdMessage.str();
        _score += term;
      }
      icount++;
      isl = int(n_bins * (slice / fDriftDistance));
      if (rr_spreads[isl] > 0 && _flash_r > 0) {
        term = std::abs(_flash_r - rr_means[isl]) / rr_spreads[isl];
        if (term > fTermThreshold) std::cout << "\nBig term R:\t" << term << ",\tisl:\t" << isl << "\n" << thresholdMessage.str();
        _score += term;
      }
      icount++;
      if (fDetector == "SBND" && fUseUncoatedPMT) {
        isl = int(n_bins * (slice / fDriftDistance));
        double myratio = 100.0 * _flash_unpe;
        if (pe_spreads[isl] > 0 && _flash_pe > 0) {
          myratio /= _flash_pe;
          term = std::abs(myratio - pe_means[isl]) / pe_spreads[isl];
          if (term > fTermThreshold) std::cout << "\nBig term RATIO:\t" << term << ",\tisl:\t" << isl << "\n" << thresholdMessage.str();
          _score += term;
          icount++;
        }
      }
      //      _score/=icount;
      if (_flash_pe > 0 ) { // TODO: is this really the best condition?
        mscore[itpc] = _score;
        if (fMakeTree) {_flashmatch_nuslice_tree->Fill();}
        fill_counter++;
        filled_tree++;
        if(filled_tree > 1) {
          multiple_fill_counter++;
          mf::LogWarning("FlashPredict") << "Event has produced " << filled_tree << " scores";
        }
        if (!pfpneutrino) {
          mf::LogError("FlashPredict") << "No pfpneutrino on event, but filling!.";
        }
      }
    }  // end loop over TPCs


    double this_score = 0.; double score_count = 0.; // double totc = 0; //TODO: Use this
    for (size_t itpc=0; itpc<fNTPC; ++itpc) {
      this_score += mscore[itpc];
      // totc += charge[itpc];
      if (mscore[itpc] > 0) score_count++;
    }
    if (score_count > 0) {
      this_score /= score_count;

      // create t0 and pfp-t0 association here
      T0_v->push_back(anab::T0(_flash_time, icountPE, p, 0, this_score));
      util::CreateAssn(*this, e, *T0_v, pfp_ptr, *pfp_t0_assn_v);
    }
  } // over all PFparticles

  if(filled_tree == 0 &&
     std::none_of(qInTPC.begin(), qInTPC.end(),
                  [](bool v){return v;})){
    mf::LogWarning("FlashPredict") << "No charge associated to pfpneutrino. Skipping...";
    no_charge_counter++;
    updateBookKeeping();
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  updateBookKeeping();

  e.put(std::move(T0_v));
  e.put(std::move(pfp_t0_assn_v));

}// end of producer module

void FlashPredict::computeFlashMetrics(size_t itpc, std::vector<recob::OpHit> const& OpHitSubset)
{
  // store PMT photon counts in the tree as well
  double PMTxyz[3];
  double unpe_tot = 0;
  double pnorm = 0;
  double sum =    0;
  double sumy =   0; double sumz =   0;
  double sum_Ay = 0; double sum_Az = 0;
  double sum_By = 0; double sum_Bz = 0;
  double sum_Cy = 0; double sum_Cz = 0;
  double sum_D =  0;
  // TODO: change this next loop, such that it only loops
  // through channels in the current fCryostat
  for(auto const& oph : OpHitSubset) {
    std::string op_type = "pmt"; // the label ICARUS has
    if (fDetector == "SBND") op_type = fPDMapAlgPtr->pdType(oph.OpChannel());
    geometry->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
    // check cryostat and tpc

    // TODO BUG!: SBND has 2 TPC, ICARUS 4 per cryo, the next functions breaks ICARUS
    // if (!isPDInCryoTPC(PMTxyz[0], fCryostat, itpc, fDetector)) continue;

    // only use PMTs for SBND
    if (op_type == "pmt_coated" || op_type == "pmt") {
      // Add up the position, weighting with PEs
      _flash_x = PMTxyz[0];
      sum     += 1.0;
      pnorm   += oph.PE();
      sumy    += oph.PE() * PMTxyz[1];
      sumz    += oph.PE() * PMTxyz[2];
      sum_By  += PMTxyz[1];
      sum_Bz  += PMTxyz[2];
      sum_Ay  += oph.PE() * PMTxyz[1] * oph.PE() * PMTxyz[1];
      sum_Az  += oph.PE() * PMTxyz[2] * oph.PE() * PMTxyz[2];
      sum_D   += oph.PE() * oph.PE();
      sum_Cy  += oph.PE() * oph.PE() * PMTxyz[1];
      sum_Cz  += oph.PE() * oph.PE() * PMTxyz[2];
    }
    else if ( op_type == "pmt_uncoated") {
      unpe_tot += oph.PE();
    }
    else if ( (op_type == "arapuca_vuv" || op_type == "arapuca_vis") ) {
      //TODO: Use ARAPUCA
      // arape_tot+=oph.PE();
      continue;
    }
    else if ( op_type == "xarapuca_vuv" || op_type == "xarapuca_vis")  {
      //TODO: Use XARAPUCA
      // xarape_tot+=oph.PE();
      continue;
    }
  }

  if (pnorm > 0) {
    _flash_pe = pnorm * fPEscale;
    _flash_y = sum_Cy / sum_D;
    _flash_z = sum_Cz / sum_D;
    sum_By = _flash_y;
    sum_Bz = _flash_z;
    _flash_r = sqrt((sum_Ay - 2.0 * sum_By * sum_Cy + sum_By * sum_By * sum_D + sum_Az - 2.0 * sum_Bz * sum_Cz + sum_Bz * sum_Bz * sum_D) / sum_D);
    _flash_unpe = unpe_tot * fPEscale;
    icountPE = (unsigned int)(std::round(_flash_pe));
    //   std::cout << "itpc:\t" << itpc << "\n";
    //   std::cout << "_flash_pe:\t" << _flash_pe << "\n";
    //   std::cout << "_flash_y:\t" << _flash_y << "\n";
    //   std::cout << "_flash_z:\t" << _flash_z << "\n";
  }
  else {
    mf::LogWarning("FlashPredict") << "Really odd that I landed here, this shouldn't had happen.\n"
                                   << "pnorm:\t" << pnorm << "\n"
                                   << "OpHitSubset.size():\t" << OpHitSubset.size() << "\n";
    _flash_y = 0;
    _flash_z = 0;
    _flash_r = 0;
    _flash_pe = 0;
    _flash_unpe = 0;
  }
}


::flashmatch::Flash_t FlashPredict::GetFlashPESpectrum(const recob::OpFlash& opflash)
{
  // prepare container to store flash
  ::flashmatch::Flash_t flash;
  flash.time = opflash.Time();
  // geometry service
  const art::ServiceHandle<geo::Geometry> geometry;
  uint nOpDets(geometry->NOpDets());
  std::vector<double> PEspectrum;
  PEspectrum.resize(nOpDets);
  // apply gain to OpDets
  for (uint OpChannel=0; OpChannel<nOpDets; ++OpChannel) {
    uint opdet = geometry->OpDetFromOpChannel(OpChannel);
    PEspectrum[opdet] = opflash.PEs().at(OpChannel);
  }
  _pe_reco_v = PEspectrum;

  // Reset variables
  flash.x = flash.y = flash.z = 0;
  flash.x_err = flash.y_err = flash.z_err = 0;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;
  for (unsigned int opdet=0; opdet<PEspectrum.size(); opdet++) {
    double PMTxyz[3];
    geometry->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);
    // Add up the position, weighting with PEs
    sumy    += PEspectrum[opdet] * PMTxyz[1];
    sumy2   += PEspectrum[opdet] * PMTxyz[1] * PMTxyz[1];
    sumz    += PEspectrum[opdet] * PMTxyz[2];
    sumz2   += PEspectrum[opdet] * PMTxyz[2] * PMTxyz[2];
    totalPE += PEspectrum[opdet];
  }

  flash.y = sumy / totalPE;
  flash.z = sumz / totalPE;
  // This is just sqrt(<x^2> - <x>^2)
  if ((sumy2 * totalPE - sumy * sumy) > 0.)
    flash.y_err = std::sqrt(sumy2 * totalPE - sumy * sumy) / totalPE;
  if ((sumz2 * totalPE - sumz * sumz) > 0.)
    flash.z_err = std::sqrt(sumz2 * totalPE - sumz * sumz) / totalPE;
  // Set the flash properties
  flash.pe_v.resize(nOpDets);
  flash.pe_err_v.resize(nOpDets);
  // Fill the flash with the PE spectrum
  for (unsigned int i=0; i<nOpDets; ++i) {
    const auto PE(PEspectrum.at(i));
    flash.pe_v.at(i) = PE;
    flash.pe_err_v.at(i) = std::sqrt(PE);
  }
  if (flash.pe_v.size() != nOpDets)
    throw cet::exception("FlashNeutrinoId") << "Number of channels in beam flash doesn't match the number of OpDets!" << std::endl;
  return flash;

}// ::flashmatch::Flash_t FlashPredict::GetFlashPESpectrum


void FlashPredict::CollectDownstreamPFParticles(
  const lar_pandora::PFParticleMap &pfParticleMap,
  const art::Ptr<recob::PFParticle> &particle,
  lar_pandora::PFParticleVector &downstreamPFParticles) const
{
  if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end()){
    downstreamPFParticles.push_back(particle);
  }
  for (const auto &daughterId : particle->Daughters()) {
    const auto iter(pfParticleMap.find(daughterId));
    if (iter == pfParticleMap.end()){
      throw cet::exception("FlashNeutrinoId") << "Scrambled PFParticle IDs" << std::endl;
    }
    this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
  }
} // void FlashPredict::CollectDownstreamPFParticles


void FlashPredict::CollectDownstreamPFParticles(
  const lar_pandora::PFParticleMap &pfParticleMap,
  const lar_pandora::PFParticleVector &parentPFParticles,
  lar_pandora::PFParticleVector &downstreamPFParticles) const
{
  for (const auto &particle : parentPFParticles){
    this->CollectDownstreamPFParticles(pfParticleMap, particle, downstreamPFParticles);
  }
} // void FlashPredict::CollectDownstreamPFParticles


void FlashPredict::AddDaughters(
  const art::Ptr<recob::PFParticle>& pfp_ptr,
  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h,
  std::vector<art::Ptr<recob::PFParticle> > &pfp_v)
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


// TODO: no hardcoding
// TODO: collapse with the next
bool FlashPredict::isPDInCryoTPC(double pd_x, int icryo, size_t itpc, std::string detector)
{
  // check whether this optical detector views the light inside this tpc.
  std::ostringstream lostPDMessage;
  lostPDMessage << "\nThere's an " << detector << "photo detector that belongs nowhere. \n"
                << "icryo: " << icryo << "\n"
                << "itpc:  " << itpc <<  "\n"
                << "pd_x:  " << pd_x <<  std::endl;

  if (detector == "ICARUS") {
    if (icryo == 0) {
      if (itpc == 0 && -400 < pd_x && pd_x < -300 ) return true;
      else if (itpc == 1 && -100 < pd_x && pd_x < 0) return true;
      // else {std::cout << lostPDMessage.str(); return false;}
    }
    else if (icryo == 1) {
      if (itpc == 0 && 0 < pd_x && pd_x < 100) return true;
      else if (itpc == 1 && 300 < pd_x && pd_x < 400) return true;
      // else {std::cout << lostPDMessage.str(); return false;}
    }
  }
  else if (detector == "SBND") {
    if ((itpc == 0 && -213. < pd_x && pd_x < 0) || (itpc == 1 && 0 < pd_x && pd_x < 213) ) return true;
    //    else {std::cout << lostPDMessage.str(); return false;}
    else {
      return false;
    }
  }
  return false;
}

bool FlashPredict::isPDInCryoTPC(int pdChannel, int icryo, size_t itpc, std::string detector)
{
  // check whether this optical detector views the light inside this tpc.
  double dummy[3];
  geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter(dummy);
  double pd_x = dummy[0];
  return isPDInCryoTPC(pd_x, icryo, itpc, detector);
}

// TODO: no hardcoding
// TODO: collapse with the previous
// TODO: figure out what to do with the charge that falls into the crevices
bool FlashPredict::isChargeInCryoTPC(double qp_x, int icryo, int itpc, std::string detector)
{
  std::ostringstream lostChargeMessage;
  lostChargeMessage << "\nThere's " << detector << " charge that belongs nowhere. \n"
                    << "icryo: " << icryo << "\n"
                    << "itpc: "  << itpc << "\n"
                    << "qp_x: " << qp_x << std::endl;

  if (detector == "ICARUS") {
    if (icryo == 0) {
      if (itpc == 0 && -368.49 <= qp_x && qp_x <= -220.29 ) return true;
      else if (itpc == 1 && -220.14 <= qp_x && qp_x <= -71.94) return true;
      // else {std::cout << lostChargeMessage.str(); return false;}
    }
    else if (icryo == 1) {
      if (itpc == 0 && 71.94 <= qp_x && qp_x <= 220.14) return true;
      else if (itpc == 1 && 220.29 <= qp_x && qp_x <= 368.49) return true;
      // else {std::cout << lostChargeMessage.str(); return false;}
    }
  }
  else if (detector == "SBND") {
    if ((itpc == 0 && qp_x < 0) || (itpc == 1 && qp_x > 0) ) return true;
    else {
      return false;
    }
    //    else {std::cout << lostChargeMessage.str(); return false;}
  }
  return false;
}

void FlashPredict::printBookKeeping(std::string stream="info")
{
  std::ostringstream message;
  message
    // << "Bookkeeping\n"
    << "event_counter:        \t  " << event_counter << "\n"
    << "nopfpneutrino_counter:\t -" << nopfpneutrino_counter << "\n"
    << "nullophittime_counter:\t -" << nullophittime_counter << "\n"
    << "nonvalidophit_counter:\t -" << nonvalidophit_counter << "\n"
    << "no_ophit_counter:     \t -" << no_ophit_counter << "\n"
    << "no_charge_counter:    \t -" << no_charge_counter << "\n"
    << "multiple_fill_counter:\t +" << multiple_fill_counter << "\n"
    << "fill_counter:\t" << fill_counter << "\n"
    << "bookkeeping:\t" << bookkeeping << "\n";
  if(stream == "debug")
    mf::LogDebug("FlashPredict") << message.str();
  else if(stream == "info")
    mf::LogInfo("FlashPredict") << message.str();
  else if(stream == "warning")
    mf::LogWarning("FlashPredict") << message.str();
  else
    mf::LogVerbatim("FlashPredict") << message.str();
}

void FlashPredict::updateBookKeeping()
{
  bookkeeping = event_counter
    - nopfpneutrino_counter - nullophittime_counter
    - nonvalidophit_counter - no_ophit_counter
    - no_charge_counter
    + multiple_fill_counter;

  if(fill_counter != bookkeeping) printBookKeeping("warning");
}

void FlashPredict::beginJob()
{
  // Implementation of optional member function here.




}

void FlashPredict::endJob()
{
  printBookKeeping("warning");
}

DEFINE_ART_MODULE(FlashPredict)
