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
  // , fCaloProducer(p.get<std::string>("CaloProducer"))
  // , fTrackProducer(p.get<std::string>("TrackProducer"))
  , fBeamWindowStart(p.get<double>("BeamWindowStart")) // TODO: should come from service
  , fBeamWindowEnd(p.get<double>("BeamWindowEnd"))// in ns // TODO: should come from service
  , fLightWindowStart(p.get<double>("LightWindowStart")) // in us w.r.t. flash time
  , fLightWindowEnd(p.get<double>("LightWindowEnd"))  // in us w.r.t flash time
  , fTimeBins(unsigned(500 * (fBeamWindowEnd - fBeamWindowStart)))
  , fSelectNeutrino(p.get<bool>("SelectNeutrino", true))
  , fUseUncoatedPMT(p.get<bool>("UseUncoatedPMT", false))
  // , fUseCalo(p.get<bool>("UseCalo", false))
  , fInputFilename(p.get<std::string>("InputFileName")) // root file with score metrics
  , fNoAvailableMetrics(p.get<bool>("NoAvailableMetrics", false))
  , fMakeTree(p.get<bool>("MakeTree", false))
  , fMinFlashPE(p.get<double>("MinFlashPE", 0.0))
  , fPEscale(p.get<double>("PEscale", 1.0))
  , fChargeToNPhotonsShower(p.get<double>("ChargeToNPhotonsShower", 1.0))  // ~40000/1600
  , fChargeToNPhotonsTrack(p.get<double>("ChargeToNPhotonsTrack", 1.0))  // ~40000/1600
  , fCryostat(p.get<int>("Cryostat", 0)) //set =0 ot =1 for ICARUS to match reco chain selection
  , fDriftDistance(p.get<double>("DriftDistance"))// TODO: should come from geometry
  , fVUVToVIS(p.get<unsigned>("VUVToVIS", 4))
  , fTermThreshold(p.get<double>("ThresholdTerm", 30.))
{
  produces< std::vector< anab::T0 > >();
  produces< art::Assns <recob::PFParticle, anab::T0> >();
  // fFlashProducer         = p.get<art::InputTag>("FlashProducer");

  fPDMapAlgPtr = art::make_tool<opdet::PDMapAlg>(p.get<fhicl::ParameterSet>("PDMapAlg"));
  fGeoCryo = std::make_unique<geo::CryostatGeo>(geometry->Cryostat(fCryostat));
  fNTPC = geometry->NTPC();

  fDetector = geometry->DetectorName();
  if(fDetector.find("sbnd") != std::string::npos) {
    fDetector = "SBND";
    fSBND = true;
    fICARUS = false;
  }
  else if (fDetector.find("icarus") != std::string::npos) {
    fDetector = "ICARUS";
    fSBND = false;
    fICARUS = true;
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
  _score         = -9999.;
  bk.events++;

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
  copyOpHitsInWindow(opHits, ophit_h);

  {// TODO: pack this into a function
  // get flash time
  TH1D ophittime("ophittime", "ophittime", fTimeBins, fBeamWindowStart, fBeamWindowEnd); // in us
  ophittime.SetOption("HIST");
  TH1D ophittime2("ophittime2", "ophittime2", 5 * fTimeBins, -5.0, +10.0); // in us
  ophittime2.SetOption("HIST");

  for(auto const& oph : opHits) {
    auto ch = oph.OpChannel();
    auto opDetXYZ = geometry->OpDetGeoFromOpChannel(ch).GetCenter();
    if (fSBND && fPDMapAlgPtr->isPDType(ch, "pmt_uncoated"))
      ophittime2.Fill(oph.PeakTime(), fPEscale * oph.PE());
    if (fSBND && !fPDMapAlgPtr->isPDType(ch, "pmt_coated")) continue; // use only coated PMTs for SBND for flash_time
    if (fICARUS && !fGeoCryo->ContainsPosition(opDetXYZ)) continue; // use only PMTs in the specified cryostat for ICARUS
    ophittime.Fill(oph.PeakTime(), fPEscale * oph.PE());
  }

  if (ophittime.GetEntries() <= 0 || ophittime.Integral() < fMinFlashPE) {
    mf::LogWarning("FlashPredict")
      << "\nOpHitTime has no entries: " << ophittime.GetEntries()
      << "\nor the integral: " << ophittime.Integral()
      << " is less than " << fMinFlashPE
      << "\nSkipping...";
    bk.nullophittime++;
    updateBookKeeping();
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }

  auto ibin =  ophittime.GetMaximumBin();
  _flash_time = (ibin * 0.002) + fBeamWindowStart; // in us // TODO: hardcoding
  double lowedge = _flash_time + fLightWindowStart;
  double highedge = _flash_time + fLightWindowEnd;
  mf::LogDebug("FlashPredict") << "light window " << lowedge << " " << highedge << std::endl;

  auto peakOutsideEdges =
    [lowedge, highedge](const recob::OpHit& oph)-> bool
      { return ((oph.PeakTime() < lowedge) || (oph.PeakTime() > highedge)); };
  // only use optical hits around the flash time
  opHits.erase(
    std::remove_if(opHits.begin(), opHits.end(),
                   peakOutsideEdges),
    opHits.end());
  }// TODO: pack this into a function

  std::set<unsigned> tpcWithOpH;
  for(auto& op: opHits){
    tpcWithOpH.insert(sbndPDinTPC(op.OpChannel()));
    if(tpcWithOpH.size() == fNTPC) break;
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
    flashmatch::QCluster_t qClustsGl;
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
          const auto& position(SP->XYZ());
          auto itpc = wid.TPC;
          tpcWithHits.insert(itpc);
          const auto charge(hit->Integral());
          auto wXYZ = geometry->WireIDToWireGeo(wid).GetCenter();
          double wires_distance_X = std::abs(position[0] - wXYZ.X());
          qClusters.emplace_back(wires_distance_X, position[1], position[2], charge * chargeToNPhotons);
          qClustsGl.emplace_back(position[0], position[1], position[2], charge * chargeToNPhotons);
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
        continue;
      }
    }

    if(!computeChargeMetrics(qClusters, qClustsGl)){
      mf::LogWarning("FlashPredict") << "Clusters with No Charge. Skipping...";
      bk.no_charge++;
      continue;
    }

    if(!computeFlashMetrics(tpcWithHits, opHits)){
      printMetrics("ERROR", pfp.PdgCode(), tpcWithHits, 0, mf::LogError("FlashPredict"));
      bk.no_flash_pe++;
      continue;
    }

    if(computeScore(tpcWithHits, pfp.PdgCode())){
      if (fMakeTree) {_flashmatch_nuslice_tree->Fill();}
      bk.scored_pfp++;
      mf::LogDebug("FlashPredict") << "Creating T0 and PFP-T0 association";
      T0_v->push_back(anab::T0(_flash_time, icountPE, p, 0, _score));
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
  _flashmatch_nuslice_tree->Branch("score", &_score, "score/D");
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
  if (fSBND) {
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
  else if (fICARUS ) {
    n_bins = 1;
    pe_means.push_back(0);
    pe_spreads.push_back(0.001);
  }
  //
  infile->Close();
  delete infile;

}


bool FlashPredict::computeChargeMetrics(flashmatch::QCluster_t& qClusters,
                                        flashmatch::QCluster_t& qClustsGl)
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
  double xGl = 0.;
  for (auto& qp : qClustsGl) {
    xGl += scale * qp.q * qp.x;
  }
  if (norm > 0) {
    _charge_x_gl = xGl  / norm;
    _charge_x = xave / norm;
    _charge_y = yave / norm;
    _charge_z = zave / norm;
    return true;
  }
  return false;
}


bool FlashPredict::computeFlashMetrics(std::set<unsigned>& tpcWithHits,
                                       std::vector<recob::OpHit> const& opHits)
{
  double sum = 0.;
  double sum_PE = 0.;
  double sum_PE2 = 0.;
  double sum_unPE = 0.;
  double sum_PE2Y  = 0.; double sum_PE2Z  = 0.;
  double sum_PE2Y2 = 0.; double sum_PE2Z2 = 0.;

  double maxPE = -9999.;
  for(auto const& oph : opHits) {
    if(!isPDRelevant(oph.OpChannel(), tpcWithHits)) continue;
    auto opDet = geometry->OpDetGeoFromOpChannel(oph.OpChannel());
    auto opDetXYZ = opDet.GetCenter();

    std::string op_type = "pmt"; // the label ICARUS has
    if (fSBND) op_type = fPDMapAlgPtr->pdType(oph.OpChannel());

    if (op_type == "pmt_coated" || op_type == "pmt") {
      // _flash_x is the X coord of the opdet where most PE are deposited
      if (oph.PE()>maxPE){
        _flash_x = opDetXYZ.X();
        maxPE = oph.PE();
      }
      double ophPE2 = oph.PE() * oph.PE();
      sum       += 1.0;
      sum_PE    += oph.PE();
      sum_PE2   += ophPE2;
      sum_PE2Y  += ophPE2 * opDetXYZ.Y();
      sum_PE2Z  += ophPE2 * opDetXYZ.Z();
      sum_PE2Y2 += ophPE2 * opDetXYZ.Y() * opDetXYZ.Y();
      sum_PE2Z2 += ophPE2 * opDetXYZ.Z() * opDetXYZ.Z();
    }
    else if (op_type == "pmt_uncoated") {
      sum_unPE += oph.PE();
    }
    else if (op_type == "arapuca_vuv" || op_type == "xarapuca_vuv" ) {
      //TODO: Use ARAPUCA
      // arape_tot+=oph.PE();
      continue;
    }
    else if (op_type == "arapuca_vis" || op_type == "xarapuca_vis")  {
      //TODO: Use XARAPUCA
      // xarape_tot+=oph.PE();
      continue;
    }
  }

  if (sum_PE > 0) {
    _flash_pe    = sum_PE   * fPEscale;
    _flash_unpe  = sum_unPE * fPEscale;
    _flash_ratio = fVUVToVIS * _flash_unpe / _flash_pe;
    _flash_y  = sum_PE2Y / sum_PE2;
    _flash_z  = sum_PE2Z / sum_PE2;
    _flash_r = std::sqrt(
      (sum_PE2Y2 - 2.0 * _flash_y * sum_PE2Y + _flash_y * _flash_y * sum_PE2 +
       sum_PE2Z2 - 2.0 * _flash_z * sum_PE2Z + _flash_z * _flash_z * sum_PE2) / sum_PE2);
    if (fSBND && fUseUncoatedPMT) icountPE = std::round(_flash_pe + _flash_unpe);
    else icountPE = std::round(_flash_pe);
    return true;
  }
  else {
    std::string channels;
    for(auto& op: opHits) channels += std::to_string(op.OpChannel()) + ' ';
    std::string tpcs;
    for(auto itpc: tpcWithHits) tpcs += std::to_string(itpc) + ' ';
    mf::LogError("FlashPredict")
      << "Really odd that I landed here, this shouldn't had happen.\n"
      << "sum:          \t" << sum << "\n"
      << "sum_PE:       \t" << sum_PE << "\n"
      << "sum_unPE:     \t" << sum_unPE << "\n"
      << "tpcWithHits:  \t" << tpcs << "\n"
      << "opHits.size():\t" << opHits.size() << "\n"
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
  int isl = int(n_bins * (_charge_x / fDriftDistance));
  auto out = mf::LogWarning("FlashPredict");

  if (dy_spreads[isl] > 0) {
    double term = scoreTerm(_flash_y, _charge_y, dy_means[isl], dy_spreads[isl]);
    if (term > fTermThreshold) printMetrics("Y", pdgc, tpcWithHits, term, out);
    _score += term;
    tcount++;
  }
  if (dz_spreads[isl] > 0) {
    double term = scoreTerm(_flash_z, _charge_z, dz_means[isl], dz_spreads[isl]);
    if (term > fTermThreshold) printMetrics("Z", pdgc, tpcWithHits, term, out);
    _score += term;
    tcount++;
  }
  if (rr_spreads[isl] > 0 && _flash_r > 0) {
    double term = scoreTerm(_flash_r, rr_means[isl], rr_spreads[isl]);
    if (term > fTermThreshold) printMetrics("R", pdgc, tpcWithHits, term, out);
    _score += term;
    tcount++;
  }
  if (fSBND && fUseUncoatedPMT) {
    if (pe_spreads[isl] > 0 && _flash_ratio > 0) {
      double term = scoreTerm(_flash_ratio, pe_means[isl], pe_spreads[isl]);
      if (term > fTermThreshold) printMetrics("RATIO", pdgc, tpcWithHits, term, out);
      _score += term;
      tcount++;
    }
  }
  mf::LogDebug("FlashPredict")
    << "score:\t" << _score << "using " << tcount << " terms";
  if(_score > 0.) return true;
  else return false;
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
  const lar_pandora::PFParticleMap& pfParticleMap,
  const art::Ptr<recob::PFParticle>& particle,
  lar_pandora::PFParticleVector& downstreamPFParticles) const
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
  const lar_pandora::PFParticleMap& pfParticleMap,
  const lar_pandora::PFParticleVector& parentPFParticles,
  lar_pandora::PFParticleVector& downstreamPFParticles) const
{
  for (const auto &particle : parentPFParticles){
    this->CollectDownstreamPFParticles(pfParticleMap, particle, downstreamPFParticles);
  }
} // void FlashPredict::CollectDownstreamPFParticles


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

void FlashPredict::copyOpHitsInWindow(std::vector<recob::OpHit>& opHits,
                                      art::Handle<std::vector<recob::OpHit>>& ophit_h)
{
  double s = fBeamWindowStart;
  double e = fBeamWindowEnd;
  // copy ophits that are inside the time window and with PEs
  auto peakInWindow =
    [s, e](const recob::OpHit& oph)-> bool
      {return ((oph.PeakTime() > s) &&
               (oph.PeakTime() < e) &&
               (oph.PE() > 0)); };
  auto it = std::copy_if(ophit_h->begin(), ophit_h->end(), opHits.begin(),
                         peakInWindow);
  opHits.resize(std::distance(opHits.begin(), it));
}


bool FlashPredict::isPDRelevant(int pdChannel,
                                std::set<unsigned>& tpcWithHits)
{
  if (fICARUS) {
    auto& p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
    // if the channel is in the Cryostat is relevant
    return fGeoCryo->ContainsPosition(p);
  }
  else if (fSBND) {
    // if there's hits on all TPCs all channels are relevant
    if(tpcWithHits.size() == fNTPC) return true;
    for(auto itpc: tpcWithHits){
      if(itpc == sbndPDinTPC(pdChannel)) return true;
    }
  }
  return false;
}


// TODO: find better, less hacky solution
unsigned FlashPredict::sbndPDinTPC(int pdChannel)
{
  auto p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  p.SetX(p.X()/2.);//OpDets are outside the TPCs
  return (geometry->PositionToTPCID(p)).TPC;
}


// TODO: no hardcoding
// TODO: collapse with the next
bool FlashPredict::isPDInCryoTPC(double pd_x, size_t itpc)
{
  // check whether this optical detector views the light inside this tpc.
  std::ostringstream lostPDMessage;
  lostPDMessage << "\nThere's a " << fDetector << " photo detector that belongs nowhere. \n"
                << "icryo: " << fCryostat << "\n"
                << "itpc:  " << itpc <<  "\n"
                << "pd_x:  " << pd_x <<  std::endl;
  if (fICARUS) {
    if (fCryostat == 0) {
      if (itpc == 0 && -400. < pd_x && pd_x < -300. ) return true;
      else if (itpc == 1 && -100. < pd_x && pd_x < 0.) return true;
      else {std::cout << lostPDMessage.str(); return false;}
    }
    else if (fCryostat == 1) {
      if (itpc == 0 && 0. < pd_x && pd_x < 100.) return true;
      else if (itpc == 1 && 300. < pd_x && pd_x < 400.) return true;
      else {std::cout << lostPDMessage.str(); return false;}
    }
  }
  else if (fSBND) {
    if (itpc == 0 && -213. < pd_x && pd_x < 0.) return true;
    else if (itpc == 1 && 0. < pd_x && pd_x < 213.) return true;
    else {std::cout << lostPDMessage.str(); return false;}
  }
  return false;
}

bool FlashPredict::isPDInCryoTPC(int pdChannel, size_t itpc)
{
  // check whether this optical detector views the light inside this tpc.
  auto p = geometry->OpDetGeoFromOpChannel(pdChannel).GetCenter();
  return isPDInCryoTPC(p.X(), itpc);
}

// TODO: no hardcoding
// TODO: collapse with the previous
// TODO: figure out what to do with the charge that falls into the crevices
bool FlashPredict::isChargeInCryoTPC(double qp_x, int icryo, int itpc)
{
  std::ostringstream lostChargeMessage;
  lostChargeMessage << "\nThere's " << fDetector << " charge that belongs nowhere. \n"
                    << "icryo: " << fCryostat << "\n"
                    << "itpc: "  << itpc << "\n"
                    << "qp_x: " << qp_x << std::endl;

  if (fICARUS) {
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
  else if (fSBND) {
    if ((itpc == 0 && qp_x < 0) || (itpc == 1 && qp_x > 0) ) return true;
    else {
      return false;
    }
    //    else {std::cout << lostChargeMessage.str(); return false;}
  }
  return false;
}


template <typename Stream>
void FlashPredict::printBookKeeping(Stream&& out)
{
  std::ostringstream m;
  m << "Book Keeping\n";
  m << "-----------------------------------\n"
    << "Job tally\n"
    << "\tevents:       \t  " << bk.events << "\n";
  if(bk.nopfpneutrino) m << "\tnopfpneutrino:\t -" << bk.nopfpneutrino << "\n";
  if(bk.nonvalidophit) m << "\tnonvalidophit:\t -" << bk.nonvalidophit << "\n";
  if(bk.nullophittime) m << "\tnullophittime:\t -" << bk.nullophittime << "\n";
  m << "\t-------------------\n";
  if(bk.job_bookkeeping != bk.events_processed)
    m << "\tjob_bookkeeping:  \t" << bk.job_bookkeeping << "\n";
  m << "\tevents_processed: \t" << bk.events_processed << "\n"
    << "-----------------------------------\n"
    << "pfp tally\n"
    << "\tpfp to score: \t  " << bk.pfp_to_score << "\n";
  if(bk.no_oph_hits) m << "\tno_oph_hits:  \t -" << bk.no_oph_hits << "\n";
  if(bk.no_charge)   m << "\tno_charge:    \t -" << bk.no_charge << "\n";
  if(bk.no_flash_pe) m << "\tno_flash_pe:  \t -" << bk.no_flash_pe << "ERROR!\n";
  m << "\t-------------------\n";
  if(bk.pfp_bookkeeping != bk.scored_pfp)
    m << "\tpfp_bookkeeping:  \t" << bk.pfp_bookkeeping << "\n";
  m << "\tscored_pfp_:      \t" << bk.scored_pfp << "\n"
    << "-----------------------------------";
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
  int isl = int(n_bins * (_charge_x / fDriftDistance));
  std::string tpcs;
  for(auto itpc: tpcWithHits) tpcs += std::to_string(itpc) + ' ';
  out
    << "Big term " << metric << ":\t" << term << "\n"
    << std::left << std::setw(12) << std::setfill(' ')
    << "isl:        \t" << isl << "\n"
    << "pfp.PdgCode:\t" << pdgc << "\n"
    << "_run:       \t" << _run << "\n"
    << "_sub:       \t" << _sub << "\n"
    << "_evt:       \t" << _evt << "\n"
    << "tpc:        \t" << tpcs << "\n"
    << "_flash_time:\t" << std::setw(8) << _flash_time << "\n"
    << "_charge_q:  \t" << std::setw(8) << _charge_q   << "\n"
    << "_flash_pe:  \t" << std::setw(8) << _flash_pe   << ",\t"
    << "_flash_unpe:\t" << std::setw(8) << _flash_unpe << "\n"
    << "_flash_y:   \t" << std::setw(8) << _flash_y    << ",\t"
    << "_charge_y:  \t" << std::setw(8) << _charge_y   << "\n"
    << "_flash_z:   \t" << std::setw(8) << _flash_z    << ",\t"
    << "_charge_z:  \t" << std::setw(8) << _charge_z   << "\n"
    << "_flash_x:   \t" << std::setw(8) << _flash_x    << ",\t"
    << "_charge_x:  \t" << std::setw(8) << _charge_x   << "\n"
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
