namespace single_photon
{

  void SinglePhoton::AnalyzeEventWeight(art::Event const & e){
    std::cout<<"SinglePhoton::AnalyzeEventWeight- starting"<<std::endl;
    art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
    e.getByLabel("generator", mcFluxHandle);
    if (!mcFluxHandle.isValid()) return;
    std::vector< art::Ptr<simb::MCFlux> > mcFluxVec;
    art::fill_ptr_vector(mcFluxVec, mcFluxHandle);
    if (mcFluxVec.size() == 0){
      std::cout << ">> No MCFlux information" << std::endl;
      return;
    }
    
    art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
    e.getByLabel("generator", mcTruthHandle);
    if (!mcTruthHandle.isValid()) return;
    std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;
    art::fill_ptr_vector(mcTruthVec, mcTruthHandle);
    if (mcTruthVec.size() == 0){
      std::cout << ">> No MCTruth information" << std::endl;
      return;
    }
    
    art::Handle< std::vector< simb::GTruth > > gTruthHandle;
    e.getByLabel("generator", gTruthHandle);
    if (!gTruthHandle.isValid()) return;
    std::vector< art::Ptr<simb::GTruth> > gTruthVec;
    art::fill_ptr_vector(gTruthVec, gTruthHandle);
    if (gTruthVec.size() == 0){
      std::cout << ">> No GTruth information" << std::endl;
      return;
    }

  const art::Ptr<simb::MCFlux> mcFlux = mcFluxVec.at(0);
  const art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);
  const simb::MCParticle& nu = mcTruth->GetNeutrino().Nu();
  const art::Ptr<simb::GTruth> gTruth = gTruthVec.at(0);

  m_run_number_eventweight = e.run();
  m_subrun_number_eventweight = e.subRun();
  m_event_number_eventweight = e.event();

  // possibly the wrong variables, but let's see for now...
  //m_mcflux_evtno     = mcFlux->fevtno;
  m_mcflux_nu_pos_x    = nu.Vx();
  m_mcflux_nu_pos_y    = nu.Vy();
  m_mcflux_nu_pos_z    = nu.Vz();
  m_mcflux_nu_mom_x    = nu.Px(); 
  m_mcflux_nu_mom_y    = nu.Py(); 
  m_mcflux_nu_mom_z    = nu.Pz(); 
  m_mcflux_nu_mom_E    = nu.E();
  m_mcflux_ntype     = mcFlux->fntype;
  m_mcflux_ptype     = mcFlux->fptype;
  m_mcflux_nimpwt    = mcFlux->fnimpwt;
  m_mcflux_dk2gen    = mcFlux->fdk2gen;
  m_mcflux_nenergyn  = mcFlux->fnenergyn;
  m_mcflux_tpx       = mcFlux->ftpx;
  m_mcflux_tpy       = mcFlux->ftpy;
  m_mcflux_tpz       = mcFlux->ftpz;
  m_mcflux_tptype    = mcFlux->ftptype;
  m_mcflux_vx        = mcFlux->fvx;
  m_mcflux_vy        = mcFlux->fvy;
  m_mcflux_vz        = mcFlux->fvz;

  // loop MCParticle info for m_mctruth object

  m_mctruth_nparticles = mcTruth->NParticles();

  for (int i = 0; i < m_mctruth_nparticles; i++){

    const simb::MCParticle& mcParticle = mcTruth->GetParticle(i);

    m_mctruth_particles_track_Id[i] = mcParticle.TrackId();
    m_mctruth_particles_pdg_code[i] = mcParticle.PdgCode();
    m_mctruth_particles_mother[i]  = mcParticle.Mother();
    m_mctruth_particles_status_code[i] = mcParticle.StatusCode();
    m_mctruth_particles_num_daughters[i] = mcParticle.NumberDaughters();

    for (int j = 0; j < m_mctruth_particles_num_daughters[i]; j++){

      const simb::MCParticle& daughterMcParticle = mcTruth->GetParticle(j);
      m_mctruth_particles_daughters[i][j] = daughterMcParticle.TrackId();

    }

    m_mctruth_particles_Gvx[i] = mcParticle.Gvx();
    m_mctruth_particles_Gvy[i] = mcParticle.Gvy();
    m_mctruth_particles_Gvz[i] = mcParticle.Gvz();
    m_mctruth_particles_Gvt[i] = mcParticle.Gvt();
    m_mctruth_particles_px0[i] = mcParticle.Px(0);
    m_mctruth_particles_py0[i] = mcParticle.Py(0);
    m_mctruth_particles_pz0[i] = mcParticle.Pz(0);
    m_mctruth_particles_e0[i] = mcParticle.E(0);
    m_mctruth_particles_rescatter[i] = mcParticle.Rescatter();
    m_mctruth_particles_polx[i] = mcParticle.Polarization().X();
    m_mctruth_particles_poly[i] = mcParticle.Polarization().Y();
    m_mctruth_particles_polz[i] = mcParticle.Polarization().Z();
  }

  const simb::MCNeutrino& mcNeutrino = mcTruth->GetNeutrino();

  m_mctruth_neutrino_ccnc = mcNeutrino.CCNC();
  m_mctruth_neutrino_mode = mcNeutrino.Mode();
  m_mctruth_neutrino_interaction_type = mcNeutrino.InteractionType();
  m_mctruth_neutrino_target = mcNeutrino.Target();
  m_mctruth_neutrino_nucleon = mcNeutrino.HitNuc();
  m_mctruth_neutrino_quark = mcNeutrino.HitQuark();
  m_mctruth_neutrino_w = mcNeutrino.W();
  m_mctruth_neutrino_x = mcNeutrino.X();
  m_mctruth_neutrino_y = mcNeutrino.Y();
  m_mctruth_neutrino_qsqr = mcNeutrino.QSqr();

  m_gtruth_is_sea_quark = gTruth->fIsSeaQuark;
  m_gtruth_tgt_pdg = gTruth->ftgtPDG;
  m_gtruth_tgt_A = gTruth->ftgtA;
  m_gtruth_tgt_Z = gTruth->ftgtZ;
  m_gtruth_tgt_p4_x = gTruth->fTgtP4.X();
  m_gtruth_tgt_p4_y = gTruth->fTgtP4.Y();
  m_gtruth_tgt_p4_z = gTruth->fTgtP4.Z();
  m_gtruth_tgt_p4_E = gTruth->fTgtP4.E();

  m_gtruth_weight = gTruth->fweight;
  m_gtruth_probability = gTruth->fprobability;
  m_gtruth_xsec = gTruth->fXsec;
  m_gtruth_diff_xsec = gTruth->fDiffXsec;
  m_gtruth_gphase_space = gTruth->fGPhaseSpace;

  m_gtruth_vertex_x = gTruth->fVertex.X();
  m_gtruth_vertex_y = gTruth->fVertex.Y();
  m_gtruth_vertex_z = gTruth->fVertex.Z();
  m_gtruth_vertex_T = gTruth->fVertex.T();
  m_gtruth_gscatter = gTruth->fGscatter;
  m_gtruth_gint = gTruth->fGint;
  m_gtruth_res_num = gTruth->fResNum;
  m_gtruth_num_piplus = gTruth->fNumPiPlus;
  m_gtruth_num_pi0 = gTruth->fNumPi0;
  m_gtruth_num_piminus = gTruth->fNumPiMinus;
  m_gtruth_num_proton = gTruth->fNumProton;
  m_gtruth_num_neutron = gTruth->fNumNeutron;
  m_gtruth_is_charm = gTruth->fIsCharm;
  m_gtruth_is_strange = gTruth->fIsStrange;
  m_gtruth_decay_mode = gTruth->fDecayMode;
  m_gtruth_strange_hadron_pdg = gTruth->fStrangeHadronPdg;
  m_gtruth_charm_hadron_pdg = gTruth->fCharmHadronPdg;
  m_gtruth_gx = gTruth->fgX;
  m_gtruth_gy = gTruth->fgY;
  m_gtruth_gt = gTruth->fgT;
  m_gtruth_gw = gTruth->fgW;
  m_gtruth_gQ2 = gTruth->fgQ2;
  m_gtruth_gq2 = gTruth->fgq2;
  m_gtruth_probe_pdg = gTruth->fProbePDG;
  m_gtruth_probe_p4_x = gTruth->fProbeP4.X();
  m_gtruth_probe_p4_y = gTruth->fProbeP4.Y();
  m_gtruth_probe_p4_z = gTruth->fProbeP4.Z();
  m_gtruth_probe_p4_E = gTruth->fProbeP4.E();
  m_gtruth_hit_nuc_p4_x = gTruth->fHitNucP4.X();
  m_gtruth_hit_nuc_p4_y = gTruth->fHitNucP4.Y();
  m_gtruth_hit_nuc_p4_z = gTruth->fHitNucP4.Z();
  m_gtruth_hit_nuc_p4_E = gTruth->fHitNucP4.E();
  m_gtruth_hit_nuc_pos = gTruth->fHitNucPos;
  m_gtruth_fs_had_syst_p4_x = gTruth->fFShadSystP4.X();
  m_gtruth_fs_had_syst_p4_y = gTruth->fFShadSystP4.Y();
  m_gtruth_fs_had_syst_p4_z = gTruth->fFShadSystP4.Z();
  m_gtruth_fs_had_syst_p4_E = gTruth->fFShadSystP4.E();
 
  //moved to inside singlphoontmodule.cc for filter reasons
  //eventweight_tree->Fill();
  std::cout<<"SinglePhoton::AnalyzeEventWeight-eventweight_tree filled"<<std::endl;
  }
}
