namespace single_photon
{
  void SinglePhoton::ClearEventWeightBranches(){
    m_mcflux_nu_pos_x=-9999;
    m_mcflux_nu_pos_y=-9999;
    m_mcflux_nu_pos_z=-9999;
    m_mcflux_nu_mom_x=-9999;
    m_mcflux_nu_mom_y=-9999;
    m_mcflux_nu_mom_z=-9999;
    m_mcflux_nu_mom_z=-9999;
    m_mcflux_nu_mom_E=-9999;
    m_mcflux_ntype=0;
    m_mcflux_ptype=0;
    m_mcflux_nimpwt=-9999;
    m_mcflux_dk2gen=-9999;
    m_mcflux_nenergyn=-9999;
    m_mcflux_tpx=-9999;
    m_mcflux_tpy=-9999;
    m_mcflux_tpz=-9999;
    m_mcflux_vx=-9999;
    m_mcflux_vy=-9999;
    m_mcflux_vz=-9999;
    m_mcflux_tptype=0;
    m_mctruth_nparticles=0;
    /*
    //  m_mctruth_particles_track_ID[];
    //  m_mctruth_particles_pdg_code[];
    //  m_mctruth_particles_mother[];
    //  m_mctruth_particles_status_code[];
    //  m_mctruth_particles_num_daughters[]; //other similar variables
    //  m_mctruth_particles_daughters[];
    //m_mctruth_particles_Gvx.clear();
    //m_mctruth_particles_Gvy.clear();
    m_mctruth_particles_Gvz.clear();
    m_mctruth_particles_Gvt.clear();
    m_mctruth_particles_px0.clear();
    m_mctruth_particles_py0.clear();
    m_mctruth_particles_pz0.clear();
    m_mctruth_particles_e0.clear();
    //int m_mctruth_particles_rescatter.clear();
    m_mctruth_particles_polx.clear();
    m_mctruth_particles_poly.clear();
    m_mctruth_particles_polz.clear();
    
    //int m_mctruth_neutrino_CCNC;
    //int m_mctruth_neutrino_mode: "m_mctruth_neutrino_mode" //declared in mctruth vars
    //m_mctruth_neutrino_interactionType: "m_mctruth_neutrino_interactionType" 
	    int m_mctruth_neutrino_target.clear();
	    int m_mctruth_neutrino_nucleon.clear();
	    int m_mctruth_neutrino_quark.clear();
	    m_mctruth_neutrino_w.clear();
	    m_mctruth_neutrino_x.clear();
	    m_mctruth_neutrino_y.clear();
    */
	    //m_mctruth_neutrino_QSqr: "m_mctruth_neutrino_QSqr"
	    m_gtruth_is_sea_quark=false;
	    m_gtruth_tgt_pdg=0;
	    m_gtruth_tgt_Z = -9999;
        m_gtruth_tgt_A = -9999;
        m_gtruth_tgt_p4_x = -9999;
        m_gtruth_tgt_p4_y = -9999;
        m_gtruth_tgt_p4_z = -9999;
        m_gtruth_tgt_p4_E = -9999;
        m_gtruth_weight=-9999;
	    m_gtruth_probability=-9999;
	    m_gtruth_xsec=-9999;
	    m_gtruth_diff_xsec=-9999;
	    m_gtruth_gphase_space=-9999;
	    m_gtruth_vertex_x=-9999;
	    m_gtruth_vertex_y=-9999;
	    m_gtruth_vertex_z=-9999;
	    m_gtruth_vertex_T=-9999;
	    m_gtruth_gscatter=-9999;
	    m_gtruth_gint=-9999;
	    m_gtruth_res_num=-9999;
	    m_gtruth_num_piplus=-9999;
	    m_gtruth_num_pi0=-9999;
	    m_gtruth_num_piminus=-9999;
	    m_gtruth_num_proton=-9999;
	    m_gtruth_num_neutron=-9999;
	    m_gtruth_is_charm=false;
	    m_gtruth_is_strange=false;
        m_gtruth_charm_hadron_pdg = -9999;
        m_gtruth_strange_hadron_pdg = -9999;
        m_gtruth_decay_mode = -9999;
	    m_gtruth_gx=-9999;
	    m_gtruth_gy=-9999;
	    m_gtruth_gy=-9999;
	    m_gtruth_gt=-9999;
	    m_gtruth_gw=-9999;
	    m_gtruth_gQ2=-9999;
	    m_gtruth_gq2=-9999;
	    m_gtruth_probe_pdg=0;
	    m_gtruth_probe_p4_x=-9999;
	    m_gtruth_probe_p4_y=-9999;
	    m_gtruth_probe_p4_z=-9999;
	    m_gtruth_probe_p4_E=-9999;
	    m_gtruth_hit_nuc_p4_x=-9999;
	    m_gtruth_hit_nuc_p4_y=-9999;
	    m_gtruth_hit_nuc_p4_z=-9999;
	    m_gtruth_hit_nuc_p4_E=-9999;
	    m_gtruth_hit_nuc_pos=-9999;
	    m_gtruth_fs_had_syst_p4_x=-9999;
	    m_gtruth_fs_had_syst_p4_y=-9999;
	    m_gtruth_fs_had_syst_p4_z=-9999;
	    m_gtruth_fs_had_syst_p4_E=-9999;
    
  }
  void SinglePhoton::CreateEventWeightBranches(){
    std::cout<<"SinglePhoton:analyze_Eventweigh:eventweight_tree make branches start"<<std::endl;
    //-----------------run info
    eventweight_tree->Branch("run", &m_run_number_eventweight); 
    eventweight_tree->Branch("subrun", &m_subrun_number_eventweight);
    eventweight_tree->Branch("event",  &m_event_number_eventweight);
    //------------------mcflux
    eventweight_tree->Branch("MCFlux_NuPosX",  &m_mcflux_nu_pos_x );
    eventweight_tree->Branch("MCFlux_NuPosY",  &m_mcflux_nu_pos_y );
    eventweight_tree->Branch("MCFlux_NuPosZ",  &m_mcflux_nu_pos_z );
    eventweight_tree->Branch("MCFlux_NuMomX",  &m_mcflux_nu_mom_x);
    eventweight_tree->Branch("MCFlux_NuMomY",  &m_mcflux_nu_mom_y );
    eventweight_tree->Branch("MCFlux_NuMomZ",  &m_mcflux_nu_mom_z);
    eventweight_tree->Branch("MCFlux_NuMomE",  &m_mcflux_nu_mom_E);
    eventweight_tree->Branch("MCFlux_ntype",  &m_mcflux_ntype );
    eventweight_tree->Branch("MCFlux_ptype",  &m_mcflux_ptype );
    eventweight_tree->Branch("MCFlux_nimpwt",  &m_mcflux_nimpwt );
    eventweight_tree->Branch("MCFlux_dk2gen",  &m_mcflux_dk2gen );
    eventweight_tree->Branch("MCFlux_nenergyn",  &m_mcflux_nenergyn); 
    eventweight_tree->Branch("MCFlux_tpx",  &m_mcflux_tpx );
    eventweight_tree->Branch("MCFlux_tpy",  &m_mcflux_tpy );
    eventweight_tree->Branch("MCFlux_tpz",  &m_mcflux_tpz );
    eventweight_tree->Branch("MCFlux_vx",  &m_mcflux_vx);
    eventweight_tree->Branch("MCFlux_vy",  &m_mcflux_vy );
    eventweight_tree->Branch("MCFlux_vz",  &m_mcflux_vz );
    eventweight_tree->Branch("MCFlux_tptype",  & m_mcflux_tptype );
    //---------------mctruth
    eventweight_tree->Branch("MCTruth_NParticles",  &m_mctruth_nparticles); 
    //eventweight_tree->Branch("MCTruth_particles_TrackId",  &m_mctruth_particles_track_Id, "m_mctruth_particles_track_Id[m_mctruth_nparticles]/I" );
    eventweight_tree->Branch("MCTruth_particles_TrackId",  &m_mctruth_particles_track_Id, "MCTruth_particles_TrackId[MCTruth_NParticles]/I" );
    //eventweight_tree->Branch("MCTruth_particles_TrackId",  &m_mctruth_particles_track_Id);
    eventweight_tree->Branch("MCTruth_particles_PdgCode",  &m_mctruth_particles_pdg_code, "MCTruth_particles_PdgCode[MCTruth_NParticles]/I" );
    eventweight_tree->Branch("MCTruth_particles_Mother",  &m_mctruth_particles_mother, "MCTruth_particles_Mother[MCTruth_NParticles]/I" );
    eventweight_tree->Branch("MCTruth_particles_StatusCode",  &m_mctruth_particles_status_code, "MCTruth_particles_StatusCode[MCTruth_NParticles]/I");
    eventweight_tree->Branch("MCTruth_particles_NumberDaughters",  &m_mctruth_particles_num_daughters ,"MCTruth_particles_NumberDaughters[MCTruth_NParticles]/I" );
    eventweight_tree->Branch("MCTruth_particles_Daughters",  &m_mctruth_particles_daughters,"MCTruth_particles_Daughters[MCTruth_NParticles][100]" );
    eventweight_tree->Branch("MCTruth_particles_Gvx",  &m_mctruth_particles_Gvx,"MCTruth_particles_Gvx[MCTruth_NParticles]/D");
    eventweight_tree->Branch("MCTruth_particles_Gvy",  &m_mctruth_particles_Gvy,"MCTruth_particles_Gvy[MCTruth_NParticles]/D" );
    eventweight_tree->Branch("MCTruth_particles_Gvz",  &m_mctruth_particles_Gvz,"MCTruth_particles_Gvz[MCTruth_NParticles]/D" );
    eventweight_tree->Branch("MCTruth_particles_Gvt",  &m_mctruth_particles_Gvt,"MCTruth_particles_Gvt[MCTruth_NParticles]/D" );
    eventweight_tree->Branch("MCTruth_particles_px0",  &m_mctruth_particles_px0, "MCTruth_particles_px0[MCTruth_NParticles]/D"  );
    eventweight_tree->Branch("MCTruth_particles_py0",  &m_mctruth_particles_py0, "MCTruth_particles_py0[MCTruth_NParticles]/D"  );
    eventweight_tree->Branch("MCTruth_particles_pz0",  &m_mctruth_particles_pz0,  "MCTruth_particles_pz0[MCTruth_NParticles]/D"  );
    eventweight_tree->Branch("MCTruth_particles_e0",  &m_mctruth_particles_e0, "MCTruth_particles_e0[MCTruth_NParticles]/D" );
    eventweight_tree->Branch("MCTruth_particles_Rescatter",  &m_mctruth_particles_rescatter,"MCTruth_particles_Rescatter[MCTruth_NParticles]/I" );
    eventweight_tree->Branch("MCTruth_particles_polx",  &m_mctruth_particles_polx, "MCTruth_particles_polx[MCTruth_NParticles]/D" );
    eventweight_tree->Branch("MCTruth_particles_poly",  &m_mctruth_particles_poly, "MCTruth_particles_poly[MCTruth_NParticles]/D" );
    eventweight_tree->Branch("MCTruth_particles_polz",  &m_mctruth_particles_polz, "MCTruth_particles_polz[MCTruth_NParticles]/D");
    eventweight_tree->Branch("MCTruth_neutrino_CCNC",  &m_mctruth_neutrino_ccnc );
    eventweight_tree->Branch("MCTruth_neutrino_mode",  &m_mctruth_neutrino_mode );
    eventweight_tree->Branch("MCTruth_neutrino_interactionType",  &m_mctruth_neutrino_interaction_type );
    eventweight_tree->Branch("MCTruth_neutrino_target",  &m_mctruth_neutrino_target );
    eventweight_tree->Branch("MCTruth_neutrino_nucleon",  &m_mctruth_neutrino_nucleon );
    eventweight_tree->Branch("MCTruth_neutrino_quark",  &m_mctruth_neutrino_quark );
    eventweight_tree->Branch("MCTruth_neutrino_W",  &m_mctruth_neutrino_w );
    eventweight_tree->Branch("MCTruth_neutrino_X",  &m_mctruth_neutrino_x );
    eventweight_tree->Branch("MCTruth_neutrino_Y",  &m_mctruth_neutrino_y );
    eventweight_tree->Branch("MCTruth_neutrino_QSqr",  &m_mctruth_neutrino_qsqr );

     //---------------------gtruth
    eventweight_tree->Branch("GTruth_IsSeaQuark",  &m_gtruth_is_sea_quark );
    eventweight_tree->Branch("GTruth_tgtPDG",  &m_gtruth_tgt_pdg );
    eventweight_tree->Branch("GTruth_tgtA",  &m_gtruth_tgt_A );
    eventweight_tree->Branch("GTruth_tgtZ",  &m_gtruth_tgt_Z );
    eventweight_tree->Branch("GTruth_TgtP4x",  &m_gtruth_tgt_p4_x );
    eventweight_tree->Branch("GTruth_TgtP4y",  &m_gtruth_tgt_p4_y );
    eventweight_tree->Branch("GTruth_TgtP4z",  &m_gtruth_tgt_p4_z );
    eventweight_tree->Branch("GTruth_TgtP4E",  &m_gtruth_tgt_p4_E );
    eventweight_tree->Branch("GTruth_weight",  &m_gtruth_weight );
    eventweight_tree->Branch("GTruth_probability",  &m_gtruth_probability );
    eventweight_tree->Branch("GTruth_Xsec",  &m_gtruth_xsec );
    eventweight_tree->Branch("GTruth_DiffXsec", &m_gtruth_diff_xsec );
    eventweight_tree->Branch("GTruth_GPhaseSpace", &m_gtruth_gphase_space );
    eventweight_tree->Branch("GTruth_vertexX",  &m_gtruth_vertex_x );
    eventweight_tree->Branch("GTruth_vertexY",  &m_gtruth_vertex_y );
    eventweight_tree->Branch("GTruth_vertexZ",  &m_gtruth_vertex_z );
    eventweight_tree->Branch("GTruth_vertexT",  &m_gtruth_vertex_T );
    eventweight_tree->Branch("GTruth_Gscatter", &m_gtruth_gscatter );
    eventweight_tree->Branch("GTruth_Gint",  &m_gtruth_gint );
    eventweight_tree->Branch("GTruth_ResNum", &m_gtruth_res_num); 
    eventweight_tree->Branch("GTruth_NumPiPlus",  &m_gtruth_num_piplus); 
    eventweight_tree->Branch("GTruth_NumPi0",  &m_gtruth_num_pi0);
    eventweight_tree->Branch("GTruth_NumPiMinus",  &m_gtruth_num_piminus);  
    eventweight_tree->Branch("GTruth_NumProton",  &m_gtruth_num_proton );
    eventweight_tree->Branch("GTruth_NumNeutron", &m_gtruth_num_neutron );
    eventweight_tree->Branch("GTruth_IsCharm",  &m_gtruth_is_charm );
    eventweight_tree->Branch("GTruth_IsStrange",  &m_gtruth_is_strange );
    eventweight_tree->Branch("GTruth_StrangeHadronPDG",  &m_gtruth_strange_hadron_pdg );
    eventweight_tree->Branch("GTruth_CharmHadronPDG",  &m_gtruth_charm_hadron_pdg );
    eventweight_tree->Branch("GTruth_DecayMode",&m_gtruth_decay_mode);
    eventweight_tree->Branch("GTruth_gX",  &m_gtruth_gx );
    eventweight_tree->Branch("GTruth_gY", &m_gtruth_gy );
    eventweight_tree->Branch("GTruth_gT", &m_gtruth_gt );
    eventweight_tree->Branch("GTruth_gW", &m_gtruth_gw );
    eventweight_tree->Branch("GTruth_gQ2", &m_gtruth_gQ2 );
    eventweight_tree->Branch("GTruth_gq2",  &m_gtruth_gq2 );
    eventweight_tree->Branch("GTruth_ProbePDG",  &m_gtruth_probe_pdg );
    eventweight_tree->Branch("GTruth_ProbeP4x",  &m_gtruth_probe_p4_x );
    eventweight_tree->Branch("GTruth_ProbeP4y",  &m_gtruth_probe_p4_y );
    eventweight_tree->Branch("GTruth_ProbeP4z",  &m_gtruth_probe_p4_z );
    eventweight_tree->Branch("GTruth_ProbeP4E",  &m_gtruth_probe_p4_E );
    eventweight_tree->Branch("GTruth_HitNucP4x", &m_gtruth_hit_nuc_p4_x );
    eventweight_tree->Branch("GTruth_HitNucP4y", &m_gtruth_hit_nuc_p4_y );
    eventweight_tree->Branch("GTruth_HitNucP4z", &m_gtruth_hit_nuc_p4_z );
    eventweight_tree->Branch("GTruth_HitNucP4E", &m_gtruth_hit_nuc_p4_E );
    eventweight_tree->Branch("GTruth_HitNucPos", &m_gtruth_hit_nuc_pos );
    eventweight_tree->Branch("GTruth_FShadSystP4x", &m_gtruth_fs_had_syst_p4_x );
    eventweight_tree->Branch("GTruth_FShadSystP4y", &m_gtruth_fs_had_syst_p4_y );
    eventweight_tree->Branch("GTruth_FShadSystP4z",  &m_gtruth_fs_had_syst_p4_z );
    eventweight_tree->Branch("GTruth_FShadSystP4E",  &m_gtruth_fs_had_syst_p4_E );
    std::cout<<"SinglePhoton:analyze_Eventweigh:eventweight_tree make branches end"<<std::endl;
  }

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
