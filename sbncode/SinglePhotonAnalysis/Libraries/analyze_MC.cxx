#include "sbncode/SinglePhotonAnalysis/Libraries/analyze_MC.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/variables.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/Processors.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/fiducial_volume.h"
#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_gadget.h"

namespace single_photon
{

  //analyze_Geant4.h
  void AnalyzeGeant4( const    std::vector<art::Ptr<simb::MCParticle>> &mcParticleVector, var_all& vars){


    std::vector<int> spacers = Printer_header({"#MCP","  pdg", " Status"," trkID"," Mother"," Process", "      Process_End","   Energy", "    Vertex(x,  ","     y,     ","       z  )"});
    for(size_t j=0;j< mcParticleVector.size();j++){

      const art::Ptr<simb::MCParticle> mcp = mcParticleVector[j];
      //            std::cout<<"PARG: "<<j<<" PDG "<<mcp->PdgCode()<<" Status "<<mcp->StatusCode()<<" trackid: "<<mcp->TrackId()<<" Mothe "<<mcp->Mother()<<" Process "<<mcp->Process()<<" EndProcess "<<mcp->EndProcess()<<" Energy "<<mcp->E()<<" start ("<<mcp->Vx()<<","<<mcp->Vy()<<","<<mcp->Vz()<<")"<<std::endl;

      Printer_content({
          std::to_string(j),
          std::to_string(mcp->PdgCode()),
          std::to_string(mcp->StatusCode()),
          std::to_string(mcp->TrackId()),
          std::to_string(mcp->Mother()),
          mcp->Process(),
          mcp->EndProcess(),
          std::to_string(mcp->E()),
          std::to_string(mcp->Vx()),
          std::to_string(mcp->Vy()),
          std::to_string(mcp->Vz())
          },
          spacers);

      vars.m_geant4_pdg.push_back(mcp->PdgCode());
      vars.m_geant4_trackid.push_back(mcp->TrackId());
      vars.m_geant4_statuscode.push_back(mcp->StatusCode());
      vars.m_geant4_mother.push_back(mcp->Mother());
      vars.m_geant4_E.push_back(mcp->E());
      vars.m_geant4_mass.push_back(mcp->Mass());
      vars.m_geant4_px.push_back(mcp->Px());
      vars.m_geant4_py.push_back(mcp->Py());
      vars.m_geant4_pz.push_back(mcp->Pz());
      vars.m_geant4_vx.push_back(mcp->Vx());
      vars.m_geant4_vy.push_back(mcp->Vy());
      vars.m_geant4_vz.push_back(mcp->Vz());
      vars.m_geant4_end_process.push_back(mcp->EndProcess());
      vars.m_geant4_process.push_back(mcp->Process());
      vars.m_geant4_costheta.push_back(vars.m_geant4_pz.back()/sqrt(pow(vars.m_geant4_pz.back(),2)+pow(vars.m_geant4_px.back(),2)+pow(vars.m_geant4_py.back(),2)));
      vars.m_geant4_dx.push_back(mcp->Px()/sqrt(pow(vars.m_geant4_pz.back(),2)+pow(vars.m_geant4_px.back(),2)+pow(vars.m_geant4_py.back(),2)));
      vars.m_geant4_dy.push_back(mcp->Py()/sqrt(pow(vars.m_geant4_pz.back(),2)+pow(vars.m_geant4_px.back(),2)+pow(vars.m_geant4_py.back(),2)));
      vars.m_geant4_dz.push_back(mcp->Pz()/sqrt(pow(vars.m_geant4_pz.back(),2)+pow(vars.m_geant4_px.back(),2)+pow(vars.m_geant4_py.back(),2)));

      if(j>2)break;
    }

  }

  //analyze_EventWeight.h
  void AnalyzeEventWeight(art::Event const & e, var_all& vars){

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

    vars.m_run_number_eventweight = e.run();
    vars.m_subrun_number_eventweight = e.subRun();
    vars.m_event_number_eventweight = e.event();

    // possibly the wrong variables, but let's see for now...
    //vars.m_mcflux_evtno     = mcFlux->fevtno;
    vars.m_mcflux_nu_pos_x    = nu.Vx();
    vars.m_mcflux_nu_pos_y    = nu.Vy();
    vars.m_mcflux_nu_pos_z    = nu.Vz();
    vars.m_mcflux_nu_mom_x    = nu.Px(); 
    vars.m_mcflux_nu_mom_y    = nu.Py(); 
    vars.m_mcflux_nu_mom_z    = nu.Pz(); 
    vars.m_mcflux_nu_mom_E    = nu.E();
    vars.m_mcflux_ntype     = mcFlux->fntype;
    vars.m_mcflux_ptype     = mcFlux->fptype;
    vars.m_mcflux_nimpwt    = mcFlux->fnimpwt;
    vars.m_mcflux_dk2gen    = mcFlux->fdk2gen;
    vars.m_mcflux_nenergyn  = mcFlux->fnenergyn;
    vars.m_mcflux_tpx       = mcFlux->ftpx;
    vars.m_mcflux_tpy       = mcFlux->ftpy;
    vars.m_mcflux_tpz       = mcFlux->ftpz;
    vars.m_mcflux_tptype    = mcFlux->ftptype;
    vars.m_mcflux_vx        = mcFlux->fvx;
    vars.m_mcflux_vy        = mcFlux->fvy;
    vars.m_mcflux_vz        = mcFlux->fvz;

    // loop MCParticle info for vars.m_mctruth object

    vars.m_mctruth_nparticles = mcTruth->NParticles();

    for (int i = 0; i < vars.m_mctruth_nparticles; i++){

      const simb::MCParticle& mcParticle = mcTruth->GetParticle(i);

      vars.m_mctruth_particles_track_Id[i] = mcParticle.TrackId();
      vars.m_mctruth_particles_pdg_code[i] = mcParticle.PdgCode();
      vars.m_mctruth_particles_mother[i]  = mcParticle.Mother();
      vars.m_mctruth_particles_status_code[i] = mcParticle.StatusCode();
      vars.m_mctruth_particles_num_daughters[i] = mcParticle.NumberDaughters();

      for (int j = 0; j < vars.m_mctruth_particles_num_daughters[i]; j++){

        const simb::MCParticle& daughterMcParticle = mcTruth->GetParticle(j);
        vars.m_mctruth_particles_daughters[i][j] = daughterMcParticle.TrackId();

      }

      vars.m_mctruth_particles_Gvx[i] = mcParticle.Gvx();
      vars.m_mctruth_particles_Gvy[i] = mcParticle.Gvy();
      vars.m_mctruth_particles_Gvz[i] = mcParticle.Gvz();
      vars.m_mctruth_particles_Gvt[i] = mcParticle.Gvt();
      vars.m_mctruth_particles_px0[i] = mcParticle.Px(0);
      vars.m_mctruth_particles_py0[i] = mcParticle.Py(0);
      vars.m_mctruth_particles_pz0[i] = mcParticle.Pz(0);
      vars.m_mctruth_particles_e0[i] = mcParticle.E(0);
      vars.m_mctruth_particles_rescatter[i] = mcParticle.Rescatter();
      vars.m_mctruth_particles_polx[i] = mcParticle.Polarization().X();
      vars.m_mctruth_particles_poly[i] = mcParticle.Polarization().Y();
      vars.m_mctruth_particles_polz[i] = mcParticle.Polarization().Z();
    }

    const simb::MCNeutrino& mcNeutrino = mcTruth->GetNeutrino();

    vars.m_mctruth_neutrino_ccnc = mcNeutrino.CCNC();
    vars.m_mctruth_neutrino_mode = mcNeutrino.Mode();
    vars.m_mctruth_neutrino_interaction_type = mcNeutrino.InteractionType();
    vars.m_mctruth_neutrino_target = mcNeutrino.Target();
    vars.m_mctruth_neutrino_nucleon = mcNeutrino.HitNuc();
    vars.m_mctruth_neutrino_quark = mcNeutrino.HitQuark();
    vars.m_mctruth_neutrino_w = mcNeutrino.W();
    vars.m_mctruth_neutrino_x = mcNeutrino.X();
    vars.m_mctruth_neutrino_y = mcNeutrino.Y();
    vars.m_mctruth_neutrino_qsqr = mcNeutrino.QSqr();

    vars.m_gtruth_is_sea_quark = gTruth->fIsSeaQuark;
    vars.m_gtruth_tgt_pdg = gTruth->ftgtPDG;
    vars.m_gtruth_tgt_A = gTruth->ftgtA;
    vars.m_gtruth_tgt_Z = gTruth->ftgtZ;
    vars.m_gtruth_tgt_p4_x = gTruth->fTgtP4.X();
    vars.m_gtruth_tgt_p4_y = gTruth->fTgtP4.Y();
    vars.m_gtruth_tgt_p4_z = gTruth->fTgtP4.Z();
    vars.m_gtruth_tgt_p4_E = gTruth->fTgtP4.E();

    vars.m_gtruth_weight = gTruth->fweight;
    vars.m_gtruth_probability = gTruth->fprobability;
    vars.m_gtruth_xsec = gTruth->fXsec;
    vars.m_gtruth_diff_xsec = gTruth->fDiffXsec;
    vars.m_gtruth_gphase_space = gTruth->fGPhaseSpace;

    vars.m_gtruth_vertex_x = gTruth->fVertex.X();
    vars.m_gtruth_vertex_y = gTruth->fVertex.Y();
    vars.m_gtruth_vertex_z = gTruth->fVertex.Z();
    vars.m_gtruth_vertex_T = gTruth->fVertex.T();
    vars.m_gtruth_gscatter = gTruth->fGscatter;
    vars.m_gtruth_gint = gTruth->fGint;
    vars.m_gtruth_res_num = gTruth->fResNum;
    vars.m_gtruth_num_piplus = gTruth->fNumPiPlus;
    vars.m_gtruth_num_pi0 = gTruth->fNumPi0;
    vars.m_gtruth_num_piminus = gTruth->fNumPiMinus;
    vars.m_gtruth_num_proton = gTruth->fNumProton;
    vars.m_gtruth_num_neutron = gTruth->fNumNeutron;
    vars.m_gtruth_is_charm = gTruth->fIsCharm;
    vars.m_gtruth_is_strange = gTruth->fIsStrange;
    vars.m_gtruth_decay_mode = gTruth->fDecayMode;
    vars.m_gtruth_strange_hadron_pdg = gTruth->fStrangeHadronPdg;
    vars.m_gtruth_charm_hadron_pdg = gTruth->fCharmHadronPdg;
    vars.m_gtruth_gx = gTruth->fgX;
    vars.m_gtruth_gy = gTruth->fgY;
    vars.m_gtruth_gt = gTruth->fgT;
    vars.m_gtruth_gw = gTruth->fgW;
    vars.m_gtruth_gQ2 = gTruth->fgQ2;
    vars.m_gtruth_gq2 = gTruth->fgq2;
    vars.m_gtruth_probe_pdg = gTruth->fProbePDG;
    vars.m_gtruth_probe_p4_x = gTruth->fProbeP4.X();
    vars.m_gtruth_probe_p4_y = gTruth->fProbeP4.Y();
    vars.m_gtruth_probe_p4_z = gTruth->fProbeP4.Z();
    vars.m_gtruth_probe_p4_E = gTruth->fProbeP4.E();
    vars.m_gtruth_hit_nuc_p4_x = gTruth->fHitNucP4.X();
    vars.m_gtruth_hit_nuc_p4_y = gTruth->fHitNucP4.Y();
    vars.m_gtruth_hit_nuc_p4_z = gTruth->fHitNucP4.Z();
    vars.m_gtruth_hit_nuc_p4_E = gTruth->fHitNucP4.E();
    vars.m_gtruth_hit_nuc_pos = gTruth->fHitNucPos;
    vars.m_gtruth_fs_had_syst_p4_x = gTruth->fFShadSystP4.X();
    vars.m_gtruth_fs_had_syst_p4_y = gTruth->fFShadSystP4.Y();
    vars.m_gtruth_fs_had_syst_p4_z = gTruth->fFShadSystP4.Z();
    vars.m_gtruth_fs_had_syst_p4_E = gTruth->fFShadSystP4.E();

    //moved to inside singlphoontmodule.cc for filter reasons
    //eventweight_tree->Fill();
    std::cout<<"SinglePhoton::AnalyzeEventWeight-eventweight_tree filled"<<std::endl;
  }


  //analyze_MCTruth.h
  //if 1g1p
  //is there at least 1 reco track and shower
  //if there yes, are they in the same slice
  //if yes, what is the slice score?
  //if yes, was it a cosmic slice?
  //if yes is that the neutrino slice?
  //if they are in different slices
  //was the shower/track:
  //in the neutrino slice?
  //in a cosmic slice?

  //if 1g0p 
  //if there is 1 shower
  //what is the slice id? slice score? cosmic?
  //if missing, not-recoed         

  //if missing atleast 1shower, not recoed


  //for a given signal def (`ncdelta` for now) , finds the MCParticles in event
  //loops over association between reco tracks/showers to get associated slice(s)
  //can also look at things like shower energy, conversion length, etc.
  void AnalyzeRecoMCSlices(std::string signal_def, 
      std::vector<PandoraPFParticle> all_PPFPs,
      std::map<int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIDMap,
      std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,
      std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > &trackToMCParticleMap,
      var_all& vars,
      para_all& paras){

    for(size_t index=0; index< all_PPFPs.size(); ++index){
		PandoraPFParticle* temp_ppfp = &all_PPFPs[index];
		if(!temp_ppfp->get_IsNuSlice()) continue;
		vars.m_reco_slice_num_pfps[temp_ppfp->get_SliceID()]++;
      vars.m_reco_slice_num_showers[temp_ppfp->get_SliceID()]+=temp_ppfp->get_HasShower();
      vars.m_reco_slice_num_tracks [temp_ppfp->get_SliceID()]+=temp_ppfp->get_HasTrack();
    }

    //first check if in the event there's a match to a given signal
    if(signal_def == "ncdelta"){
      //@para updated in the AnalyzeMCTruths function @ analyze_MCTruth.h
      std::cout<<"AnalyzeSlice()\t||\t looking for signal def "<<signal_def<<", vars.m_mctruth_is_delta_radiative = "<<vars.m_mctruth_is_delta_radiative<<std::endl; 

      std::vector<int> matched_shower_ids;
      //first look for sim showers
      for (unsigned int j = 0; j< vars.m_sim_shower_parent_pdg.size(); j++){
        int parent= vars.m_sim_shower_parent_pdg[j];
        int pdg =   vars.m_sim_shower_pdg[j];

        //if this sim shower is a photon and it's primary (parent pdg is -1)
        if(parent == -1 && pdg ==22){
          //first check that this particle isn't alread saved
          //use map from track ID to get MCP
          //if this shower is matched to a recob:shower
          if (vars.m_sim_shower_matched[j] > 0 &&  vars.m_reco_shower_energy_max[j] >20){
            int matched_shower_id = vars.m_sim_shower_trackID[j];


            //if this shower isn't already stored
            if (std::find(matched_shower_ids.begin(), matched_shower_ids.end(), matched_shower_id) == matched_shower_ids.end()){
              matched_shower_ids.push_back(matched_shower_id);
              vars.m_matched_signal_shower_overlay_fraction.push_back(vars.m_sim_shower_overlay_fraction[j]);
              //vars.m_matched_signal_shower_conversion_length;
              vars.m_matched_signal_shower_true_E.push_back(vars.m_sim_shower_energy[j]);
              vars.m_matched_signal_shower_nuscore.push_back( vars.m_reco_shower_nuscore[j]);
              int id = vars.m_reco_shower_sliceId[j];
              std::cout<<"found matched photon shower in slice "<<id<<" with vars.m_sim_shower_energy[j] = "<<vars.m_sim_shower_energy[j]<<std::endl;

              vars.m_matched_signal_shower_sliceId.push_back(id);


              vars.m_matched_signal_shower_is_clearcosmic.push_back( vars.m_reco_shower_isclearcosmic[j]);
              vars.m_matched_signal_shower_is_nuslice.push_back(vars.m_reco_shower_is_nuslice[j]);
              // num track/shower in slice here is in reverse order
              vars.m_matched_signal_shower_tracks_in_slice.push_back( vars.m_reco_slice_num_showers[id]);
              vars.m_matched_signal_shower_showers_in_slice.push_back(vars.m_reco_slice_num_tracks[id]);
              // std::cout<<"found signal photon shower pdg"<< vars.m_sim_shower_pdg[j]<<"and is in neutrino slice =  "<< vars.m_sim_shower_is_nuslice[j]<<std::endl;

            }//if not already stored
          }//if matched to a reco shower >20MeV
        }//if it's a photon from the neutrino interaction
      }//for all sim showers

      vars.m_matched_signal_shower_num = vars.m_matched_signal_shower_true_E.size();

      //NEXT, same procedure for tracks
      std::vector<int> matched_track_ids;
      for (unsigned int k = 0; k< vars.m_sim_track_parent_pdg.size(); k++){
        int parent= vars.m_sim_track_parent_pdg[k];
        int pdg =  vars.m_sim_track_pdg[k];

        int matched_track_id = vars.m_sim_track_trackID[k];

        //if this sim track is a proton and it's primary (parent pdg is -1)
        if((parent == -1 ||parent == 12 || parent ==14 ) && pdg == 2212){

          if (vars.m_sim_track_matched[k] > 0){

            if (std::find(matched_track_ids.begin(), matched_track_ids.end(), matched_track_id) == matched_track_ids.end()){
              matched_track_ids.push_back(matched_track_id);


              // vars.m_matched_signal_track_overlay_fraction.push_back(vars.m_sim_track_overlay_fraction[j]);
              vars.m_matched_signal_track_true_E.push_back(vars.m_sim_track_energy[k]);
              vars.m_matched_signal_track_nuscore.push_back( vars.m_reco_track_nuscore[k]);
              vars.m_matched_signal_track_sliceId.push_back(vars.m_reco_track_sliceId[k]);
              vars.m_matched_signal_track_is_clearcosmic.push_back( vars.m_reco_track_isclearcosmic[k]);
              vars.m_matched_signal_track_is_nuslice.push_back(vars.m_reco_track_is_nuslice[k]);

              int id = vars.m_reco_track_sliceId[k];
              vars.m_matched_signal_track_tracks_in_slice.push_back(vars.m_reco_slice_num_tracks[ id]);
              vars.m_matched_signal_track_showers_in_slice.push_back(vars.m_reco_slice_num_showers[ id]);

            }

          }//if matched
        }//if proton from neutrino interaction
      }//for all sim tracks
      vars.m_matched_signal_track_num = vars.m_matched_signal_track_true_E.size();
    }//end of "ncdelta" scenario

    //brief summary
    if (vars.m_matched_signal_shower_num > 1) vars.m_multiple_matched_showers = true;
    if (vars.m_matched_signal_track_num > 1) vars.m_multiple_matched_tracks = true;
    if (vars.m_matched_signal_shower_num == 0)  vars.m_no_matched_showers = true;

    //check if either 1g1p or 1g0p topology
    if (vars.m_matched_signal_shower_num ==1  && vars.m_matched_signal_track_num ==1){//1g1p
      //check if same slice
      vars.m_is_matched_1g1p = true;
      if ( vars.m_matched_signal_track_sliceId[0] == vars.m_matched_signal_shower_sliceId[0]){
        vars.m_reco_1g1p_is_same_slice = true;
        vars.m_reco_1g1p_is_nuslice = vars.m_matched_signal_shower_is_nuslice[0];
        vars.m_reco_1g1p_nuscore = vars.m_matched_signal_shower_nuscore[0];
      } else{
        vars.m_reco_1g1p_is_multiple_slices = true;
      }
    }else if(vars.m_matched_signal_shower_num ==1  && vars.m_matched_signal_track_num ==0){//1g0p
      vars.m_reco_1g0p_is_nuslice = vars.m_matched_signal_shower_is_nuslice[0];
      vars.m_reco_1g0p_nuscore =  vars.m_matched_signal_shower_nuscore[0];
      vars.m_is_matched_1g0p = true;

    }
  }//findslice

  //This only look at MCTruch info. Reco matching create sivars.m_shower/track for pairing up MCTruth to Reco objects;
  void AnalyzeMCTruths(std::vector<art::Ptr<simb::MCTruth>> & mcTruthVector , std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector,  var_all& vars, para_all& paras){

    std::map<int,std::string> is_delta_map = {
      {2224,"Delta++"},
      {2214,"Delta+"},
      {1114,"Delta-"},
      {2114,"Delta0"},
      {-2224,"Anti-Delta++"},
      {-2214,"Anti-Delta+"},
      {-1114,"Anti-Delta-"},
      {-2114,"Anti-Delta0"}};

    vars.m_mctruth_num = mcTruthVector.size();
    if(g_is_verbose) std::cout<<"# of simb::MCTruth: "<<vars.m_mctruth_num<<std::endl;
    if(vars.m_mctruth_num >1){
      std::cout<<"AnalyzeMCTruths()\t||\t WARNING There is more than 1 MCTruth neutrino interaction. Just running over the first simb::MCTruth."<<std::endl;
    }else if(vars.m_mctruth_num==0){
      std::cout<<"AnalyzeMCTruths()\t||\t WARNING There is 0 MCTruth neutrino interaction. Break simb::MCTruth."<<std::endl;
    }

    //one mctruth per event.  contains list of all particles 

    std::cout<<std::endl;
    std::vector<int> spacers = Printer_header({" NuPdg"," CC=0"," TruthVertex(x,","    y,      ",",       z  )"});
    for(int i=0; i<std::min(1,vars.m_mctruth_num); i++){
      const art::Ptr<simb::MCTruth> truth = mcTruthVector[i];


      vars.m_mctruth_origin = truth->Origin();
      //            if(g_is_verbose) std::cout<<"Getting origin "<<truth->Origin()<<std::endl;

      if(!truth->NeutrinoSet()){
        if(g_is_verbose) std::cout<<"Warning, no neutrino set skipping. "<<std::endl;
      }else{
        //                if(g_is_verbose) std::cout<<"Getting origin "<<truth->Origin()<<std::endl;
        vars.m_mctruth_ccnc = truth->GetNeutrino().CCNC();
        //                if(g_is_verbose) std::cout<<"Getting ccnc "<<truth->GetNeutrino().CCNC()<<std::endl;
        vars.m_mctruth_mode = truth->GetNeutrino().Mode();
        //                if(g_is_verbose) std::cout<<"Getting Mode"<<std::endl;
        vars.m_mctruth_interaction_type = truth->GetNeutrino().InteractionType();
        //                if(g_is_verbose) std::cout<<"Getting Type"<<std::endl;
        vars.m_mctruth_qsqr = truth->GetNeutrino().QSqr();
        //                if(g_is_verbose) std::cout<<"Getting Q"<<std::endl;
        vars.m_mctruth_nu_pdg = truth->GetNeutrino().Nu().PdgCode();
        //                if(g_is_verbose) std::cout<<"Getting E"<<std::endl;
        vars.m_mctruth_nu_E = truth->GetNeutrino().Nu().E();
        //                if(g_is_verbose) std::cout<<"Getting pdg"<<std::endl;
        vars.m_mctruth_lepton_pdg = truth->GetNeutrino().Lepton().PdgCode();
        //                if(g_is_verbose) std::cout<<"Getting pdg lepton"<<std::endl;
        vars.m_mctruth_lepton_E = truth->GetNeutrino().Lepton().E();
        //                if(g_is_verbose) std::cout<<"Getting lepton E"<<std::endl;

        //                if(g_is_verbose) std::cout<<"Getting SC corrected vertex position"<<std::endl;
        std::vector<double> corrected(3);
        // get corrected lepton position
        spacecharge_correction( truth->GetNeutrino().Lepton(),corrected);

        vars.m_mctruth_nu_vertex_x = corrected[0];
        vars.m_mctruth_nu_vertex_y = corrected[1];
        vars.m_mctruth_nu_vertex_z = corrected[2];
        vars.m_mctruth_reco_vertex_dist = sqrt(pow (vars.m_mctruth_nu_vertex_x-vars.m_vertex_pos_x,2)+pow (vars.m_mctruth_nu_vertex_y-vars.m_vertex_pos_y,2)+pow (vars.m_mctruth_nu_vertex_z-vars.m_vertex_pos_z,2));

        //std::vector<int> spacers = Printer_header({"NuPdg","CC=0","TruthVertex(x,","   y,      ",",      z  )"});
        Printer_content(
            {std::to_string(vars.m_mctruth_nu_pdg),
            std::to_string(vars.m_mctruth_ccnc),
            std::to_string(corrected[0]),
            std::to_string(corrected[1]),
            std::to_string(corrected[2])
            },spacers);

      }



      vars.m_mctruth_num_daughter_particles = truth->NParticles(); //MCTruth_NParticles


      if(g_is_verbose) std::cout<<"\nThis MCTruth has "<<truth->NParticles()<<" daughters "<<std::endl;
      std::vector<int> spacers = Printer_header({"           pdg"," TrkID"," MotherID","         Status","        Energy"});




      //some temp variables to see if its 1g1p or 1g1n
      int tmp_n_photons_from_delta = 0; 
      int tmp_n_protons_from_delta = 0; 
      int tmp_n_neutrons_from_delta = 0; 


      vars.m_mctruth_leading_exiting_proton_energy = -9999;

      for(int j=0; j< vars.m_mctruth_num_daughter_particles; j++){

        const simb::MCParticle par = truth->GetParticle(j);
        vars.m_mctruth_daughters_pdg[j] = par.PdgCode();
        vars.m_mctruth_daughters_E[j] = par.E();

        vars.m_mctruth_daughters_status_code[j] = par.StatusCode();
        vars.m_mctruth_daughters_trackID[j] = par.TrackId();
        vars.m_mctruth_daughters_mother_trackID[j] = par.Mother();
        vars.m_mctruth_daughters_px[j] = par.Px();
        vars.m_mctruth_daughters_py[j] = par.Py();
        vars.m_mctruth_daughters_pz[j] = par.Pz();
        vars.m_mctruth_daughters_startx[j] = par.Vx();
        vars.m_mctruth_daughters_starty[j] = par.Vy();
        vars.m_mctruth_daughters_startz[j] = par.Vz();
        vars.m_mctruth_daughters_time[j] = par.T();
        vars.m_mctruth_daughters_endx[j] = par.EndX();
        vars.m_mctruth_daughters_endy[j] = par.EndY();
        vars.m_mctruth_daughters_endz[j] = par.EndZ();
        vars.m_mctruth_daughters_endtime[j] = par.EndT();
        vars.m_mctruth_daughters_process[j] = par.Process();  //Process() and EndProcess() return string
        vars.m_mctruth_daughters_end_process[j] = par.EndProcess();

        if(paras.s_is_textgen) continue; //quick hack, fix in files

        switch(vars.m_mctruth_daughters_pdg[j]){
          case(22): // if it's a gamma
            {
              if(par.StatusCode() == 1){
                vars.m_mctruth_num_exiting_photons++;
                vars.m_mctruth_exiting_photon_mother_trackID.push_back(par.Mother());
                vars.m_mctruth_exiting_photon_trackID.push_back(par.TrackId());
                vars.m_mctruth_exiting_photon_energy.push_back(par.E());
                vars.m_mctruth_exiting_photon_px.push_back(par.Px());
                vars.m_mctruth_exiting_photon_py.push_back(par.Py());
                vars.m_mctruth_exiting_photon_pz.push_back(par.Pz());
              }
              //                         if(g_is_verbose)   std::cout<<"AnalyzeMCTruths()\t||\t Photon "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and photon energy "<<par.E()<<std::endl;

              //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
              if((par.StatusCode()==1 || par.StatusCode()==14 )){
                const  simb::MCParticle mother = truth->GetParticle(par.Mother());

                if(is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                  vars.m_mctruth_delta_photon_energy = par.E();
                  tmp_n_photons_from_delta ++;
                  vars.m_mctruth_is_delta_radiative++;
                }
              }
            }
            break;
          case(111): // if it's a pi0
            {
              // Make sure the pi0 actually exits the nucleus
              if (par.StatusCode() == 1) {
                vars.m_mctruth_exiting_pi0_E.push_back(par.E());
                vars.m_mctruth_exiting_pi0_mom.push_back(sqrt(pow(par.Px(),2)+pow(par.Py(),2)+pow(par.Pz(),2)));
                vars.m_mctruth_exiting_pi0_px.push_back(par.Px());
                vars.m_mctruth_exiting_pi0_py.push_back(par.Py());
                vars.m_mctruth_exiting_pi0_pz.push_back(par.Pz());
                vars.m_mctruth_num_exiting_pi0++;
              }
              break;
            }
          case(211):
          case(-211):  // it's pi+ or pi-
            if (par.StatusCode() == 1) {
              vars.m_mctruth_num_exiting_pipm++;
            }
            break;
          case(2212):  // if it's a proton
            {
              if(par.StatusCode() == 1){
                vars.m_mctruth_num_exiting_protons++;
                vars.m_mctruth_exiting_proton_mother_trackID.push_back(par.Mother());
                vars.m_mctruth_exiting_proton_trackID.push_back(par.TrackId());
                vars.m_mctruth_exiting_proton_energy.push_back(par.E());
                vars.m_mctruth_exiting_proton_px.push_back(par.Px());
                vars.m_mctruth_exiting_proton_py.push_back(par.Py());
                vars.m_mctruth_exiting_proton_pz.push_back(par.Pz());
              }


              //                         if(g_is_verbose)   std::cout<<"SingleProton::AnalyzeMCTruths()\t||\t Proton "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and proton energy "<<par.E()<<std::endl;


              //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
              if(par.StatusCode()==14 ){

                const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                if(is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                  vars.m_mctruth_delta_proton_energy = par.E();
                  tmp_n_protons_from_delta ++;
                }
              }


              break;
            }
          case(2112): // if it's a neutron
            {

              vars.m_mctruth_num_exiting_neutrons++;  // Guanqun: neutron always exits the nucleus? should check it
              vars.m_mctruth_exiting_neutron_mother_trackID.push_back(par.Mother());
              vars.m_mctruth_exiting_neutron_trackID.push_back(par.TrackId());
              vars.m_mctruth_exiting_neutron_energy.push_back(par.E());
              vars.m_mctruth_exiting_neutron_px.push_back(par.Px());
              vars.m_mctruth_exiting_neutron_py.push_back(par.Py());
              vars.m_mctruth_exiting_neutron_pz.push_back(par.Pz());

              //                      if(g_is_verbose)      std::cout<<"SingleProton::AnalyzeMCTruths()\t||\t Neutron "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and neutron energy "<<par.E()<<std::endl;

              //if its mother is a delta with statuscode 3, and it has status code 14, then its the internal product of the delta decay.
              if(par.StatusCode()==14){
                const  simb::MCParticle mother = truth->GetParticle(par.Mother());
                if(is_delta_map.count(mother.PdgCode())>0 && mother.StatusCode()==3){
                  vars.m_mctruth_delta_neutron_energy = par.E();
                  tmp_n_neutrons_from_delta ++;
                }
              }
            }

            break;
          case(-2224):
          case(2224):
            if(par.StatusCode() == 1){  vars.m_mctruth_num_exiting_deltapp++; }
            break;
          case(-2214):
          case(2214)://delta +
          case(-1114):
          case(1114): // if it's delta-
            if(par.StatusCode() == 1){ vars.m_mctruth_num_exiting_deltapm++; }
            break;
          case(-2114):
          case(2114): // if it's delta0
            if(par.StatusCode() == 1){ 
              vars.m_mctruth_num_exiting_delta0++;
              vars.m_mctruth_exiting_delta0_num_daughters.push_back(par.NumberDaughters());
              //                         if(g_is_verbose)   std::cout<<"AnalyzeMCTruths()\t||\t Delta0 "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with "<<vars.m_mctruth_exiting_delta0_num_daughters.back()<<" daughters. StatusCode "<<par.StatusCode()<<std::endl;
            }
            break;
          default:
            break;
        }

        //      if(g_is_verbose)      std::cout<<"SingleProton::AnalyzeMCTruths()\t||\t Neutron "<<par.PdgCode()<<" (id: "<<par.TrackId()<<") with mother trackID: "<<par.Mother()<<". Status Code: "<<par.StatusCode()<<" and neutron energy "<<par.E()<<std::endl;
        Printer_content(
            {std::to_string(par.PdgCode()),
            std::to_string(par.TrackId()),
            std::to_string(par.Mother()),
            std::to_string(par.StatusCode()),
            std::to_string(par.E())
            },spacers);
      } // end of vars.m_mctruth_num_daughter_particles loop

      if(paras.s_is_textgen) continue; //quick hack, fix in files

      for(size_t p=0; p< vars.m_mctruth_exiting_proton_energy.size(); p++){
        if( vars.m_mctruth_exiting_proton_energy[p] > vars.m_mctruth_leading_exiting_proton_energy ){
          vars.m_mctruth_leading_exiting_proton_energy = vars.m_mctruth_exiting_proton_energy[p];
        }
      }



      std::cout<<"AnalyzeMCTruths()\t||\t This event is  ";
      if(tmp_n_photons_from_delta==1 && tmp_n_protons_from_delta==1){
        vars.m_mctruth_delta_radiative_1g1p_or_1g1n = 1;
        std::cout<<"a 1g1p delta radiative event"<<std::endl;
      }else if(tmp_n_photons_from_delta==1 && tmp_n_neutrons_from_delta==1){
        vars.m_mctruth_delta_radiative_1g1p_or_1g1n = 0;
        std::cout<<"a 1g1n delta radiative event"<<std::endl;
      }else{
        std::cout<<"NOT a 1g1p or 1g1n delta radiative decay"<<std::endl;;
      }

      //Now for FSI exiting particles!
      vars.m_mctruth_exiting_photon_from_delta_decay.resize(vars.m_mctruth_num_exiting_photons,0);
      vars.m_mctruth_exiting_proton_from_delta_decay.resize(vars.m_mctruth_num_exiting_protons,0);


      //second loop for some dauhter info
      // status codes!
      // 0 initial state
      // 1 stable final state
      // 2 intermediate state
      // 3 decayed state
      // 11 Nucleon target
      // 14 hadron in the nucleas 

      // So if a  final_state_particle has a status(3) delta in its history its "from" a delta.
      //first we loop over all 14's to see which have a direct mother delta. [done above]
      //so first we loop over all state 1 (exiting) to see what a LArTPC sees (post FSI)
      for (unsigned int p = 0; p <  vars.m_mctruth_exiting_photon_energy.size(); p++){
        // paras.s_exiting_photon_energy_threshold is read from pset
        if ( vars.m_mctruth_exiting_photon_energy[p] > paras.s_exiting_photon_energy_threshold){
          vars.m_mctruth_num_reconstructable_protons++;

        }//if g above threshold
      }

      //if it's a true delta radiative event, check the energies



      if (vars.m_mctruth_is_delta_radiative==true){//if ncdelta
        for (unsigned int p = 0; p <  vars.m_mctruth_exiting_photon_energy.size(); p++){
          std::cout<<"AnalyzeMCTruths()\t||\tLooking at exiting photon with energy "<<vars.m_mctruth_exiting_photon_energy[p]<<std::endl;
          if ( vars.m_mctruth_exiting_photon_energy[p] > paras.s_exiting_photon_energy_threshold){
            vars.m_mctruth_is_reconstructable_1g0p = true; // Guanqun: just means now we have a reconstructable shower, but we don't know if there is a reconstructed nucleon or not yet.

          }//if g above threshold
        }//for all exiting g
        for(unsigned int pr = 0; pr <  vars.m_mctruth_exiting_proton_energy.size(); pr++){
          if ( vars.m_mctruth_exiting_proton_energy[pr]> paras.s_exiting_proton_energy_threshold){
            //if it's already 1g1p then we've found a 1g2p which we aren't counting
            // Guanqun: limit to only 1 reconstructable proton?
            if( vars.m_mctruth_is_reconstructable_1g1p == true && vars.m_mctruth_is_reconstructable_1g0p == false){
              vars.m_mctruth_is_reconstructable_1g1p = false;
            }
            //if there's a photon then it's actually a 1g1p
            if( vars.m_mctruth_is_reconstructable_1g0p == true &&  vars.m_mctruth_is_reconstructable_1g1p == false){
              vars.m_mctruth_is_reconstructable_1g1p = true;
              vars.m_mctruth_is_reconstructable_1g0p = false;
            } 
            std::cout<<"AnalyzeMCTruths()\t||\tChecking proton with energy "<<vars.m_mctruth_exiting_proton_energy[pr]<<", is 1g1p/1g0p= "<< vars.m_mctruth_is_reconstructable_1g1p<<"/"<< vars.m_mctruth_is_reconstructable_1g0p<<std::endl;
          }//if p above threshold
        }//for all exiting p

      }//if ncdelta


      //So for all photons that have status code 1 i.e all exiting ones...
      for(int p =0; p < vars.m_mctruth_num_exiting_photons; ++p){
        const simb::MCParticle mother = truth->GetParticle(vars.m_mctruth_exiting_photon_mother_trackID[p]);

        std::cout<<"AnalyzeMCTruths()\t||\t -- gamma ("<<vars.m_mctruth_exiting_photon_trackID[p]<<") of status_code 1.. "<<std::endl;
        std::cout<<"AnalyzeMCTruths()\t||\t ---- with mother "<<mother.PdgCode()<<" ("<<vars.m_mctruth_exiting_photon_mother_trackID[p]<<") status_code "<<mother.StatusCode()<<std::endl;
        simb::MCParticle nth_mother = mother;
        int n_generation = 2;

        // Guanqun: why not consider its first-generation mother?
        // for a photon exiting nucleus, its first mother is always also a photon (photon exits the nucleus, it becomes another photon..)
        while(nth_mother.StatusCode() != 0 || n_generation < 4){

          if(nth_mother.Mother()<0) break;
          nth_mother = truth->GetParticle(nth_mother.Mother()); 
          std::cout<<"AnalyzeMCTruths()\t||\t ---- and "<<n_generation<<"-mother "<<nth_mother.PdgCode()<<" ("<<nth_mother.TrackId()<<") and status_code "<<nth_mother.StatusCode()<<std::endl;
          if( is_delta_map.count(nth_mother.PdgCode())>0 && nth_mother.StatusCode()==3){
            std::cout<<"AnalyzeMCTruths()\t||\t ------ Defintely From a Delta Decay! : "<<is_delta_map[nth_mother.PdgCode()]<<std::endl;
            vars.m_mctruth_exiting_photon_from_delta_decay[p] = 1;
          }
          n_generation++;
        }
      }

      //So for all protons that have status code 1 i.e all exiting ones...
      for(int p =0; p < vars.m_mctruth_num_exiting_protons; ++p){
        const simb::MCParticle mother = truth->GetParticle(vars.m_mctruth_exiting_proton_mother_trackID[p]);

        if(g_is_verbose){
          std::cout<<"AnalyzeMCTruths()\t||\t -- proton ("<<vars.m_mctruth_exiting_proton_trackID[p]<<") of status_code 1.. "<<std::endl;
          std::cout<<"AnalyzeMCTruths()\t||\t ---- with mother "<<mother.PdgCode()<<" ("<<vars.m_mctruth_exiting_proton_mother_trackID[p]<<") status_code "<<mother.StatusCode()<<std::endl;
        }
        simb::MCParticle nth_mother = mother;
        int n_generation = 2;

        while(nth_mother.StatusCode() != 0 && n_generation < 4){
          if(nth_mother.Mother()<0) break;
          nth_mother = truth->GetParticle(nth_mother.Mother()); 
          if(g_is_verbose) std::cout<<"AnalyzeMCTruths()\t||\t ---- and "<<n_generation<<"-mother "<<nth_mother.PdgCode()<<" ("<<nth_mother.TrackId()<<") and status_code "<<nth_mother.StatusCode()<<std::endl;
          if(is_delta_map.count(nth_mother.PdgCode())>0 && nth_mother.StatusCode()==3){
            if(g_is_verbose) std::cout<<"AnalyzeMCTruths()\t||\t ------ Defintely From a Delta Decay! : "<<is_delta_map[nth_mother.PdgCode()]<<std::endl;
            vars.m_mctruth_exiting_proton_from_delta_decay[p] = 1;
          } 
          n_generation++;
        }


      }


      if(g_is_verbose){
        std::cout<<"AnalyzeMCTruths()\t||\t This is a CCNC: "<<vars.m_mctruth_ccnc<<" event with a nu_pdg: "<<vars.m_mctruth_nu_pdg<<" and "<<vars.m_mctruth_num_daughter_particles<<" exiting particles."<<std::endl;
        std::cout<<"AnalyzeMCTruths()\t||\t With  "<<vars.m_mctruth_num_exiting_pi0<<" Pi0, "<<vars.m_mctruth_num_exiting_pipm<<" Pi+/-, "<<vars.m_mctruth_num_exiting_protons<<" Protons, "<<vars.m_mctruth_num_exiting_neutrons<<" neutrons and "<<vars.m_mctruth_num_exiting_delta0<<" delta0, "<<vars.m_mctruth_num_exiting_deltapm<<" deltapm, "<<vars.m_mctruth_num_exiting_deltapp<<" Deltas++"<<std::endl;
      }

    }// end of MCtruth loo

    //make a stupid temp map
    std::map<size_t,size_t> mymap;
    for(size_t k = 0; k < mcParticleVector.size(); k++){
      const art::Ptr<simb::MCParticle> mcp = mcParticleVector[k];
      mymap[mcp->TrackId()]       = k;
    }


    //Just some VERY hacky pi^0 photon stuff
    int npi0check = 0;
    for(size_t k = 0; k < mcParticleVector.size(); k++){
      const art::Ptr<simb::MCParticle> mcp = mcParticleVector[k];

      if(false) std::cout << k << " Mother:"<< mcp->Mother() << " pdgcode: " << mcp->PdgCode() << " trkid: " << mcp->TrackId() << " statuscode: " << mcp->StatusCode() << std::endl;

      // if it's a pi0, its mother trackID is 0 and it has two daughters
      if(mcp->PdgCode() == 111 && mcp->Mother() == 0 && mcp->NumberDaughters()==2 ){
        npi0check++;
        // get its two daughters
        const art::Ptr<simb::MCParticle> dau1 = mcParticleVector[mymap[mcp->Daughter(0)]];
        const art::Ptr<simb::MCParticle> dau2 = mcParticleVector[mymap[mcp->Daughter(1)]];

        if(false)  std::cout<<"On Dau1: "<<" Mother:"<< dau1->Mother()<<" pdgcode: "<<dau1->PdgCode()<<" trkid: "<<dau1->TrackId()<<" statuscode: "<<dau1->StatusCode()<<std::endl;
        if(false)  std::cout<<"On Dau2: "<<" Mother:"<< dau2->Mother()<<" pdgcode: "<<dau2->PdgCode()<<" trkid: "<<dau2->TrackId()<<" statuscode: "<<dau2->StatusCode()<<std::endl;

        double e1 = dau1->E();
        double e2 = dau2->E();

        std::vector<double> raw_1_End  ={dau1->EndX(), dau1->EndY(), dau1->EndZ()};
        std::vector<double> raw_1_Start  ={dau1->Vx(), dau1->Vy(), dau1->Vz()};

        std::vector<double> raw_2_End  ={dau2->EndX(), dau2->EndY(), dau2->EndZ()};
        std::vector<double> raw_2_Start  ={dau2->Vx(), dau2->Vy(), dau2->Vz()};

        std::vector<double> corrected_1_start(3), corrected_2_start(3);
        std::vector<double> corrected_1_end(3), corrected_2_end(3);

        spacecharge_correction(dau1, corrected_1_start, raw_1_Start);
        spacecharge_correction(dau1, corrected_1_end, raw_1_End);

        spacecharge_correction(dau2, corrected_2_start, raw_2_Start);
        spacecharge_correction(dau2, corrected_2_end, raw_2_End);

        for(int p1=0; p1<dau1->NumberDaughters();p1++){
          auto dd = mcParticleVector[mymap[dau1->Daughter(p1)]];
          std::cout<<"Post1 "<<dd->PdgCode()<<" "<<dd->TrackId()<<" "<<dd->StatusCode()<<" "<<dd->EndProcess()<<std::endl;
        }

        for(int p1=0; p1<dau2->NumberDaughters();p1++){
          auto dd = mcParticleVector[mymap[dau2->Daughter(p1)]];
          std::cout<<"Post2 "<<dd->PdgCode()<<" "<<dd->TrackId()<<" "<<dd->StatusCode()<<" "<<dd->EndProcess()<<" "<<dd->E()<<std::endl;
        }

        int exit1 = isInTPCActive(corrected_1_end, paras);
        int exit2 = isInTPCActive(corrected_2_end, paras);

        if(e2<e1){
          vars.m_mctruth_pi0_leading_photon_energy = e1;
          vars.m_mctruth_pi0_subleading_photon_energy = e2;
          vars.m_mctruth_pi0_leading_photon_end_process = dau1->EndProcess();
          vars.m_mctruth_pi0_subleading_photon_end_process = dau2->EndProcess();
          vars.m_mctruth_pi0_leading_photon_start = corrected_1_start;
          vars.m_mctruth_pi0_leading_photon_end = corrected_1_end;
          vars.m_mctruth_pi0_subleading_photon_start = corrected_2_start;
          vars.m_mctruth_pi0_subleading_photon_end = corrected_2_end;
          //note: the order of subleading/leading photon is reversed// Fixed as of 2022 reprocess!
          vars.m_mctruth_pi0_leading_photon_exiting_TPC =exit1; 
          vars.m_mctruth_pi0_subleading_photon_exiting_TPC = exit2;
          vars.m_mctruth_pi0_leading_photon_mom = {dau1->Px(),dau1->Py(),dau1->Pz()};
          vars.m_mctruth_pi0_subleading_photon_mom = {dau2->Px(),dau2->Py(),dau2->Pz()};

        }else{
          vars.m_mctruth_pi0_leading_photon_energy = e2;
          vars.m_mctruth_pi0_subleading_photon_energy = e1;
          vars.m_mctruth_pi0_leading_photon_end_process = dau2->EndProcess();
          vars.m_mctruth_pi0_subleading_photon_end_process = dau1->EndProcess();
          vars.m_mctruth_pi0_leading_photon_start = corrected_2_start;
          vars.m_mctruth_pi0_leading_photon_end = corrected_2_end;
          vars.m_mctruth_pi0_subleading_photon_start = corrected_1_start;
          vars.m_mctruth_pi0_subleading_photon_end = corrected_1_end;
          vars.m_mctruth_pi0_leading_photon_exiting_TPC = exit2;
          vars.m_mctruth_pi0_subleading_photon_exiting_TPC = exit1;
          vars.m_mctruth_pi0_subleading_photon_mom = {dau1->Px(),dau1->Py(),dau1->Pz()};
          vars.m_mctruth_pi0_leading_photon_mom = {dau2->Px(),dau2->Py(),dau2->Pz()};

        }

      }
    }

    if(npi0check>1) std::cout<<"WARNING WARNING!!!! there are "<<npi0check<<" Pi0's in this event in geant4 that come from the nucleas"<<std::endl;

  }//end of analyze this


}
