#include "sbncode/SinglePhotonAnalysis/Libraries/init_branches.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/fiducial_volume.h"

namespace single_photon
{
  //Process of initialize branches:
  //ClearBranches, ResizeBranches, CreateBranches

  void ClearMeta(){
    //------------ Event related Variables -------------
    m_event_number = -99;
    m_subrun_number = -99;
    m_run_number = -99;
    m_test_matched_hits = 0;

    m_pot_per_event = 0;
    m_pot_per_subrun = m_subrun_pot;
    m_number_of_events_in_subrun = 0;

    m_genie_spline_weight = 1.0;

    //------------ Vertex related Variables -------------
    m_reco_vertex_size = 0;
    m_vertex_pos_x=-99999;
    m_vertex_pos_y=-99999;
    m_vertex_pos_z=-99999;
    m_vertex_pos_tick=-9999;
    m_vertex_pos_wire_p0=-9999;
    m_vertex_pos_wire_p1=-9999;
    m_vertex_pos_wire_p2=-9999;
    m_reco_vertex_in_SCB = -9999;
    m_reco_vertex_dist_to_SCB = -9999;
    m_reco_vertex_dist_to_active_TPC= -9999;
    m_reco_vertex_dist_to_CPA= -9999;

//    m_reco_vertex_to_nearest_dead_wire_plane0=-99999;
//    m_reco_vertex_to_nearest_dead_wire_plane1=-99999;
//    m_reco_vertex_to_nearest_dead_wire_plane2=-99999;

    m_reco_slice_objects = 0;
  }


  void CreateMetaBranches(){

    //true_eventweight_tree
    true_eventweight_tree->Branch("mcweight", "std::map<std::string, std::vector<double>>",&fmcweight);

    //run_subrun_tree
    run_subrun_tree->Branch("run",&m_run,"run/I");
    run_subrun_tree->Branch("subrun",&m_subrun,"subrun/I");
    run_subrun_tree->Branch("subrun_pot",&m_subrun_pot,"subrun_pot/D");
    run_subrun_tree->Branch("subrun_counts",&m_subrun_counts,"subrun_counts/I");

    //pot_tree
    pot_tree->Branch("number_of_events",&m_number_of_events,"number_of_events/I");
    pot_tree->Branch("number_of_vertices",&m_number_of_vertices,"number_of_vertices/I");
    pot_tree->Branch("POT",&m_pot_count,"POT/D");

    //vertex_tree -- part of it
    // --------------------- Event Related variables ------------
    vertex_tree->Branch("run_number", &m_run_number, "run_number/I");
    vertex_tree->Branch("subrun_number", &m_subrun_number, "subrun_number/I");
    vertex_tree->Branch("event_number", &m_event_number, "event_number/I");

    vertex_tree->Branch("pot_per_event",&m_pot_per_event,"pot_per_event/D");
    vertex_tree->Branch("pot_per_subrun",&m_pot_per_subrun,"pot_per_subrun/D");
    vertex_tree->Branch("number_of_events_in_subrun",&m_number_of_events_in_subrun,"number_of_events_in_subrun/D");


    vertex_tree->Branch("genie_spline_weight", &m_genie_spline_weight, "genie_spline_weight/D");
    vertex_tree->Branch("genie_CV_tune_weight", &m_genie_CV_tune_weight, "genie_CV_tune_weight/D");

    vertex_tree->Branch("photonu_weight_low", &m_photonu_weight_low, "photonu_weight_low/D");
    vertex_tree->Branch("photonu_weight_high", &m_photonu_weight_high, "photonu_weight_high/D");

    vertex_tree->Branch("test_matched_hits", &m_test_matched_hits, "test_matched_hits/I");

    // --------------------- Vertex Related variables ------------
    vertex_tree->Branch("reco_vertex_size", &m_reco_vertex_size);
    vertex_tree->Branch("reco_vertex_x", &m_vertex_pos_x);
    vertex_tree->Branch("reco_vertex_y", &m_vertex_pos_y);
    vertex_tree->Branch("reco_vertex_z", &m_vertex_pos_z);
    vertex_tree->Branch("reco_vertex_in_SCB", &m_reco_vertex_in_SCB);
    vertex_tree->Branch("reco_vertex_dist_to_SCB",&m_reco_vertex_dist_to_SCB);
    vertex_tree->Branch("reco_vertex_dist_to_active_TPC",&m_reco_vertex_dist_to_active_TPC);
    vertex_tree->Branch("reco_vertex_dist_to_CPA",&m_reco_vertex_dist_to_CPA);
//    vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane0",&m_reco_vertex_to_nearest_dead_wire_plane0);
//    vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane1",&m_reco_vertex_to_nearest_dead_wire_plane1);
//    vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane2",&m_reco_vertex_to_nearest_dead_wire_plane2);

    vertex_tree->Branch("reco_slice_objects", &m_reco_slice_objects, "reco_slice_objects/I");

    vertex_tree->Branch("m_flash_optfltr_pe_beam",&m_flash_optfltr_pe_beam);
    vertex_tree->Branch("m_flash_optfltr_pe_veto",&m_flash_optfltr_pe_veto);
    vertex_tree->Branch("m_flash_optfltr_pe_veto_tot",&m_flash_optfltr_pe_veto_tot);
    vertex_tree->Branch("m_flash_optfltr_pe_beam_tot",&m_flash_optfltr_pe_beam_tot);
  }

  //isolation.h
  void ClearIsolation(){
    m_isolation_min_dist_trk_shr.clear();
    m_isolation_min_dist_trk_unassoc.clear();

    m_isolation_num_shr_hits_win_1cm_trk.clear();
    m_isolation_num_shr_hits_win_2cm_trk.clear();
    m_isolation_num_shr_hits_win_5cm_trk.clear();
    m_isolation_num_shr_hits_win_10cm_trk.clear();

    m_isolation_num_unassoc_hits_win_1cm_trk.clear();
    m_isolation_num_unassoc_hits_win_2cm_trk.clear();
    m_isolation_num_unassoc_hits_win_5cm_trk.clear();
    m_isolation_num_unassoc_hits_win_10cm_trk.clear();

    m_isolation_nearest_shr_hit_to_trk_wire.clear();
    m_isolation_nearest_shr_hit_to_trk_time.clear();

    m_isolation_nearest_unassoc_hit_to_trk_wire.clear();
    m_isolation_nearest_unassoc_hit_to_trk_time.clear();
  }

  void CreateIsolationBranches(){
    vertex_tree->Branch("isolation_min_dist_trk_shr", &m_isolation_min_dist_trk_shr);
    vertex_tree->Branch("isolation_min_dist_trk_unassoc", &m_isolation_min_dist_trk_unassoc);

    vertex_tree->Branch("isolation_num_shr_hits_win_1cm_trk", &m_isolation_num_shr_hits_win_1cm_trk);
    vertex_tree->Branch("isolation_num_shr_hits_win_2cm_trk", &m_isolation_num_shr_hits_win_2cm_trk);
    vertex_tree->Branch("isolation_num_shr_hits_win_5cm_trk", &m_isolation_num_shr_hits_win_5cm_trk);
    vertex_tree->Branch("isolation_num_shr_hits_win_10cm_trk", &m_isolation_num_shr_hits_win_10cm_trk);

    vertex_tree->Branch("isolation_num_unassoc_hits_win_1cm_trk", &m_isolation_num_unassoc_hits_win_1cm_trk);
    vertex_tree->Branch("isolation_num_unassoc_hits_win_2cm_trk", &m_isolation_num_unassoc_hits_win_2cm_trk);
    vertex_tree->Branch("isolation_num_unassoc_hits_win_5cm_trk", &m_isolation_num_unassoc_hits_win_5cm_trk);
    vertex_tree->Branch("isolation_num_unassoc_hits_win_10cm_trk", &m_isolation_num_unassoc_hits_win_10cm_trk);


    vertex_tree->Branch("isolation_nearest_shr_hit_to_trk_wire", &m_isolation_nearest_shr_hit_to_trk_wire);
    vertex_tree->Branch("isolation_nearest_shr_hit_to_trk_time", &m_isolation_nearest_shr_hit_to_trk_time);

    vertex_tree->Branch("isolation_nearest_unassoc_hit_to_trk_wire", &m_isolation_nearest_unassoc_hit_to_trk_wire);
    vertex_tree->Branch("isolation_nearest_unassoc_hit_to_trk_time", &m_isolation_nearest_unassoc_hit_to_trk_time);

  }

  //second_shower_search.h
  void ClearSecondShowers(){
    m_sss_num_unassociated_hits=0;
    m_sss_num_unassociated_hits_below_threshold=0;
    m_sss_num_associated_hits=0;

    m_sss_num_candidates = 0;

    m_sss_candidate_in_nu_slice.clear();
    m_sss_candidate_num_hits.clear();
    m_sss_candidate_num_wires.clear();
    m_sss_candidate_num_ticks.clear();
    m_sss_candidate_plane.clear();
    m_sss_candidate_PCA.clear();
    m_sss_candidate_mean_ADC.clear();
    m_sss_candidate_ADC_RMS.clear();
    m_sss_candidate_impact_parameter.clear();
    m_sss_candidate_fit_slope.clear();
    m_sss_candidate_veto_score.clear();
    m_sss_candidate_fit_constant.clear();
    m_sss_candidate_mean_tick.clear();
    m_sss_candidate_max_tick.clear();
    m_sss_candidate_min_tick.clear();
    m_sss_candidate_min_wire.clear();
    m_sss_candidate_max_wire.clear();
    m_sss_candidate_mean_wire.clear();
    m_sss_candidate_min_dist.clear();
    m_sss_candidate_wire_tick_based_length.clear();
    m_sss_candidate_energy.clear();
    m_sss_candidate_angle_to_shower.clear();
    m_sss_candidate_closest_neighbour.clear();
    m_sss_candidate_matched.clear();
    m_sss_candidate_matched_energy_fraction_best_plane.clear();
    m_sss_candidate_pdg.clear();
    m_sss_candidate_parent_pdg.clear();
    m_sss_candidate_trackid.clear();
    m_sss_candidate_true_energy.clear();
    m_sss_candidate_overlay_fraction.clear();
    m_sss_candidate_remerge.clear();
  }

  void ClearSecondShowers3D(){

    m_sss3d_num_showers = 0;
    m_sss3d_shower_start_x.clear();
    m_sss3d_shower_start_y.clear();
    m_sss3d_shower_start_z.clear();
    m_sss3d_shower_dir_x.clear();
    m_sss3d_shower_dir_y.clear();
    m_sss3d_shower_dir_z.clear();
    m_sss3d_shower_length.clear();
    m_sss3d_shower_conversion_dist.clear();
    m_sss3d_shower_invariant_mass.clear();
    m_sss3d_shower_implied_invariant_mass.clear();
    m_sss3d_shower_impact_parameter.clear();
    m_sss3d_shower_energy_max.clear();
    m_sss3d_shower_score.clear();
    m_sss3d_slice_nu.clear();
    m_sss3d_slice_clear_cosmic.clear();
    m_sss3d_shower_ioc_ratio.clear();
  }

  void ClearStubs(){
    m_trackstub_num_unassociated_hits = 0; 
    m_trackstub_unassociated_hits_below_threshold = 0; 
    m_trackstub_associated_hits=0; 
    m_trackstub_num_candidates=0; 
    m_trackstub_candidate_in_nu_slice.clear();
    m_trackstub_candidate_num_hits.clear();
    m_trackstub_candidate_num_wires.clear(); 
    m_trackstub_candidate_num_ticks.clear();
    m_trackstub_candidate_plane.clear(); 
    m_trackstub_candidate_PCA.clear();
    m_trackstub_candidate_mean_ADC.clear();
    m_trackstub_candidate_ADC_RMS.clear();
    m_trackstub_candidate_veto_score.clear();
    m_trackstub_candidate_mean_tick.clear();
    m_trackstub_candidate_max_tick.clear();
    m_trackstub_candidate_min_tick.clear();
    m_trackstub_candidate_min_wire.clear();
    m_trackstub_candidate_max_wire.clear();
    m_trackstub_candidate_mean_wire.clear();
    m_trackstub_candidate_min_dist.clear();  
    m_trackstub_candidate_min_impact_parameter_to_shower.clear(); 
    m_trackstub_candidate_min_conversion_dist_to_shower_start.clear();  
    m_trackstub_candidate_min_ioc_to_shower_start.clear();        
    m_trackstub_candidate_ioc_based_length.clear();         
    m_trackstub_candidate_wire_tick_based_length.clear();           
    m_trackstub_candidate_mean_ADC_first_half.clear();              
    m_trackstub_candidate_mean_ADC_second_half.clear();
    m_trackstub_candidate_mean_ADC_first_to_second_ratio.clear(); 
    m_trackstub_candidate_track_angle_wrt_shower_direction.clear();   
    m_trackstub_candidate_linear_fit_chi2.clear();          
    m_trackstub_candidate_energy.clear();
    m_trackstub_candidate_remerge.clear(); 
    m_trackstub_candidate_matched.clear(); 
    m_trackstub_candidate_matched_energy_fraction_best_plane.clear(); 
    m_trackstub_candidate_pdg.clear();   
    m_trackstub_candidate_parent_pdg.clear();
    m_trackstub_candidate_trackid.clear(); 
    m_trackstub_candidate_true_energy.clear();
    m_trackstub_candidate_overlay_fraction.clear(); 

    m_trackstub_num_candidate_groups = 0;                
    m_grouped_trackstub_candidate_indices.clear(); 
    m_trackstub_candidate_group_timeoverlap_fraction.clear();   
  }

  void CreateSecondShowerBranches(){
    vertex_tree->Branch("sss_num_unassociated_hits",&m_sss_num_unassociated_hits,"sss_num_unassociated_hits/I");
    vertex_tree->Branch("sss_num_unassociated_hits_below_threshold",&m_sss_num_unassociated_hits_below_threshold,"sss_num_unassociated_hits_below_threshold/I");
    vertex_tree->Branch("sss_num_associated_hits",&m_sss_num_associated_hits,"sss_num_associated_hits/I");

    vertex_tree->Branch("sss_num_candidates",&m_sss_num_candidates,"sss_num_candidates/I");
    vertex_tree->Branch("sss_candidate_veto_score",&m_sss_candidate_veto_score);
    vertex_tree->Branch("sss_candidate_in_nu_slice", &m_sss_candidate_in_nu_slice);
    vertex_tree->Branch("sss_candidate_num_hits",&m_sss_candidate_num_hits);
    vertex_tree->Branch("sss_candidate_num_wires",&m_sss_candidate_num_wires);
    vertex_tree->Branch("sss_candidate_num_ticks",&m_sss_candidate_num_ticks);
    vertex_tree->Branch("sss_candidate_plane",&m_sss_candidate_plane);
    vertex_tree->Branch("sss_candidate_PCA",&m_sss_candidate_PCA);
    vertex_tree->Branch("sss_candidate_mean_ADC",&m_sss_candidate_mean_ADC);
    vertex_tree->Branch("sss_candidate_ADC_RMS", &m_sss_candidate_ADC_RMS);
    vertex_tree->Branch("sss_candidate_impact_parameter",&m_sss_candidate_impact_parameter); 
    vertex_tree->Branch("sss_candidate_fit_slope",&m_sss_candidate_fit_slope);
    vertex_tree->Branch("sss_candidate_fit_constant",&m_sss_candidate_fit_constant);
    vertex_tree->Branch("sss_candidate_mean_tick",&m_sss_candidate_mean_tick);
    vertex_tree->Branch("sss_candidate_max_tick",&m_sss_candidate_max_tick);
    vertex_tree->Branch("sss_candidate_min_tick",&m_sss_candidate_min_tick);
    vertex_tree->Branch("sss_candidate_mean_wire",&m_sss_candidate_mean_wire);
    vertex_tree->Branch("sss_candidate_max_wire",&m_sss_candidate_max_wire);
    vertex_tree->Branch("sss_candidate_min_wire",&m_sss_candidate_min_wire);
    vertex_tree->Branch("sss_candidate_min_dist",&m_sss_candidate_min_dist);
    vertex_tree->Branch("sss_candidate_wire_tick_based_length", &m_sss_candidate_wire_tick_based_length);
    vertex_tree->Branch("sss_candidate_energy",&m_sss_candidate_energy);
    vertex_tree->Branch("sss_candidate_angle_to_shower",&m_sss_candidate_angle_to_shower);
    vertex_tree->Branch("sss_candidate_closest_neighbour",&m_sss_candidate_closest_neighbour);
    vertex_tree->Branch("sss_candidate_remerge",&m_sss_candidate_remerge);

    vertex_tree->Branch("sss_candidate_matched",&m_sss_candidate_matched);
    vertex_tree->Branch("sss_candidate_pdg",&m_sss_candidate_pdg);
    vertex_tree->Branch("sss_candidate_parent_pdg",&m_sss_candidate_parent_pdg);
    vertex_tree->Branch("sss_candidate_trackid",&m_sss_candidate_trackid);
    vertex_tree->Branch("sss_candidate_true_energy", &m_sss_candidate_true_energy);
    vertex_tree->Branch("sss_candidate_overlay_fraction",&m_sss_candidate_overlay_fraction);
    vertex_tree->Branch("sss_candidate_matched_energy_fraction_best_plane", &m_sss_candidate_matched_energy_fraction_best_plane);


    vertex_tree->Branch("sss3d_ioc_ranked_en",&m_sss3d_ioc_ranked_en);
    vertex_tree->Branch("sss3d_ioc_ranked_conv",&m_sss3d_ioc_ranked_conv);
    vertex_tree->Branch("sss3d_ioc_ranked_invar",&m_sss3d_ioc_ranked_invar);
    vertex_tree->Branch("sss3d_ioc_ranked_implied_invar",&m_sss3d_ioc_ranked_implied_invar);
    vertex_tree->Branch("sss3d_ioc_ranked_ioc",&m_sss3d_ioc_ranked_ioc);
    vertex_tree->Branch("sss3d_ioc_ranked_opang",&m_sss3d_ioc_ranked_opang);
    vertex_tree->Branch("sss3d_ioc_ranked_implied_opang",&m_sss3d_ioc_ranked_implied_opang);
    vertex_tree->Branch("sss3d_ioc_ranked_id",&m_sss3d_ioc_ranked_id);

    vertex_tree->Branch("sss3d_invar_ranked_en",&m_sss3d_invar_ranked_en);
    vertex_tree->Branch("sss3d_invar_ranked_conv",&m_sss3d_invar_ranked_conv);
    vertex_tree->Branch("sss3d_invar_ranked_invar",&m_sss3d_invar_ranked_invar);
    vertex_tree->Branch("sss3d_invar_ranked_implied_invar",&m_sss3d_invar_ranked_implied_invar);
    vertex_tree->Branch("sss3d_invar_ranked_ioc",&m_sss3d_invar_ranked_ioc);
    vertex_tree->Branch("sss3d_invar_ranked_opang",&m_sss3d_invar_ranked_opang);
    vertex_tree->Branch("sss3d_invar_ranked_implied_opang",&m_sss3d_invar_ranked_implied_opang);
    vertex_tree->Branch("sss3d_invar_ranked_id",&m_sss3d_invar_ranked_id);


    vertex_tree->Branch("sss2d_ioc_ranked_en",&m_sss2d_ioc_ranked_en);
    vertex_tree->Branch("sss2d_ioc_ranked_conv",&m_sss2d_ioc_ranked_conv);
    vertex_tree->Branch("sss2d_ioc_ranked_ioc",&m_sss2d_ioc_ranked_ioc);
    vertex_tree->Branch("sss2d_ioc_ranked_pca",&m_sss2d_ioc_ranked_pca);
    vertex_tree->Branch("sss2d_ioc_ranked_invar",&m_sss2d_ioc_ranked_invar);
    vertex_tree->Branch("sss2d_ioc_ranked_angle_to_shower",&m_sss2d_ioc_ranked_angle_to_shower);
    vertex_tree->Branch("sss2d_ioc_ranked_num_planes",&m_sss2d_ioc_ranked_num_planes);

    vertex_tree->Branch("sss2d_invar_ranked_en",&m_sss2d_invar_ranked_en);
    vertex_tree->Branch("sss2d_invar_ranked_conv",&m_sss2d_invar_ranked_conv);
    vertex_tree->Branch("sss2d_invar_ranked_ioc",&m_sss2d_invar_ranked_ioc);
    vertex_tree->Branch("sss2d_invar_ranked_pca",&m_sss2d_invar_ranked_pca);
    vertex_tree->Branch("sss2d_invar_ranked_invar",&m_sss2d_invar_ranked_invar);
    vertex_tree->Branch("sss2d_invar_ranked_angle_to_shower",&m_sss2d_invar_ranked_angle_to_shower);
    vertex_tree->Branch("sss2d_invar_ranked_num_planes",&m_sss2d_invar_ranked_num_planes);

    vertex_tree->Branch("sss2d_conv_ranked_en",&m_sss2d_conv_ranked_en);
    vertex_tree->Branch("sss2d_conv_ranked_conv",&m_sss2d_conv_ranked_conv);
    vertex_tree->Branch("sss2d_conv_ranked_ioc",&m_sss2d_conv_ranked_ioc);
    vertex_tree->Branch("sss2d_conv_ranked_pca",&m_sss2d_conv_ranked_pca);
    vertex_tree->Branch("sss2d_conv_ranked_invar",&m_sss2d_conv_ranked_invar);
    vertex_tree->Branch("sss2d_conv_ranked_angle_to_shower",&m_sss2d_conv_ranked_angle_to_shower);
    vertex_tree->Branch("sss2d_conv_ranked_num_planes",&m_sss2d_conv_ranked_num_planes);

  }

  void CreateSecondShowerBranches3D(){
    vertex_tree->Branch("sss3d_num_showers",&m_sss3d_num_showers,"sss3d_num_showers/I");

    vertex_tree->Branch("sss3d_shower_start_x",&m_sss3d_shower_start_x);
    vertex_tree->Branch("sss3d_shower_start_y",&m_sss3d_shower_start_y);
    vertex_tree->Branch("sss3d_shower_start_z",&m_sss3d_shower_start_z);
    vertex_tree->Branch("sss3d_shower_dir_x",&m_sss3d_shower_dir_x);
    vertex_tree->Branch("sss3d_shower_dir_y",&m_sss3d_shower_dir_y);
    vertex_tree->Branch("sss3d_shower_dir_z",&m_sss3d_shower_dir_z);

    vertex_tree->Branch("sss3d_shower_length",&m_sss3d_shower_length);
    vertex_tree->Branch("sss3d_shower_conversion_dist",&m_sss3d_shower_conversion_dist);
    vertex_tree->Branch("sss3d_shower_invariant_mass",&m_sss3d_shower_invariant_mass);
    vertex_tree->Branch("sss3d_shower_implied_invariant_mass",&m_sss3d_shower_implied_invariant_mass);
    vertex_tree->Branch("sss3d_shower_impact_parameter",&m_sss3d_shower_impact_parameter);
    vertex_tree->Branch("sss3d_shower_ioc_ratio",&m_sss3d_shower_ioc_ratio);
    vertex_tree->Branch("sss3d_shower_energy_max",&m_sss3d_shower_energy_max);
    vertex_tree->Branch("sss3d_shower_score",&m_sss3d_shower_score);
    vertex_tree->Branch("sss3d_slice_nu",&m_sss3d_slice_nu);
    vertex_tree->Branch("sss3d_slice_clear_cosmic",&m_sss3d_slice_clear_cosmic);
  }

  void CreateStubBranches(){

    vertex_tree->Branch("trackstub_num_unassociated_hits",&m_trackstub_num_unassociated_hits,"trackstub_num_unassociated_hits/I");
    vertex_tree->Branch("trackstub_unassociated_hits_below_threshold",&m_trackstub_unassociated_hits_below_threshold,"trackstub_unassociated_hits_below_threshold/I");
    vertex_tree->Branch("trackstub_associated_hits",&m_trackstub_associated_hits,"trackstub_associated_hits/I");
    vertex_tree->Branch("trackstub_num_candidates", &m_trackstub_num_candidates, "trackstub_num_candidates/I");
    vertex_tree->Branch("trackstub_candidate_in_nu_slice", &m_trackstub_candidate_in_nu_slice);
    vertex_tree->Branch("trackstub_candidate_num_hits", &m_trackstub_candidate_num_hits);
    vertex_tree->Branch("trackstub_candidate_num_wires", &m_trackstub_candidate_num_wires);
    vertex_tree->Branch("trackstub_candidate_num_ticks", &m_trackstub_candidate_num_ticks);
    vertex_tree->Branch("trackstub_candidate_plane", &m_trackstub_candidate_plane);
    vertex_tree->Branch("trackstub_candidate_PCA", &m_trackstub_candidate_PCA);
    vertex_tree->Branch("trackstub_candidate_mean_ADC", &m_trackstub_candidate_mean_ADC);
    vertex_tree->Branch("trackstub_candidate_ADC_RMS", &m_trackstub_candidate_ADC_RMS);
    vertex_tree->Branch("trackstub_candidate_veto_score", &m_trackstub_candidate_veto_score);
    vertex_tree->Branch("trackstub_candidate_mean_tick", &m_trackstub_candidate_mean_tick);
    vertex_tree->Branch("trackstub_candidate_max_tick", &m_trackstub_candidate_max_tick);
    vertex_tree->Branch("trackstub_candidate_min_tick", &m_trackstub_candidate_min_tick);
    vertex_tree->Branch("trackstub_candidate_min_wire", &m_trackstub_candidate_min_wire);
    vertex_tree->Branch("trackstub_candidate_max_wire", &m_trackstub_candidate_max_wire);
    vertex_tree->Branch("trackstub_candidate_mean_wire", &m_trackstub_candidate_mean_wire);
    vertex_tree->Branch("trackstub_candidate_min_dist", &m_trackstub_candidate_min_dist);
    vertex_tree->Branch("trackstub_candidate_min_impact_parameter_to_shower", &m_trackstub_candidate_min_impact_parameter_to_shower);
    vertex_tree->Branch("trackstub_candidate_min_conversion_dist_to_shower_start", &m_trackstub_candidate_min_conversion_dist_to_shower_start);
    vertex_tree->Branch("trackstub_candidate_min_ioc_to_shower_start", &m_trackstub_candidate_min_ioc_to_shower_start);
    vertex_tree->Branch("trackstub_candidate_ioc_based_length", &m_trackstub_candidate_ioc_based_length);
    vertex_tree->Branch("trackstub_candidate_wire_tick_based_length", &m_trackstub_candidate_wire_tick_based_length);
    vertex_tree->Branch("trackstub_candidate_mean_ADC_first_half", &m_trackstub_candidate_mean_ADC_first_half);
    vertex_tree->Branch("trackstub_candidate_mean_ADC_second_half", &m_trackstub_candidate_mean_ADC_second_half);
    vertex_tree->Branch("trackstub_candidate_mean_ADC_first_to_second_ratio", &m_trackstub_candidate_mean_ADC_first_to_second_ratio);
    vertex_tree->Branch("trackstub_candidate_track_angle_wrt_shower_direction", &m_trackstub_candidate_track_angle_wrt_shower_direction);
    vertex_tree->Branch("trackstub_candidate_linear_fit_chi2", &m_trackstub_candidate_linear_fit_chi2);
    vertex_tree->Branch("trackstub_candidate_energy", &m_trackstub_candidate_energy);
    vertex_tree->Branch("trackstub_candidate_remerge", &m_trackstub_candidate_remerge);
    vertex_tree->Branch("trackstub_candidate_matched", &m_trackstub_candidate_matched);
    vertex_tree->Branch("trackstub_candidate_matched_energy_fraction_best_plane", &m_trackstub_candidate_matched_energy_fraction_best_plane);
    vertex_tree->Branch("trackstub_candidate_pdg", &m_trackstub_candidate_pdg);
    vertex_tree->Branch("trackstub_candidate_parent_pdg", &m_trackstub_candidate_parent_pdg);
    vertex_tree->Branch("trackstub_candidate_trackid", &m_trackstub_candidate_trackid);
    vertex_tree->Branch("trackstub_candidate_true_energy", &m_trackstub_candidate_true_energy);
    vertex_tree->Branch("trackstub_candidate_overlay_fraction", &m_trackstub_candidate_overlay_fraction);


    vertex_tree->Branch("trackstub_num_candidate_groups", &m_trackstub_num_candidate_groups, "trackstub_num_candidate_groups/I");
    vertex_tree->Branch("grouped_trackstub_candidate_indices", &m_grouped_trackstub_candidate_indices);
    vertex_tree->Branch("trackstub_candidate_group_timeoverlap_fraction", &m_trackstub_candidate_group_timeoverlap_fraction);

  }

  void ResizeSecondShowers(size_t size){}

  //analyze_OpFlashes.h
  void ClearFlashes(){
    m_reco_num_flashes =0;
    m_reco_num_flashes_in_beamgate =0;
    m_reco_flash_total_pe.clear();
    m_reco_flash_time.clear();
    m_reco_flash_time_width.clear();
    m_reco_flash_abs_time.clear();
    m_reco_flash_frame.clear();
    m_reco_flash_ycenter.clear();
    m_reco_flash_ywidth.clear();
    m_reco_flash_zcenter.clear();
    m_reco_flash_zwidth.clear();
    m_reco_flash_total_pe_in_beamgate.clear();
    m_reco_flash_time_in_beamgate.clear();
    m_reco_flash_ycenter_in_beamgate.clear();
    m_reco_flash_zcenter_in_beamgate.clear();
    m_CRT_veto_nhits = -999;
    m_CRT_veto_hit_PE.clear();
    m_CRT_min_hit_time = -999;
    m_CRT_min_hit_PE = -999;
    m_CRT_min_hit_x = -999;
    m_CRT_min_hit_y = -999;
    m_CRT_min_hit_z = -999;
    m_CRT_hits_time.clear();
    m_CRT_hits_PE.clear();
    m_CRT_hits_x.clear(); 
    m_CRT_hits_y.clear();
    m_CRT_hits_z.clear();
    m_CRT_dt = -999;

  }

  void CreateFlashBranches(){
    vertex_tree->Branch("beamgate_flash_start",&m_beamgate_flash_start,"beamgate_flash_start/D");
    vertex_tree->Branch("beamgate_flash_end",&m_beamgate_flash_end,"beamgate_flash_end/D");
    vertex_tree->Branch("reco_num_flashes",&m_reco_num_flashes,"reco_num_flashes/I");
    vertex_tree->Branch("reco_num_flashes_in_beamgate",&m_reco_num_flashes_in_beamgate,"reco_num_flashes_in_beamgate/I");
    vertex_tree->Branch("reco_flash_total_pe", &m_reco_flash_total_pe);
    vertex_tree->Branch("reco_flash_time", &m_reco_flash_time);
    vertex_tree->Branch("reco_flash_time_width",&m_reco_flash_time_width);
    vertex_tree->Branch("reco_flash_abs_time",&m_reco_flash_abs_time);
    vertex_tree->Branch("reco_flash_frame",&m_reco_flash_frame);
    vertex_tree->Branch("reco_flash_ycenter",&m_reco_flash_ycenter);
    vertex_tree->Branch("reco_flash_ywidth",&m_reco_flash_ywidth);
    vertex_tree->Branch("reco_flash_zcenter",&m_reco_flash_zcenter);
    vertex_tree->Branch("reco_flash_zwidth",&m_reco_flash_zwidth);
    vertex_tree->Branch("reco_flash_total_pe_in_beamgate", &m_reco_flash_total_pe_in_beamgate);
    vertex_tree->Branch("reco_flash_time_in_beamgate", &m_reco_flash_time_in_beamgate);
    vertex_tree->Branch("reco_flash_ycenter_in_beamgate",&m_reco_flash_ycenter_in_beamgate);
    vertex_tree->Branch("reco_flash_zcenter_in_beamgate",&m_reco_flash_zcenter_in_beamgate);

    vertex_tree->Branch("CRT_veto_nhits",&m_CRT_veto_nhits,"CRT_veto_nhits/I");
    vertex_tree->Branch("CRT_veto_hit_PE",&m_CRT_veto_hit_PE);
    vertex_tree->Branch("CRT_dt",& m_CRT_dt," CRT_dt/D");
    vertex_tree->Branch("CRT_min_hit_time",&m_CRT_min_hit_time,"CRT_min_hit_time/D");
    vertex_tree->Branch("CRT_min_hit_PE",&m_CRT_min_hit_PE,"CRT_min_hit_PE/D");
    vertex_tree->Branch("CRT_min_hit_x",&m_CRT_min_hit_x,"CRT_min_hit_x/D");
    vertex_tree->Branch("CRT_min_hit_y",&m_CRT_min_hit_y,"CRT_min_hit_y/D");
    vertex_tree->Branch("CRT_min_hit_z",&m_CRT_min_hit_z,"CRT_min_hit_z/D");
    vertex_tree->Branch("CRT_hits_time",&m_CRT_hits_time);
    vertex_tree->Branch("CRT_hits_PE",&m_CRT_hits_PE);
    vertex_tree->Branch("CRT_hits_x",&m_CRT_hits_x);
    vertex_tree->Branch("CRT_hits_y",&m_CRT_hits_y);
    vertex_tree->Branch("CRT_hits_z",&m_CRT_hits_z);
  }

  void ResizeFlashes(size_t size){
    m_reco_flash_total_pe.resize(size);
    m_reco_flash_time.resize(size);
    m_reco_flash_time_width.resize(size);
    m_reco_flash_abs_time.resize(size);
    m_reco_flash_frame.resize(size);
    m_reco_flash_ycenter.resize(size);
    m_reco_flash_ywidth.resize(size);
    m_reco_flash_zcenter.resize(size);
    m_reco_flash_zwidth.resize(size);
    m_reco_flash_total_pe_in_beamgate.resize(size);
    m_reco_flash_time_in_beamgate.resize(size);
    m_reco_flash_ycenter_in_beamgate.resize(size);
    m_reco_flash_zcenter_in_beamgate.resize(size);
    m_CRT_veto_hit_PE.resize(size);

    m_CRT_hits_time.resize(size);
    m_CRT_hits_PE.resize(size);
    m_CRT_hits_x.resize(size); 
    m_CRT_hits_y.resize(size);
    m_CRT_hits_z.resize(size);
  }

  //analyze_Tracks.h
  void ClearTracks(){
    m_reco_asso_tracks=0;
    m_reco_track_length.clear();
    m_reco_track_num_daughters.clear();
    m_reco_track_daughter_trackscore.clear();
    m_reco_track_dirx.clear();
    m_reco_track_diry.clear();
    m_reco_track_dirz.clear();
    m_reco_track_startx.clear();
    m_reco_track_starty.clear();
    m_reco_track_startz.clear();
    m_reco_track_endx.clear();
    m_reco_track_endy.clear();
    m_reco_track_endz.clear();
    m_reco_track_end_dist_to_active_TPC.clear();
    m_reco_track_start_dist_to_active_TPC.clear();
    m_reco_track_end_dist_to_CPA.clear();
    m_reco_track_start_dist_to_CPA.clear();
    m_reco_track_end_dist_to_SCB.clear();
    m_reco_track_start_dist_to_SCB.clear();
    m_reco_track_end_in_SCB.clear();
    m_reco_track_start_in_SCB.clear();

    m_reco_track_theta_yz.clear();
    m_reco_track_phi_yx.clear();

    m_reco_track_calo_energy_plane0.clear();
    m_reco_track_calo_energy_plane1.clear();
    m_reco_track_calo_energy_plane2.clear();
    m_reco_track_calo_energy_max.clear();

    m_reco_track_num_trajpoints.clear();
    m_reco_track_num_spacepoints.clear();
    m_reco_track_proton_kinetic_energy.clear();
    m_reco_track_ordered_energy_index.clear();
    m_reco_track_ordered_displacement_index.clear();

    m_reco_track_spacepoint_principal0.clear();
    m_reco_track_spacepoint_principal1.clear();
    m_reco_track_spacepoint_principal2.clear();

    m_reco_track_spacepoint_chi.clear();
    m_reco_track_spacepoint_max_dist.clear();

    m_reco_track_best_calo_plane.clear();

    m_reco_track_mean_dEdx_best_plane.clear();
    m_reco_track_mean_dEdx_start_half_best_plane.clear();
    m_reco_track_mean_dEdx_end_half_best_plane.clear();
    m_reco_track_good_calo_best_plane.clear();
    m_reco_track_trunc_dEdx_best_plane.clear();
    m_reco_track_mean_trunc_dEdx_best_plane.clear();
    m_reco_track_mean_trunc_dEdx_start_half_best_plane.clear();
    m_reco_track_mean_trunc_dEdx_end_half_best_plane.clear();
    m_reco_track_trunc_PIDA_best_plane.clear();
    m_reco_track_resrange_best_plane.clear();
    m_reco_track_dEdx_best_plane.clear();


    m_reco_track_mean_dEdx_p0.clear();
    m_reco_track_mean_dEdx_start_half_p0.clear();
    m_reco_track_mean_dEdx_end_half_p0.clear();
    m_reco_track_good_calo_p0.clear();
    m_reco_track_trunc_dEdx_p0.clear();
    m_reco_track_mean_trunc_dEdx_p0.clear();
    m_reco_track_mean_trunc_dEdx_start_half_p0.clear();
    m_reco_track_mean_trunc_dEdx_end_half_p0.clear();
    m_reco_track_trunc_PIDA_p0.clear();
    m_reco_track_resrange_p0.clear();
    m_reco_track_dEdx_p0.clear();

    m_reco_track_mean_dEdx_p1.clear();
    m_reco_track_mean_dEdx_start_half_p1.clear();
    m_reco_track_mean_dEdx_end_half_p1.clear();
    m_reco_track_good_calo_p1.clear();
    m_reco_track_trunc_dEdx_p1.clear();
    m_reco_track_mean_trunc_dEdx_p1.clear();
    m_reco_track_mean_trunc_dEdx_start_half_p1.clear();
    m_reco_track_mean_trunc_dEdx_end_half_p1.clear();
    m_reco_track_trunc_PIDA_p1.clear();
    m_reco_track_resrange_p1.clear();
    m_reco_track_dEdx_p1.clear();

    m_reco_track_mean_dEdx_p2.clear();
    m_reco_track_mean_dEdx_start_half_p2.clear();
    m_reco_track_mean_dEdx_end_half_p2.clear();
    m_reco_track_good_calo_p2.clear();
    m_reco_track_trunc_dEdx_p2.clear();
    m_reco_track_mean_trunc_dEdx_p2.clear();
    m_reco_track_mean_trunc_dEdx_start_half_p2.clear();
    m_reco_track_mean_trunc_dEdx_end_half_p2.clear();
    m_reco_track_trunc_PIDA_p2.clear();
    m_reco_track_resrange_p2.clear();
    m_reco_track_dEdx_p2.clear();

    m_reco_track_num_calo_hits_p1.clear();
    m_reco_track_num_calo_hits_p0.clear();
    m_reco_track_num_calo_hits_p2.clear();

    m_sim_track_matched.clear();
    m_sim_track_overlay_fraction.clear();
    m_sim_track_energy.clear();
    m_sim_track_kinetic_energy.clear();
    m_sim_track_mass.clear();
    m_sim_track_pdg.clear();
    m_sim_track_origin.clear();
    m_sim_track_parent_pdg.clear();
    m_sim_track_process.clear();
    m_sim_track_startx.clear();
    m_sim_track_starty.clear();
    m_sim_track_startz.clear();
    m_sim_track_endx.clear();
    m_sim_track_endy.clear();
    m_sim_track_endz.clear();
    m_sim_track_length.clear();

    m_sim_track_px.clear();
    m_sim_track_py.clear();
    m_sim_track_pz.clear();
    m_sim_track_trackID.clear();

    // PID
    m_reco_track_pid_bragg_likelihood_mu_plane0.clear();
    m_reco_track_pid_bragg_likelihood_mu_plane1.clear();
    m_reco_track_pid_bragg_likelihood_mu_plane2.clear();
    m_reco_track_pid_bragg_likelihood_p_plane0.clear();
    m_reco_track_pid_bragg_likelihood_p_plane1.clear();
    m_reco_track_pid_bragg_likelihood_p_plane2.clear();
    m_reco_track_pid_bragg_likelihood_mip_plane0.clear();
    m_reco_track_pid_bragg_likelihood_mip_plane1.clear();
    m_reco_track_pid_bragg_likelihood_mip_plane2.clear();
    m_reco_track_pid_chi2_mu_plane0.clear();
    m_reco_track_pid_chi2_mu_plane1.clear();
    m_reco_track_pid_chi2_mu_plane2.clear();
    m_reco_track_pid_chi2_p_plane0.clear();
    m_reco_track_pid_chi2_p_plane1.clear();
    m_reco_track_pid_chi2_p_plane2.clear();
    m_reco_track_pid_pida_plane0.clear();
    m_reco_track_pid_pida_plane1.clear();
    m_reco_track_pid_pida_plane2.clear();
    m_reco_track_pid_three_plane_proton_pid.clear();

//    m_reco_track_end_to_nearest_dead_wire_plane0.clear();
//    m_reco_track_end_to_nearest_dead_wire_plane1.clear();
//    m_reco_track_end_to_nearest_dead_wire_plane2.clear();

    m_reco_track_sliceId.clear();
    m_reco_track_nuscore.clear();
    m_reco_track_isclearcosmic.clear();
    m_reco_track_trackscore.clear();
    m_reco_track_pfparticle_pdg.clear();
    m_reco_track_is_nuslice.clear();

    m_sim_track_sliceId.clear();
    m_sim_track_nuscore.clear();
    m_sim_track_isclearcosmic.clear();
  }

  void CreateTrackBranches(){
    vertex_tree->Branch("reco_asso_tracks",&m_reco_asso_tracks,"reco_asso_tracks/I");
    vertex_tree->Branch("reco_track_num_daughters",&m_reco_track_num_daughters);
    vertex_tree->Branch("reco_track_daughter_trackscore",&m_reco_track_daughter_trackscore);
    vertex_tree->Branch("reco_track_displacement", &m_reco_track_length);
    vertex_tree->Branch("reco_track_dirx", &m_reco_track_dirx);
    vertex_tree->Branch("reco_track_diry", &m_reco_track_diry);
    vertex_tree->Branch("reco_track_dirz", &m_reco_track_dirz);
    vertex_tree->Branch("reco_track_startx", &m_reco_track_startx);
    vertex_tree->Branch("reco_track_starty", &m_reco_track_starty);
    vertex_tree->Branch("reco_track_startz", &m_reco_track_startz);

    vertex_tree->Branch("reco_track_endx", &m_reco_track_endx);
    vertex_tree->Branch("reco_track_endy", &m_reco_track_endy);
    vertex_tree->Branch("reco_track_endz", &m_reco_track_endz);
    vertex_tree->Branch("reco_track_end_dist_to_active_TPC", &m_reco_track_end_dist_to_active_TPC);
    vertex_tree->Branch("reco_track_start_dist_to_active_TPC", &m_reco_track_start_dist_to_active_TPC);
    vertex_tree->Branch("reco_track_end_dist_to_CPA", &m_reco_track_end_dist_to_CPA);
    vertex_tree->Branch("reco_track_start_dist_to_CPA", &m_reco_track_start_dist_to_CPA);
    vertex_tree->Branch("reco_track_end_dist_to_SCB", &m_reco_track_end_dist_to_SCB);
    vertex_tree->Branch("reco_track_start_dist_to_SCB", &m_reco_track_start_dist_to_SCB);
    vertex_tree->Branch("reco_track_end_in_SCB", &m_reco_track_end_in_SCB);
    vertex_tree->Branch("reco_track_start_in_SCB", &m_reco_track_start_in_SCB);


    vertex_tree->Branch("reco_track_theta_yz", &m_reco_track_theta_yz);
    vertex_tree->Branch("reco_track_phi_yx", &m_reco_track_phi_yx);

    vertex_tree->Branch("reco_track_calo_energy_plane0", &m_reco_track_calo_energy_plane0);
    vertex_tree->Branch("reco_track_calo_energy_plane1", &m_reco_track_calo_energy_plane1);
    vertex_tree->Branch("reco_track_calo_energy_plane2", &m_reco_track_calo_energy_plane2);
    vertex_tree->Branch("reco_track_calo_energy_max", &m_reco_track_calo_energy_max);

    vertex_tree->Branch("reco_track_num_trajpoints", &m_reco_track_num_trajpoints);
    vertex_tree->Branch("reco_track_num_spacepoints", &m_reco_track_num_spacepoints);
    vertex_tree->Branch("reco_track_proton_kinetic_energy", &m_reco_track_proton_kinetic_energy);
    vertex_tree->Branch("reco_track_ordered_energy_index", &m_reco_track_ordered_energy_index);
    vertex_tree->Branch("reco_track_ordered_displacement_index", &m_reco_track_ordered_displacement_index);
    vertex_tree->Branch("i_trk", &m_reco_track_ordered_displacement_index);

    vertex_tree->Branch("reco_track_spacepoint_principal0",&m_reco_track_spacepoint_principal0);
    vertex_tree->Branch("reco_track_spacepoint_principal1",&m_reco_track_spacepoint_principal1);
    vertex_tree->Branch("reco_track_spacepoint_principal2",&m_reco_track_spacepoint_principal2);

    vertex_tree->Branch("reco_track_spacepoint_chi",&m_reco_track_spacepoint_chi);
    vertex_tree->Branch("reco_track_spacepoint_max_dist",&m_reco_track_spacepoint_max_dist);

    vertex_tree->Branch("reco_track_best_calo_plane",&m_reco_track_best_calo_plane);

    vertex_tree->Branch("reco_track_mean_dEdx_best_plane",&m_reco_track_mean_dEdx_best_plane);
    vertex_tree->Branch("reco_track_mean_dEdx_plane0",&m_reco_track_mean_dEdx_p0);
    vertex_tree->Branch("reco_track_mean_dEdx_plane1",&m_reco_track_mean_dEdx_p1);
    vertex_tree->Branch("reco_track_mean_dEdx_plane2",&m_reco_track_mean_dEdx_p2);

    vertex_tree->Branch("reco_track_mean_dEdx_start_half_best_plane",&m_reco_track_mean_dEdx_end_half_best_plane);
    vertex_tree->Branch("reco_track_mean_dEdx_start_half_plane0",&m_reco_track_mean_dEdx_end_half_p0);
    vertex_tree->Branch("reco_track_mean_dEdx_start_half_plane1",&m_reco_track_mean_dEdx_end_half_p1);
    vertex_tree->Branch("reco_track_mean_dEdx_start_half_plane2",&m_reco_track_mean_dEdx_end_half_p2);

    vertex_tree->Branch("reco_track_mean_dEdx_end_half_best_plane",&m_reco_track_mean_dEdx_start_half_best_plane);
    vertex_tree->Branch("reco_track_mean_dEdx_end_half_plane0",&m_reco_track_mean_dEdx_start_half_p0);
    vertex_tree->Branch("reco_track_mean_dEdx_end_half_plane1",&m_reco_track_mean_dEdx_start_half_p1);
    vertex_tree->Branch("reco_track_mean_dEdx_end_half_plane2",&m_reco_track_mean_dEdx_start_half_p2);

    vertex_tree->Branch("reco_track_good_calo_best_plane",&m_reco_track_good_calo_best_plane);
    vertex_tree->Branch("reco_track_good_calo_plane0",&m_reco_track_good_calo_p0);
    vertex_tree->Branch("reco_track_good_calo_plane1",&m_reco_track_good_calo_p1);
    vertex_tree->Branch("reco_track_good_calo_plane2",&m_reco_track_good_calo_p2);

    vertex_tree->Branch("reco_track_trunc_dEdx_best_plane",&m_reco_track_trunc_dEdx_best_plane);
    vertex_tree->Branch("reco_track_trunc_dEdx_plane0",&m_reco_track_trunc_dEdx_p0);
    vertex_tree->Branch("reco_track_trunc_dEdx_plane1",&m_reco_track_trunc_dEdx_p1);
    vertex_tree->Branch("reco_track_trunc_dEdx_plane2",&m_reco_track_trunc_dEdx_p2);

    vertex_tree->Branch("reco_track_mean_trunc_dEdx_best_plane",&m_reco_track_mean_trunc_dEdx_best_plane);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_plane0",&m_reco_track_mean_trunc_dEdx_p0);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_plane1",&m_reco_track_mean_trunc_dEdx_p1);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_plane2",&m_reco_track_mean_trunc_dEdx_p2);

    vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half_best_plane",&m_reco_track_mean_trunc_dEdx_end_half_best_plane);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half_plane0",&m_reco_track_mean_trunc_dEdx_end_half_p0);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half_plane1",&m_reco_track_mean_trunc_dEdx_end_half_p1);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half_plane2",&m_reco_track_mean_trunc_dEdx_end_half_p2);

    vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half_best_plane",&m_reco_track_mean_trunc_dEdx_start_half_best_plane);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half_plane0",&m_reco_track_mean_trunc_dEdx_start_half_p0);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half_plane1",&m_reco_track_mean_trunc_dEdx_start_half_p1);
    vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half_plane2",&m_reco_track_mean_trunc_dEdx_start_half_p2);

    vertex_tree->Branch("reco_track_trunc_PIDA_best_plane",&m_reco_track_trunc_PIDA_best_plane);
    vertex_tree->Branch("reco_track_trunc_PIDA_plane0",&m_reco_track_trunc_PIDA_p0);
    vertex_tree->Branch("reco_track_trunc_PIDA_plane1",&m_reco_track_trunc_PIDA_p1);
    vertex_tree->Branch("reco_track_trunc_PIDA_plane2",&m_reco_track_trunc_PIDA_p2);

    vertex_tree->Branch("reco_track_resrange_best_plane",&m_reco_track_resrange_best_plane);
    vertex_tree->Branch("reco_track_resrange_plane0",&m_reco_track_resrange_p0);
    vertex_tree->Branch("reco_track_resrange_plane1",&m_reco_track_resrange_p1);
    vertex_tree->Branch("reco_track_resrange_plane2",&m_reco_track_resrange_p2);

    vertex_tree->Branch("reco_track_dEdx_best_plane",&m_reco_track_dEdx_best_plane);
    vertex_tree->Branch("reco_track_dEdx_plane0",&m_reco_track_dEdx_p0);
    vertex_tree->Branch("reco_track_dEdx_plane1",&m_reco_track_dEdx_p1);
    vertex_tree->Branch("reco_track_dEdx_plane2",&m_reco_track_dEdx_p2);

    vertex_tree->Branch("reco_track_num_calo_hits_plane0",&m_reco_track_num_calo_hits_p0);
    vertex_tree->Branch("reco_track_num_calo_hits_plane1",&m_reco_track_num_calo_hits_p1);
    vertex_tree->Branch("reco_track_num_calo_hits_plane2",&m_reco_track_num_calo_hits_p2);



    vertex_tree->Branch("reco_track_pid_bragg_likelihood_mu_plane0",&m_reco_track_pid_bragg_likelihood_mu_plane0);
    vertex_tree->Branch("reco_track_pid_bragg_likelihood_mu_plane1",&m_reco_track_pid_bragg_likelihood_mu_plane1);
    vertex_tree->Branch("reco_track_pid_bragg_likelihood_mu_plane2",&m_reco_track_pid_bragg_likelihood_mu_plane2);
    vertex_tree->Branch("reco_track_pid_bragg_likelihood_p_plane0",&m_reco_track_pid_bragg_likelihood_p_plane0);
    vertex_tree->Branch("reco_track_pid_bragg_likelihood_p_plane1",&m_reco_track_pid_bragg_likelihood_p_plane1);
    vertex_tree->Branch("reco_track_pid_bragg_likelihood_p_plane2",&m_reco_track_pid_bragg_likelihood_p_plane2);
    vertex_tree->Branch("reco_track_pid_bragg_likelihood_mip_plane0",&m_reco_track_pid_bragg_likelihood_mip_plane0);
    vertex_tree->Branch("reco_track_pid_bragg_likelihood_mip_plane1",&m_reco_track_pid_bragg_likelihood_mip_plane1);
    vertex_tree->Branch("reco_track_pid_bragg_likelihood_mip_plane2",&m_reco_track_pid_bragg_likelihood_mip_plane2);
    vertex_tree->Branch("reco_track_pid_chi2_mu_plane0",&m_reco_track_pid_chi2_mu_plane0);
    vertex_tree->Branch("reco_track_pid_chi2_mu_plane1",&m_reco_track_pid_chi2_mu_plane1);
    vertex_tree->Branch("reco_track_pid_chi2_mu_plane2",&m_reco_track_pid_chi2_mu_plane2);
    vertex_tree->Branch("reco_track_pid_chi2_p_plane0",&m_reco_track_pid_chi2_p_plane0);
    vertex_tree->Branch("reco_track_pid_chi2_p_plane1",&m_reco_track_pid_chi2_p_plane1);
    vertex_tree->Branch("reco_track_pid_chi2_p_plane2",&m_reco_track_pid_chi2_p_plane2);
    vertex_tree->Branch("reco_track_pid_pida_plane0",&m_reco_track_pid_pida_plane0);
    vertex_tree->Branch("reco_track_pid_pida_plane1",&m_reco_track_pid_pida_plane1);
    vertex_tree->Branch("reco_track_pid_pida_plane2",&m_reco_track_pid_pida_plane2);
    vertex_tree->Branch("reco_track_pid_three_plane_proton_pid",&m_reco_track_pid_three_plane_proton_pid);

//    vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane0",&m_reco_track_end_to_nearest_dead_wire_plane0);
//    vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane1",&m_reco_track_end_to_nearest_dead_wire_plane1);
//    vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane2",&m_reco_track_end_to_nearest_dead_wire_plane2);

    vertex_tree->Branch("reco_track_sliceId",& m_reco_track_sliceId);
    vertex_tree->Branch("reco_track_nuscore",& m_reco_track_nuscore);
    vertex_tree->Branch("reco_track_isclearcosmic",& m_reco_track_isclearcosmic);
    vertex_tree->Branch("reco_track_trackscore",& m_reco_track_trackscore);
    vertex_tree->Branch("reco_track_pfparticle_pdg",& m_reco_track_pfparticle_pdg);
    vertex_tree->Branch("reco_track_is_nuslice",& m_reco_track_is_nuslice);

    vertex_tree->Branch("sim_track_matched",&m_sim_track_matched);
    vertex_tree->Branch("sim_track_overlay_fraction",&m_sim_track_overlay_fraction);
    vertex_tree->Branch("sim_track_energy",&m_sim_track_energy);
    vertex_tree->Branch("sim_track_kinetic_energy",&m_sim_track_kinetic_energy);
    vertex_tree->Branch("sim_track_mass",&m_sim_track_mass);
    vertex_tree->Branch("sim_track_pdg",&m_sim_track_pdg);
    vertex_tree->Branch("sim_track_parent_pdg",&m_sim_track_parent_pdg);
    vertex_tree->Branch("sim_track_origin",&m_sim_track_origin);
    vertex_tree->Branch("sim_track_process",&m_sim_track_process);
    vertex_tree->Branch("sim_track_startx",&m_sim_track_startx);
    vertex_tree->Branch("sim_track_starty",&m_sim_track_starty);
    vertex_tree->Branch("sim_track_startz",&m_sim_track_startz);
    vertex_tree->Branch("sim_track_px",&m_sim_track_px);
    vertex_tree->Branch("sim_track_py",&m_sim_track_py);
    vertex_tree->Branch("sim_track_pz",&m_sim_track_pz);
    vertex_tree->Branch("sim_track_endx",&m_sim_track_endx);
    vertex_tree->Branch("sim_track_endy",&m_sim_track_endy);
    vertex_tree->Branch("sim_track_endz",&m_sim_track_endz);
    vertex_tree->Branch("sim_track_length",&m_sim_track_length);

    vertex_tree->Branch("sim_track_trackID",&m_sim_track_trackID);

    vertex_tree->Branch("sim_track_sliceId",& m_sim_track_sliceId);
    vertex_tree->Branch("sim_track_nuscore",& m_sim_track_nuscore);
    vertex_tree->Branch("sim_track_isclearcosmic",& m_sim_track_isclearcosmic);
  }

  void ResizeTracks(size_t size){
    m_reco_track_length.resize(size);
    m_reco_track_dirx.resize(size);
    m_reco_track_num_daughters.resize(size);
    m_reco_track_daughter_trackscore.resize(size);

    m_reco_track_diry.resize(size);
    m_reco_track_dirz.resize(size);
    m_reco_track_endx.resize(size);
    m_reco_track_endy.resize(size);
    m_reco_track_endz.resize(size);
    m_reco_track_end_dist_to_active_TPC.resize(size);
    m_reco_track_start_dist_to_active_TPC.resize(size);
    m_reco_track_end_dist_to_CPA.resize(size);
    m_reco_track_start_dist_to_CPA.resize(size);
    m_reco_track_end_dist_to_SCB.resize(size);
    m_reco_track_start_dist_to_SCB.resize(size);
    m_reco_track_end_in_SCB.resize(size);
    m_reco_track_start_in_SCB.resize(size);

    m_reco_track_calo_energy_plane0.resize(size);
    m_reco_track_calo_energy_plane1.resize(size);
    m_reco_track_calo_energy_plane2.resize(size);
    m_reco_track_calo_energy_max.resize(size);



    m_reco_track_startx.resize(size);
    m_reco_track_starty.resize(size);
    m_reco_track_startz.resize(size);
    m_reco_track_num_trajpoints.resize(size);
    m_reco_track_num_spacepoints.resize(size);
    m_reco_track_proton_kinetic_energy.resize(size);
    m_reco_track_ordered_energy_index.resize(size);
    m_reco_track_ordered_displacement_index.resize(size);


    m_reco_track_spacepoint_principal0.resize(size);
    m_reco_track_spacepoint_principal1.resize(size);
    m_reco_track_spacepoint_principal2.resize(size);

    m_reco_track_spacepoint_chi.resize(size);
    m_reco_track_spacepoint_max_dist.resize(size);

    m_reco_track_theta_yz.resize(size);
    m_reco_track_phi_yx.resize(size);

    m_reco_track_best_calo_plane.resize(size);

    m_reco_track_mean_dEdx_best_plane.resize(size);
    m_reco_track_mean_dEdx_start_half_best_plane.resize(size);
    m_reco_track_mean_dEdx_end_half_best_plane.resize(size);
    m_reco_track_good_calo_best_plane.resize(size);
    m_reco_track_trunc_dEdx_best_plane.resize(size);
    m_reco_track_mean_trunc_dEdx_best_plane.resize(size);
    m_reco_track_mean_trunc_dEdx_start_half_best_plane.resize(size);
    m_reco_track_mean_trunc_dEdx_end_half_best_plane.resize(size);
    m_reco_track_trunc_PIDA_best_plane.resize(size);
    m_reco_track_resrange_best_plane.resize(size);
    m_reco_track_dEdx_best_plane.resize(size);


    m_reco_track_mean_dEdx_p0.resize(size);
    m_reco_track_mean_dEdx_start_half_p0.resize(size);
    m_reco_track_mean_dEdx_end_half_p0.resize(size);
    m_reco_track_good_calo_p0.resize(size);
    m_reco_track_trunc_dEdx_p0.resize(size);
    m_reco_track_mean_trunc_dEdx_p0.resize(size);
    m_reco_track_mean_trunc_dEdx_start_half_p0.resize(size);
    m_reco_track_mean_trunc_dEdx_end_half_p0.resize(size);
    m_reco_track_trunc_PIDA_p0.resize(size);
    m_reco_track_resrange_p0.resize(size);
    m_reco_track_dEdx_p0.resize(size);

    m_reco_track_mean_dEdx_p1.resize(size);
    m_reco_track_mean_dEdx_start_half_p1.resize(size);
    m_reco_track_mean_dEdx_end_half_p1.resize(size);
    m_reco_track_good_calo_p1.resize(size);
    m_reco_track_trunc_dEdx_p1.resize(size);
    m_reco_track_mean_trunc_dEdx_p1.resize(size);
    m_reco_track_mean_trunc_dEdx_start_half_p1.resize(size);
    m_reco_track_mean_trunc_dEdx_end_half_p1.resize(size);
    m_reco_track_trunc_PIDA_p1.resize(size);
    m_reco_track_resrange_p1.resize(size);
    m_reco_track_dEdx_p1.resize(size);

    m_reco_track_mean_dEdx_p2.resize(size);
    m_reco_track_mean_dEdx_start_half_p2.resize(size);
    m_reco_track_mean_dEdx_end_half_p2.resize(size);
    m_reco_track_good_calo_p2.resize(size);
    m_reco_track_trunc_dEdx_p2.resize(size);
    m_reco_track_mean_trunc_dEdx_p2.resize(size);
    m_reco_track_mean_trunc_dEdx_start_half_p2.resize(size);
    m_reco_track_mean_trunc_dEdx_end_half_p2.resize(size);
    m_reco_track_trunc_PIDA_p2.resize(size);
    m_reco_track_resrange_p2.resize(size);
    m_reco_track_dEdx_p2.resize(size);

    m_reco_track_num_calo_hits_p1.resize(size);
    m_reco_track_num_calo_hits_p0.resize(size);
    m_reco_track_num_calo_hits_p2.resize(size);



    m_sim_track_matched.resize(size);
    m_sim_track_energy.resize(size);
    m_sim_track_mass.resize(size);
    m_sim_track_kinetic_energy.resize(size);
    m_sim_track_pdg.resize(size);
    m_sim_track_parent_pdg.resize(size);
    m_sim_track_origin.resize(size);
    m_sim_track_process.resize(size);
    m_sim_track_startx.resize(size);
    m_sim_track_starty.resize(size);
    m_sim_track_startz.resize(size);
    m_sim_track_endx.resize(size);
    m_sim_track_endy.resize(size);
    m_sim_track_endz.resize(size);
    m_sim_track_length.resize(size);

    m_sim_track_px.resize(size);
    m_sim_track_py.resize(size);
    m_sim_track_pz.resize(size);
    m_sim_track_trackID.resize(size);
    m_sim_track_overlay_fraction.resize(size);

    m_reco_track_pid_bragg_likelihood_mu_plane0.resize(size);
    m_reco_track_pid_bragg_likelihood_mu_plane1.resize(size);
    m_reco_track_pid_bragg_likelihood_mu_plane2.resize(size);
    m_reco_track_pid_bragg_likelihood_p_plane0.resize(size);
    m_reco_track_pid_bragg_likelihood_p_plane1.resize(size);
    m_reco_track_pid_bragg_likelihood_p_plane2.resize(size);
    m_reco_track_pid_bragg_likelihood_mip_plane0.resize(size);
    m_reco_track_pid_bragg_likelihood_mip_plane1.resize(size);
    m_reco_track_pid_bragg_likelihood_mip_plane2.resize(size);
    m_reco_track_pid_chi2_mu_plane0.resize(size);
    m_reco_track_pid_chi2_mu_plane1.resize(size);
    m_reco_track_pid_chi2_mu_plane2.resize(size);
    m_reco_track_pid_chi2_p_plane0.resize(size);
    m_reco_track_pid_chi2_p_plane1.resize(size);
    m_reco_track_pid_chi2_p_plane2.resize(size);
    m_reco_track_pid_pida_plane0.resize(size);
    m_reco_track_pid_pida_plane1.resize(size);
    m_reco_track_pid_pida_plane2.resize(size);
    m_reco_track_pid_three_plane_proton_pid.resize(size);

//    m_reco_track_end_to_nearest_dead_wire_plane0.resize(size);
//    m_reco_track_end_to_nearest_dead_wire_plane1.resize(size);
//    m_reco_track_end_to_nearest_dead_wire_plane2.resize(size);

    m_reco_track_sliceId.resize(size);
    m_reco_track_nuscore.resize(size);
    m_reco_track_isclearcosmic.resize(size);
    m_reco_track_trackscore.resize(size);
    m_reco_track_pfparticle_pdg.resize(size);
    m_reco_track_is_nuslice.resize(size);

    m_sim_track_sliceId.resize(size);
    m_sim_track_nuscore.resize(size);
    m_sim_track_isclearcosmic.resize(size);
  }

  //analyze_Shower.h
  void ClearShowers(){
    m_reco_asso_showers=0;
    m_reco_shower_num_daughters.clear();
    m_reco_shower_daughter_trackscore.clear();

    m_reco_shower3d_exists.clear();

    m_reco_shower3d_startx.clear();
    m_reco_shower3d_starty.clear();
    m_reco_shower3d_startz.clear();
    m_reco_shower3d_dirx.clear();
    m_reco_shower3d_diry.clear();
    m_reco_shower3d_dirz.clear();
    m_reco_shower3d_theta_yz.clear();
    m_reco_shower3d_phi_yx.clear();
    m_reco_shower3d_conversion_distance.clear();
    m_reco_shower3d_impact_parameter.clear();
    m_reco_shower3d_implied_dirx.clear();
    m_reco_shower3d_implied_diry.clear();
    m_reco_shower3d_implied_dirz.clear();
    m_reco_shower3d_openingangle.clear();
    m_reco_shower3d_length.clear();

    m_reco_shower3d_energy_plane0.clear();
    m_reco_shower3d_energy_plane1.clear();
    m_reco_shower3d_energy_plane2.clear();
    m_reco_shower3d_dEdx_plane0.clear();
    m_reco_shower3d_dEdx_plane1.clear();
    m_reco_shower3d_dEdx_plane2.clear();


    m_reco_shower_startx.clear();
    m_reco_shower_starty.clear();
    m_reco_shower_start_dist_to_active_TPC.clear();
    m_reco_shower_start_dist_to_CPA.clear();
    m_reco_shower_start_dist_to_SCB.clear();
    m_reco_shower_start_in_SCB.clear();
    m_reco_shower_end_dist_to_active_TPC.clear();
    m_reco_shower_end_dist_to_SCB.clear();

    m_reco_shower_dirx.clear();
    m_reco_shower_diry.clear();
    m_reco_shower_dirz.clear();
    m_reco_shower_theta_yz.clear();
    m_reco_shower_phi_yx.clear();
    m_reco_shower_conversion_distance.clear();
    m_reco_shower_impact_parameter.clear();
    m_reco_shower_implied_dirx.clear();
    m_reco_shower_implied_diry.clear();
    m_reco_shower_implied_dirz.clear();
    m_reco_shower_openingangle.clear();
    m_reco_shower_length.clear();
    m_reco_shower_delaunay_num_triangles_plane0.clear();
    m_reco_shower_delaunay_num_triangles_plane1.clear();
    m_reco_shower_delaunay_num_triangles_plane2.clear();
    m_reco_shower_num_hits_plane0.clear();
    m_reco_shower_num_hits_plane1.clear();
    m_reco_shower_num_hits_plane2.clear();
    m_reco_shower_delaunay_area_plane0.clear();
    m_reco_shower_delaunay_area_plane1.clear();
    m_reco_shower_delaunay_area_plane2.clear();

    m_reco_shower_kalman_exists.clear();
    m_reco_shower_kalman_median_dEdx_plane0.clear();
    m_reco_shower_kalman_median_dEdx_plane1.clear();
    m_reco_shower_kalman_median_dEdx_plane2.clear();
    m_reco_shower_kalman_median_dEdx_allplane.clear();
    m_reco_shower_kalman_mean_dEdx_plane0.clear();
    m_reco_shower_kalman_mean_dEdx_plane1.clear();
    m_reco_shower_kalman_mean_dEdx_plane2.clear();

    m_sim_shower_energy.clear();
    m_sim_shower_matched.clear();
    m_sim_shower_kinetic_energy.clear();
    m_sim_shower_mass.clear();
    m_sim_shower_pdg.clear();
    m_sim_shower_trackID.clear();
    m_sim_shower_parent_pdg.clear();
    m_sim_shower_parent_trackID.clear();
    m_sim_shower_origin.clear();
    m_sim_shower_process.clear();
    m_sim_shower_end_process.clear();
    m_sim_shower_start_x.clear();
    m_sim_shower_start_y.clear();
    m_sim_shower_start_z.clear();
    m_sim_shower_vertex_x.clear();
    m_sim_shower_vertex_y.clear();
    m_sim_shower_vertex_z.clear();
    m_sim_shower_is_true_shower.clear();
    m_sim_shower_best_matched_plane.clear();
    m_sim_shower_matched_energy_fraction_plane0.clear();
    m_sim_shower_matched_energy_fraction_plane1.clear();
    m_sim_shower_matched_energy_fraction_plane2.clear();
    m_sim_shower_overlay_fraction.clear();
    m_sim_shower_px.clear();
    m_sim_shower_py.clear();
    m_sim_shower_pz.clear();
    m_sim_shower_sliceId.clear();
    m_sim_shower_nuscore.clear();
    m_sim_shower_isclearcosmic.clear();
    m_sim_shower_is_nuslice.clear();



    m_reco_shower_ordered_energy_index.clear();
    m_reco_shower_energy_max.clear();
    m_reco_shower_energy_plane0.clear();
    m_reco_shower_energy_plane1.clear();
    m_reco_shower_energy_plane2.clear();

    m_reco_shower_reclustered_energy_plane0.clear();
    m_reco_shower_reclustered_energy_plane1.clear();
    m_reco_shower_reclustered_energy_plane2.clear();
    m_reco_shower_reclustered_energy_max.clear();

    m_reco_shower_plane0_nhits.clear();
    m_reco_shower_plane1_nhits.clear();
    m_reco_shower_plane2_nhits.clear();
    m_reco_shower_plane0_meanRMS.clear();
    m_reco_shower_plane1_meanRMS.clear();
    m_reco_shower_plane2_meanRMS.clear();

    m_reco_shower_hit_tick.clear();
    m_reco_shower_hit_wire.clear();
    m_reco_shower_hit_plane.clear();
    m_reco_shower_spacepoint_x.clear();
    m_reco_shower_spacepoint_y.clear();
    m_reco_shower_spacepoint_z.clear();


    m_reco_shower_dQdx_plane0.clear();
    m_reco_shower_dQdx_plane2.clear();
    m_reco_shower_dQdx_plane2.clear();
    m_reco_shower_dEdx_plane0.clear();
    m_reco_shower_dEdx_plane1.clear();
    m_reco_shower_dEdx_plane2.clear();
    m_reco_shower_dEdx_plane0_median.clear();
    m_reco_shower_dEdx_plane1_median.clear();
    m_reco_shower_dEdx_plane2_median.clear();

    m_reco_shower_angle_wrt_wires_plane0.clear();
    m_reco_shower_angle_wrt_wires_plane1.clear();
    m_reco_shower_angle_wrt_wires_plane2.clear();

    m_reco_shower_dEdx_amalgamated.clear();
    m_reco_shower_dEdx_amalgamated_nhits.clear();


    m_reco_shower_dQdx_plane0_median.clear();
    m_reco_shower_dQdx_plane1_median.clear();
    m_reco_shower_dQdx_plane2_median.clear();

    m_reco_shower_dEdx_plane0_mean.clear();
    m_reco_shower_dEdx_plane1_mean.clear();
    m_reco_shower_dEdx_plane2_mean.clear();  
    m_reco_shower_dEdx_plane0_max.clear();
    m_reco_shower_dEdx_plane1_max.clear();
    m_reco_shower_dEdx_plane2_max.clear();  
    m_reco_shower_dEdx_plane0_min.clear();
    m_reco_shower_dEdx_plane1_min.clear();
    m_reco_shower_dEdx_plane2_min.clear();  

    m_reco_shower_dEdx_plane0_nhits.clear();
    m_reco_shower_dEdx_plane1_nhits.clear();
    m_reco_shower_dEdx_plane2_nhits.clear();  

//    m_reco_shower_start_to_nearest_dead_wire_plane0.clear();
//    m_reco_shower_start_to_nearest_dead_wire_plane1.clear();
//    m_reco_shower_start_to_nearest_dead_wire_plane2.clear();

    m_reco_shower_flash_shortest_distz.clear();
    m_reco_shower_flash_shortest_index_z.clear();
    m_reco_shower_flash_shortest_disty.clear();
    m_reco_shower_flash_shortest_index_y.clear();

    m_reco_shower_flash_shortest_distyz.clear();
    m_reco_shower_flash_shortest_index_yz.clear();

    m_reco_shower_sliceId.clear();
    m_reco_shower_nuscore.clear();
    m_reco_shower_isclearcosmic.clear();
    m_reco_shower_is_nuslice.clear();
    m_reco_shower_trackscore.clear();
    m_reco_shower_pfparticle_pdg.clear();

  }

  void CreateShowerBranches(){
    vertex_tree->Branch("reco_asso_showers",&m_reco_asso_showers,"reco_asso_showers/I");
    vertex_tree->Branch("reco_shower_num_daughters",&m_reco_shower_num_daughters);
    vertex_tree->Branch("reco_shower_daughter_trackscore",&m_reco_shower_daughter_trackscore);

    vertex_tree->Branch("reco_shower_length", &m_reco_shower_length);
    vertex_tree->Branch("reco_shower_opening_angle", &m_reco_shower_openingangle);
    vertex_tree->Branch("reco_shower_dirx", &m_reco_shower_dirx);
    vertex_tree->Branch("reco_shower_diry", &m_reco_shower_diry);
    vertex_tree->Branch("reco_shower_dirz", &m_reco_shower_dirz);
    vertex_tree->Branch("reco_shower_startx", &m_reco_shower_startx);
    vertex_tree->Branch("reco_shower_starty", &m_reco_shower_starty);
    vertex_tree->Branch("reco_shower_startz", &m_reco_shower_startz);
    vertex_tree->Branch("reco_shower_start_dist_to_active_TPC", &m_reco_shower_start_dist_to_active_TPC);
    vertex_tree->Branch("reco_shower_start_dist_to_CPA", &m_reco_shower_start_dist_to_CPA);
    vertex_tree->Branch("reco_shower_start_dist_to_SCB",  &m_reco_shower_start_dist_to_SCB);
    vertex_tree->Branch("reco_shower_start_in_SCB",   &m_reco_shower_start_in_SCB);
    vertex_tree->Branch("reco_shower_end_dist_to_active_TPC", &m_reco_shower_end_dist_to_active_TPC);
    vertex_tree->Branch("reco_shower_end_dist_to_SCB",  &m_reco_shower_end_dist_to_SCB);


    vertex_tree->Branch("reco_shower_theta_yz",&m_reco_shower_theta_yz);
    vertex_tree->Branch("reco_shower_phi_yx",&m_reco_shower_phi_yx);
    vertex_tree->Branch("reco_shower_conversion_distance",& m_reco_shower_conversion_distance);
    vertex_tree->Branch("reco_shower_impact_parameter",& m_reco_shower_impact_parameter);
    vertex_tree->Branch("reco_shower_implied_dirx", &m_reco_shower_implied_dirx);
    vertex_tree->Branch("reco_shower_implied_diry", &m_reco_shower_implied_diry);
    vertex_tree->Branch("reco_shower_implied_dirz", &m_reco_shower_implied_dirz);

    vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane0",&m_reco_shower_delaunay_num_triangles_plane0);
    vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane1",&m_reco_shower_delaunay_num_triangles_plane1);
    vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane2",&m_reco_shower_delaunay_num_triangles_plane2);
    vertex_tree->Branch("reco_shower_num_hits_plane0",&m_reco_shower_num_hits_plane0);
    vertex_tree->Branch("reco_shower_num_hits_plane1",&m_reco_shower_num_hits_plane1);
    vertex_tree->Branch("reco_shower_num_hits_plane2",&m_reco_shower_num_hits_plane2);
    vertex_tree->Branch("reco_shower_delaunay_area_plane0",&m_reco_shower_delaunay_area_plane0);
    vertex_tree->Branch("reco_shower_delaunay_area_plane1",&m_reco_shower_delaunay_area_plane1);
    vertex_tree->Branch("reco_shower_delaunay_area_plane2",&m_reco_shower_delaunay_area_plane2);
    //the calorimetry info
    vertex_tree->Branch("reco_shower_energy_max",&m_reco_shower_energy_max);
    vertex_tree->Branch("reco_shower_energy_plane0",&m_reco_shower_energy_plane0);
    vertex_tree->Branch("reco_shower_energy_plane1",&m_reco_shower_energy_plane1);
    vertex_tree->Branch("reco_shower_energy_plane2",&m_reco_shower_energy_plane2);
    vertex_tree->Branch("reco_shower_plane0_nhits",&m_reco_shower_plane0_nhits);
    vertex_tree->Branch("reco_shower_plane1_nhits",&m_reco_shower_plane1_nhits);
    vertex_tree->Branch("reco_shower_plane2_nhits",&m_reco_shower_plane2_nhits);
    vertex_tree->Branch("reco_shower_plane0_meanRMS",&m_reco_shower_plane0_meanRMS);
    vertex_tree->Branch("reco_shower_plane1_meanRMS",&m_reco_shower_plane1_meanRMS);
    vertex_tree->Branch("reco_shower_plane2_meanRMS",&m_reco_shower_plane2_meanRMS);

    vertex_tree->Branch("reco_shower_reclustered_energy_plane0",&m_reco_shower_reclustered_energy_plane0);
    vertex_tree->Branch("reco_shower_reclustered_energy_plane1",&m_reco_shower_reclustered_energy_plane1);
    vertex_tree->Branch("reco_shower_reclustered_energy_plane2",&m_reco_shower_reclustered_energy_plane2);
    vertex_tree->Branch("reco_shower_reclustered_energy_max",&m_reco_shower_reclustered_energy_max);

    vertex_tree->Branch("reco_shower_hit_tick",&m_reco_shower_hit_tick);
    vertex_tree->Branch("reco_shower_hit_wire",&m_reco_shower_hit_wire);
    vertex_tree->Branch("reco_shower_hit_plane",&m_reco_shower_hit_plane);

    vertex_tree->Branch("reco_shower_spacepoint_x",&m_reco_shower_spacepoint_x);
    vertex_tree->Branch("reco_shower_spacepoint_y",&m_reco_shower_spacepoint_y);
    vertex_tree->Branch("reco_shower_spacepoint_z",&m_reco_shower_spacepoint_z);

    vertex_tree->Branch("reco_shower_ordered_energy_index",&m_reco_shower_ordered_energy_index);
    vertex_tree->Branch("i_shr",&m_reco_shower_ordered_energy_index);
    vertex_tree->Branch("reco_shower_dQdx_plane0",&m_reco_shower_dQdx_plane0);
    vertex_tree->Branch("reco_shower_dQdx_plane1",&m_reco_shower_dQdx_plane1);
    vertex_tree->Branch("reco_shower_dQdx_plane2",&m_reco_shower_dQdx_plane2);
    vertex_tree->Branch("reco_shower_dEdx_plane0",&m_reco_shower_dEdx_plane0);
    vertex_tree->Branch("reco_shower_dEdx_plane1",&m_reco_shower_dEdx_plane1);
    vertex_tree->Branch("reco_shower_dEdx_plane2",&m_reco_shower_dEdx_plane2);
    vertex_tree->Branch("reco_shower_dEdx_plane0_median",&m_reco_shower_dEdx_plane0_median);
    vertex_tree->Branch("reco_shower_dEdx_plane1_median",&m_reco_shower_dEdx_plane1_median);
    vertex_tree->Branch("reco_shower_dEdx_plane2_median",&m_reco_shower_dEdx_plane2_median);

    vertex_tree->Branch("reco_shower_angle_wrt_wires_plane0",& m_reco_shower_angle_wrt_wires_plane0);
    vertex_tree->Branch("reco_shower_angle_wrt_wires_plane1",& m_reco_shower_angle_wrt_wires_plane1);
    vertex_tree->Branch("reco_shower_angle_wrt_wires_plane2",& m_reco_shower_angle_wrt_wires_plane2);

    vertex_tree->Branch("reco_shower_dEdx_amalgamated",&m_reco_shower_dEdx_amalgamated);
    vertex_tree->Branch("reco_shower_dEdx_amalgamated_nhits",&m_reco_shower_dEdx_amalgamated_nhits);


    vertex_tree->Branch("reco_shower_dQdx_plane0_median",&m_reco_shower_dQdx_plane0_median);
    vertex_tree->Branch("reco_shower_dQdx_plane1_median",&m_reco_shower_dQdx_plane1_median);
    vertex_tree->Branch("reco_shower_dQdx_plane2_median",&m_reco_shower_dQdx_plane2_median);

    vertex_tree->Branch("reco_shower_dEdx_plane0_mean",&m_reco_shower_dEdx_plane0_mean);
    vertex_tree->Branch("reco_shower_dEdx_plane1_mean",&m_reco_shower_dEdx_plane1_mean);
    vertex_tree->Branch("reco_shower_dEdx_plane2_mean",&m_reco_shower_dEdx_plane2_mean);
    vertex_tree->Branch("reco_shower_dEdx_plane0_max",&m_reco_shower_dEdx_plane0_max);
    vertex_tree->Branch("reco_shower_dEdx_plane1_max",&m_reco_shower_dEdx_plane1_max);
    vertex_tree->Branch("reco_shower_dEdx_plane2_max",&m_reco_shower_dEdx_plane2_max);
    vertex_tree->Branch("reco_shower_dEdx_plane0_min",&m_reco_shower_dEdx_plane0_min);
    vertex_tree->Branch("reco_shower_dEdx_plane1_min",&m_reco_shower_dEdx_plane1_min);
    vertex_tree->Branch("reco_shower_dEdx_plane2_min",&m_reco_shower_dEdx_plane2_min);
    vertex_tree->Branch("reco_shower_dEdx_plane0_nhits",&m_reco_shower_dEdx_plane0_nhits);
    vertex_tree->Branch("reco_shower_dEdx_plane1_nhits",&m_reco_shower_dEdx_plane1_nhits);
    vertex_tree->Branch("reco_shower_dEdx_plane2_nhits",&m_reco_shower_dEdx_plane2_nhits);

//    vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane0",&m_reco_shower_start_to_nearest_dead_wire_plane0);
//    vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane1",&m_reco_shower_start_to_nearest_dead_wire_plane1);
//    vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane2",&m_reco_shower_start_to_nearest_dead_wire_plane2);

    vertex_tree->Branch("reco_shower_flash_shortest_distz",&m_reco_shower_flash_shortest_distz);
    vertex_tree->Branch("reco_shower_flash_shortest_disty",&m_reco_shower_flash_shortest_disty);
    vertex_tree->Branch("reco_shower_flash_shortest_distyz",&m_reco_shower_flash_shortest_distyz);
    vertex_tree->Branch("reco_shower_flash_shortest_index_z",&m_reco_shower_flash_shortest_index_z);
    vertex_tree->Branch("reco_shower_flash_shortest_index_y",&m_reco_shower_flash_shortest_index_y);
    vertex_tree->Branch("reco_shower_flash_shortest_index_yz",&m_reco_shower_flash_shortest_index_yz);

    vertex_tree->Branch("reco_shower_sliceId",& m_reco_shower_sliceId);
    vertex_tree->Branch("reco_shower_nuscore",& m_reco_shower_nuscore);
    vertex_tree->Branch("reco_shower_isclearcosmic",& m_reco_shower_isclearcosmic);
    vertex_tree->Branch("reco_shower_is_nuslice", & m_reco_shower_is_nuslice);
    vertex_tree->Branch("reco_shower_trackscore", & m_reco_shower_trackscore);
    vertex_tree->Branch("reco_shower_pfparticle_pdg", & m_reco_shower_pfparticle_pdg);


    vertex_tree->Branch("reco_shower3d_exists", &m_reco_shower3d_exists);
    vertex_tree->Branch("reco_shower3d_length", &m_reco_shower3d_length);
    vertex_tree->Branch("reco_shower3d_opening_angle", &m_reco_shower3d_openingangle);
    vertex_tree->Branch("reco_shower3d_dirx", &m_reco_shower3d_dirx);
    vertex_tree->Branch("reco_shower3d_diry", &m_reco_shower3d_diry);
    vertex_tree->Branch("reco_shower3d_dirz", &m_reco_shower3d_dirz);
    vertex_tree->Branch("reco_shower3d_startx", &m_reco_shower3d_startx);
    vertex_tree->Branch("reco_shower3d_starty", &m_reco_shower3d_starty);
    vertex_tree->Branch("reco_shower3d_startz", &m_reco_shower3d_startz);
    vertex_tree->Branch("reco_shower3d_theta_yz",&m_reco_shower3d_theta_yz);
    vertex_tree->Branch("reco_shower3d_phi_yx",&m_reco_shower3d_phi_yx);
    vertex_tree->Branch("reco_shower3d_conversion_distance",& m_reco_shower3d_conversion_distance);
    vertex_tree->Branch("reco_shower3d_impact_parameter",& m_reco_shower3d_impact_parameter);
    vertex_tree->Branch("reco_shower3d_implied_dirx", &m_reco_shower3d_implied_dirx);
    vertex_tree->Branch("reco_shower3d_implied_diry", &m_reco_shower3d_implied_diry);
    vertex_tree->Branch("reco_shower3d_implied_dirz", &m_reco_shower3d_implied_dirz);

    vertex_tree->Branch("reco_shower3d_energy_plane0", &m_reco_shower3d_energy_plane0);
    vertex_tree->Branch("reco_shower3d_energy_plane1", &m_reco_shower3d_energy_plane1);
    vertex_tree->Branch("reco_shower3d_energy_plane2", &m_reco_shower3d_energy_plane2);
    vertex_tree->Branch("reco_shower3d_dEdx_plane0", &m_reco_shower3d_dEdx_plane0);
    vertex_tree->Branch("reco_shower3d_dEdx_plane1", &m_reco_shower3d_dEdx_plane1);
    vertex_tree->Branch("reco_shower3d_dEdx_plane2", &m_reco_shower3d_dEdx_plane2);

    vertex_tree->Branch("reco_shower_kalman_exists",&m_reco_shower_kalman_exists);
    vertex_tree->Branch("reco_shower_kalman_dEdx_plane0_median",&m_reco_shower_kalman_median_dEdx_plane0);
    vertex_tree->Branch("reco_shower_kalman_dEdx_plane1_median",&m_reco_shower_kalman_median_dEdx_plane1);
    vertex_tree->Branch("reco_shower_kalman_dEdx_plane2_median",&m_reco_shower_kalman_median_dEdx_plane2);
    vertex_tree->Branch("reco_shower_kalman_dEdx_allplane_median",&m_reco_shower_kalman_median_dEdx_allplane);

    vertex_tree->Branch("reco_shower_kalman_dEdx_plane0_mean",&m_reco_shower_kalman_mean_dEdx_plane0);
    vertex_tree->Branch("reco_shower_kalman_dEdx_plane1_mean",&m_reco_shower_kalman_mean_dEdx_plane1);
    vertex_tree->Branch("reco_shower_kalman_dEdx_plane2_mean",&m_reco_shower_kalman_mean_dEdx_plane2);


    vertex_tree->Branch("sim_shower_matched",&m_sim_shower_matched);
    vertex_tree->Branch("sim_shower_energy",&m_sim_shower_energy);
    vertex_tree->Branch("sim_shower_kinetic_energy",&m_sim_shower_kinetic_energy);
    vertex_tree->Branch("sim_shower_mass",&m_sim_shower_mass);
    vertex_tree->Branch("sim_shower_pdg",&m_sim_shower_pdg);
    vertex_tree->Branch("sim_shower_trackID",&m_sim_shower_trackID);
    vertex_tree->Branch("sim_shower_parent_pdg",&m_sim_shower_parent_pdg);
    vertex_tree->Branch("sim_shower_parent_trackID",&m_sim_shower_parent_trackID);
    vertex_tree->Branch("sim_shower_origin",&m_sim_shower_origin);
    vertex_tree->Branch("sim_shower_process",&m_sim_shower_process);
    vertex_tree->Branch("sim_shower_end_process",&m_sim_shower_end_process);
    vertex_tree->Branch("sim_shower_start_x",&m_sim_shower_start_x);
    vertex_tree->Branch("sim_shower_start_y",&m_sim_shower_start_y);
    vertex_tree->Branch("sim_shower_start_z",&m_sim_shower_start_z);
    vertex_tree->Branch("sim_shower_vertex_x",&m_sim_shower_vertex_x);
    vertex_tree->Branch("sim_shower_vertex_y",&m_sim_shower_vertex_y);
    vertex_tree->Branch("sim_shower_vertex_z",&m_sim_shower_vertex_z);
    vertex_tree->Branch("sim_shower_px",&m_sim_shower_px);
    vertex_tree->Branch("sim_shower_py",&m_sim_shower_py);
    vertex_tree->Branch("sim_shower_pz",&m_sim_shower_pz);

    vertex_tree->Branch("sim_shower_is_true_shower",&m_sim_shower_is_true_shower);
    vertex_tree->Branch("sim_shower_best_matched_plane",&m_sim_shower_best_matched_plane);
    vertex_tree->Branch("sim_shower_matched_energy_fraction_plane0",&m_sim_shower_matched_energy_fraction_plane0);
    vertex_tree->Branch("sim_shower_matched_energy_fraction_plane1",&m_sim_shower_matched_energy_fraction_plane1);
    vertex_tree->Branch("sim_shower_matched_energy_fraction_plane2",&m_sim_shower_matched_energy_fraction_plane2);
    vertex_tree->Branch("sim_shower_overlay_fraction",&m_sim_shower_overlay_fraction);
    vertex_tree->Branch("sim_shower_sliceId", & m_sim_shower_sliceId);
    vertex_tree->Branch("sim_shower_nuscore", & m_sim_shower_nuscore);
    vertex_tree->Branch("sim_shower_isclearcosmic", & m_sim_shower_isclearcosmic);
    vertex_tree->Branch("sim_shower_is_nusclice", & m_sim_shower_is_nuslice);
  }

  void ResizeShowers(size_t size){
    m_reco_shower_num_daughters.resize(size);
    m_reco_shower_daughter_trackscore.resize(size);

    m_reco_shower_kalman_exists.resize(size);
    m_reco_shower_kalman_median_dEdx_plane0.resize(size);
    m_reco_shower_kalman_median_dEdx_plane1.resize(size);
    m_reco_shower_kalman_median_dEdx_plane2.resize(size);
    m_reco_shower_kalman_median_dEdx_allplane.resize(size);
    m_reco_shower_kalman_mean_dEdx_plane0.resize(size);
    m_reco_shower_kalman_mean_dEdx_plane1.resize(size);
    m_reco_shower_kalman_mean_dEdx_plane2.resize(size);

    m_reco_shower_reclustered_energy_plane0.resize(size);
    m_reco_shower_reclustered_energy_plane1.resize(size);
    m_reco_shower_reclustered_energy_plane2.resize(size);
    m_reco_shower_reclustered_energy_max.resize(size);


    m_reco_shower3d_exists.resize(size);
    m_reco_shower3d_startx.resize(size);
    m_reco_shower3d_starty.resize(size);
    m_reco_shower3d_startz.resize(size);
    m_reco_shower3d_dirx.resize(size);
    m_reco_shower3d_diry.resize(size);
    m_reco_shower3d_dirz.resize(size);
    m_reco_shower3d_theta_yz.resize(size);
    m_reco_shower3d_phi_yx.resize(size);
    m_reco_shower3d_conversion_distance.resize(size);
    m_reco_shower3d_openingangle.resize(size);
    m_reco_shower3d_length.resize(size);
    m_reco_shower3d_impact_parameter.resize(size);
    m_reco_shower3d_implied_dirx.resize(size);
    m_reco_shower3d_implied_diry.resize(size);
    m_reco_shower3d_implied_dirz.resize(size);
    m_reco_shower3d_energy_plane0.resize(size);
    m_reco_shower3d_energy_plane1.resize(size);
    m_reco_shower3d_energy_plane2.resize(size);
    m_reco_shower3d_dEdx_plane0.resize(size);
    m_reco_shower3d_dEdx_plane1.resize(size);
    m_reco_shower3d_dEdx_plane2.resize(size);

    m_reco_shower_start_dist_to_active_TPC.resize(size);
    m_reco_shower_start_dist_to_CPA.resize(size);
    m_reco_shower_start_dist_to_SCB.resize(size);
    m_reco_shower_start_in_SCB.resize(size);

    m_reco_shower_end_dist_to_active_TPC.resize(size);
    m_reco_shower_end_dist_to_SCB.resize(size);


    m_reco_shower_startx.resize(size);
    m_reco_shower_starty.resize(size);
    m_reco_shower_startz.resize(size);
    m_reco_shower_dirx.resize(size);
    m_reco_shower_diry.resize(size);
    m_reco_shower_dirz.resize(size);
    m_reco_shower_theta_yz.resize(size);
    m_reco_shower_phi_yx.resize(size);
    m_reco_shower_conversion_distance.resize(size);
    m_reco_shower_openingangle.resize(size);
    m_reco_shower_length.resize(size);
    m_reco_shower_impact_parameter.resize(size);
    m_reco_shower_implied_dirx.resize(size);
    m_reco_shower_implied_diry.resize(size);
    m_reco_shower_implied_dirz.resize(size);
    m_reco_shower_delaunay_num_triangles_plane0.resize(size);
    m_reco_shower_delaunay_num_triangles_plane1.resize(size);
    m_reco_shower_delaunay_num_triangles_plane2.resize(size);
    m_reco_shower_num_hits_plane0.resize(size);
    m_reco_shower_num_hits_plane1.resize(size);
    m_reco_shower_num_hits_plane2.resize(size);
    m_reco_shower_delaunay_area_plane0.resize(size);
    m_reco_shower_delaunay_area_plane1.resize(size);
    m_reco_shower_delaunay_area_plane2.resize(size);

    m_reco_shower_energy_max.resize(size);
    m_reco_shower_energy_plane0.resize(size);
    m_reco_shower_energy_plane1.resize(size);
    m_reco_shower_energy_plane2.resize(size);

    m_reco_shower_plane0_nhits.resize(size);
    m_reco_shower_plane1_nhits.resize(size);
    m_reco_shower_plane2_nhits.resize(size);

    m_reco_shower_plane0_meanRMS.resize(size);
    m_reco_shower_plane1_meanRMS.resize(size);
    m_reco_shower_plane2_meanRMS.resize(size);



    m_reco_shower_ordered_energy_index.resize(size);
    m_reco_shower_dQdx_plane0.resize(size);
    m_reco_shower_dQdx_plane1.resize(size);
    m_reco_shower_dQdx_plane2.resize(size);
    m_reco_shower_dEdx_plane0.resize(size);
    m_reco_shower_dEdx_plane1.resize(size);
    m_reco_shower_dEdx_plane2.resize(size);
    m_reco_shower_dEdx_plane0_median.resize(size);
    m_reco_shower_dEdx_plane1_median.resize(size);
    m_reco_shower_dEdx_plane2_median.resize(size);

    m_reco_shower_angle_wrt_wires_plane0.resize(size);
    m_reco_shower_angle_wrt_wires_plane1.resize(size);
    m_reco_shower_angle_wrt_wires_plane2.resize(size);

    m_reco_shower_dEdx_amalgamated.resize(size);
    m_reco_shower_dEdx_amalgamated_nhits.resize(size);

    m_reco_shower_dQdx_plane0_median.resize(size);
    m_reco_shower_dQdx_plane1_median.resize(size);
    m_reco_shower_dQdx_plane2_median.resize(size);

    m_reco_shower_dEdx_plane0_min.resize(size);
    m_reco_shower_dEdx_plane1_min.resize(size);
    m_reco_shower_dEdx_plane2_min.resize(size);
    m_reco_shower_dEdx_plane0_max.resize(size);
    m_reco_shower_dEdx_plane1_max.resize(size);
    m_reco_shower_dEdx_plane2_max.resize(size);
    m_reco_shower_dEdx_plane0_mean.resize(size);
    m_reco_shower_dEdx_plane1_mean.resize(size);
    m_reco_shower_dEdx_plane2_mean.resize(size);




    m_reco_shower_dEdx_plane0_nhits.resize(size);
    m_reco_shower_dEdx_plane1_nhits.resize(size);
    m_reco_shower_dEdx_plane2_nhits.resize(size);

//    m_reco_shower_start_to_nearest_dead_wire_plane0.resize(size);
//    m_reco_shower_start_to_nearest_dead_wire_plane1.resize(size);
//    m_reco_shower_start_to_nearest_dead_wire_plane2.resize(size);

    m_reco_shower_flash_shortest_distz.resize(size);
    m_reco_shower_flash_shortest_index_z.resize(size);
    m_reco_shower_flash_shortest_disty.resize(size);
    m_reco_shower_flash_shortest_index_y.resize(size);

    m_reco_shower_flash_shortest_distyz.resize(size);
    m_reco_shower_flash_shortest_index_yz.resize(size);

    m_reco_shower_sliceId.resize(size);
    m_reco_shower_nuscore.resize(size);
    m_reco_shower_isclearcosmic.resize(size);
    m_reco_shower_is_nuslice.resize(size);
    m_reco_shower_trackscore.resize(size);
    m_reco_shower_pfparticle_pdg.resize(size);


    m_sim_shower_energy.resize(size);
    m_sim_shower_matched.resize(size);
    m_sim_shower_kinetic_energy.resize(size);
    m_sim_shower_mass.resize(size);
    m_sim_shower_pdg.resize(size);
    m_sim_shower_trackID.resize(size);
    m_sim_shower_parent_pdg.resize(size);
    m_sim_shower_parent_trackID.resize(size);
    m_sim_shower_origin.resize(size);
    m_sim_shower_process.resize(size);
    m_sim_shower_end_process.resize(size);
    m_sim_shower_start_x.resize(size);
    m_sim_shower_start_y.resize(size);
    m_sim_shower_start_z.resize(size);
    m_sim_shower_vertex_x.resize(size);
    m_sim_shower_vertex_y.resize(size);
    m_sim_shower_vertex_z.resize(size);
    m_sim_shower_is_true_shower.resize(size);
    m_sim_shower_best_matched_plane.resize(size);
    m_sim_shower_matched_energy_fraction_plane0.resize(size);
    m_sim_shower_matched_energy_fraction_plane1.resize(size);
    m_sim_shower_matched_energy_fraction_plane2.resize(size);
    m_sim_shower_overlay_fraction.resize(size);
    m_sim_shower_px.resize(size);
    m_sim_shower_py.resize(size);
    m_sim_shower_pz.resize(size);
    m_sim_shower_sliceId.resize(size);
    m_sim_shower_nuscore.resize(size);
    m_sim_shower_isclearcosmic.resize(size);
    m_sim_shower_is_nuslice.resize(size);
  }

  //analyze_MCTruth.h
  void ClearMCTruths(){
    m_mctruth_num = 0;
    m_mctruth_origin = -99;
    m_mctruth_mode = -99;
    m_mctruth_interaction_type = -99;
    m_mctruth_nu_vertex_x = -9999;
    m_mctruth_nu_vertex_y = -9999;
    m_mctruth_nu_vertex_z = -9999;
    m_mctruth_reco_vertex_dist = -9999;
    m_mctruth_ccnc = -99;
    m_mctruth_qsqr = -99;
    m_mctruth_nu_E = -99;
    m_mctruth_nu_pdg = 0;
    m_mctruth_lepton_pdg = 0;
    m_mctruth_num_daughter_particles = -99;
    m_mctruth_daughters_pdg.clear();
    m_mctruth_daughters_E.clear();

    m_mctruth_daughters_status_code.clear();
    m_mctruth_daughters_trackID.clear();
    m_mctruth_daughters_mother_trackID.clear();
    m_mctruth_daughters_px.clear();
    m_mctruth_daughters_py.clear();
    m_mctruth_daughters_pz.clear();
    m_mctruth_daughters_startx.clear();
    m_mctruth_daughters_starty.clear();
    m_mctruth_daughters_startz.clear();
    m_mctruth_daughters_time.clear();
    m_mctruth_daughters_endx.clear();
    m_mctruth_daughters_endy.clear();
    m_mctruth_daughters_endz.clear();
    m_mctruth_daughters_endtime.clear();
    m_mctruth_daughters_process.clear();
    m_mctruth_daughters_end_process.clear();


    m_mctruth_is_delta_radiative = 0;
    m_mctruth_delta_radiative_1g1p_or_1g1n = -999;

    m_mctruth_delta_photon_energy=-999;
    m_mctruth_delta_proton_energy=-999;
    m_mctruth_delta_neutron_energy=-999;

    m_mctruth_num_exiting_photons =0;
    m_mctruth_num_exiting_protons =0;
    m_mctruth_num_exiting_pi0 =0;
    m_mctruth_num_exiting_pipm =0;
    m_mctruth_num_exiting_neutrons=0;
    m_mctruth_num_exiting_delta0=0;
    m_mctruth_num_exiting_deltapm=0;
    m_mctruth_num_exiting_deltapp=0;

    m_mctruth_num_reconstructable_protons = 0;

    m_mctruth_is_reconstructable_1g1p = 0;
    m_mctruth_is_reconstructable_1g0p = 0;

    m_mctruth_leading_exiting_proton_energy = -9999;

    m_mctruth_exiting_pi0_E.clear();
    m_mctruth_exiting_pi0_mom.clear();
    m_mctruth_exiting_pi0_px.clear();
    m_mctruth_exiting_pi0_py.clear();
    m_mctruth_exiting_pi0_pz.clear();

    m_mctruth_pi0_leading_photon_energy = -9999;
    m_mctruth_pi0_subleading_photon_energy = -9999;
    m_mctruth_pi0_leading_photon_end_process = "none";
    m_mctruth_pi0_subleading_photon_end_process = "none";
    m_mctruth_pi0_leading_photon_end = {-9999,-9999,-9999};
    m_mctruth_pi0_leading_photon_start = {-9999,-9999,-9999};
    m_mctruth_pi0_subleading_photon_end = {-9999,-9999,-9999};
    m_mctruth_pi0_subleading_photon_start = {-9999,-9999,-9999};
    m_mctruth_pi0_leading_photon_exiting_TPC = -999;
    m_mctruth_pi0_subleading_photon_exiting_TPC = -999;
    m_mctruth_pi0_leading_photon_mom = {-9999,-9999,-9999};
    m_mctruth_pi0_subleading_photon_mom = {-9999,-9999,-9999};

    m_mctruth_exiting_delta0_num_daughters.clear();

    m_mctruth_exiting_photon_mother_trackID.clear();
    m_mctruth_exiting_photon_trackID.clear();
    m_mctruth_exiting_photon_from_delta_decay.clear();
    m_mctruth_exiting_photon_energy.clear();
    m_mctruth_exiting_photon_px.clear();
    m_mctruth_exiting_photon_py.clear();
    m_mctruth_exiting_photon_pz.clear();

    m_mctruth_exiting_proton_mother_trackID.clear();
    m_mctruth_exiting_proton_trackID.clear();
    m_mctruth_exiting_proton_from_delta_decay.clear();
    m_mctruth_exiting_proton_energy.clear();
    m_mctruth_exiting_proton_px.clear();
    m_mctruth_exiting_proton_py.clear();
    m_mctruth_exiting_proton_pz.clear();

    m_mctruth_exiting_neutron_mother_trackID.clear();
    m_mctruth_exiting_neutron_trackID.clear();
    m_mctruth_exiting_neutron_from_delta_decay.clear();
    m_mctruth_exiting_neutron_energy.clear();
    m_mctruth_exiting_neutron_px.clear();
    m_mctruth_exiting_neutron_py.clear();
    m_mctruth_exiting_neutron_pz.clear();
  }

  void CreateMCTruthBranches(){
    vertex_tree->Branch("mctruth_num",&m_mctruth_num);
    vertex_tree->Branch("mctruth_origin",&m_mctruth_origin);
    vertex_tree->Branch("mctruth_nu_pdg",&m_mctruth_nu_pdg);
    vertex_tree->Branch("mctruth_nu_E",&m_mctruth_nu_E);

    vertex_tree->Branch("mctruth_nu_vertex_x",&m_mctruth_nu_vertex_x);
    vertex_tree->Branch("mctruth_nu_vertex_y",&m_mctruth_nu_vertex_y);
    vertex_tree->Branch("mctruth_nu_vertex_z",&m_mctruth_nu_vertex_z);
    vertex_tree->Branch("mctruth_reco_vertex_dist",&m_mctruth_reco_vertex_dist);

    vertex_tree->Branch("mctruth_lepton_pdg",&m_mctruth_lepton_pdg);
    vertex_tree->Branch("mctruth_lepton_E",&m_mctruth_lepton_E);
    vertex_tree->Branch("mctruth_mode",&m_mctruth_mode);
    vertex_tree->Branch("mctruth_qsqr",&m_mctruth_qsqr);
    vertex_tree->Branch("mctruth_cc_or_nc",&m_mctruth_ccnc);
    vertex_tree->Branch("mctruth_interaction_type",&m_mctruth_interaction_type);

    vertex_tree->Branch("mctruth_num_daughter_particles",&m_mctruth_num_daughter_particles);
    vertex_tree->Branch("mctruth_daughters_pdg",&m_mctruth_daughters_pdg);
    vertex_tree->Branch("mctruth_daughters_E",&m_mctruth_daughters_E);
    vertex_tree->Branch("mctruth_daughters_status_code",&m_mctruth_daughters_status_code);
    vertex_tree->Branch("mctruth_daughters_trackID",&m_mctruth_daughters_trackID);
    vertex_tree->Branch("mctruth_daughters_mother_trackID",&m_mctruth_daughters_mother_trackID);
    vertex_tree->Branch("mctruth_daughters_px",&m_mctruth_daughters_px);
    vertex_tree->Branch("mctruth_daughters_py",&m_mctruth_daughters_py);
    vertex_tree->Branch("mctruth_daughters_pz",&m_mctruth_daughters_pz);
    vertex_tree->Branch("mctruth_daughters_startx",&m_mctruth_daughters_startx);
    vertex_tree->Branch("mctruth_daughters_starty",&m_mctruth_daughters_starty);
    vertex_tree->Branch("mctruth_daughters_startz",&m_mctruth_daughters_startz);
    vertex_tree->Branch("mctruth_daughters_time",&m_mctruth_daughters_time);
    vertex_tree->Branch("mctruth_daughters_endx",&m_mctruth_daughters_endx);
    vertex_tree->Branch("mctruth_daughters_endy",&m_mctruth_daughters_endy);
    vertex_tree->Branch("mctruth_daughters_endz",&m_mctruth_daughters_endz);
    vertex_tree->Branch("mctruth_daughters_endtime",&m_mctruth_daughters_endtime);
    vertex_tree->Branch("mctruth_daughters_process",&m_mctruth_daughters_process);
    vertex_tree->Branch("mctruth_daughters_end_process",&m_mctruth_daughters_end_process);




    vertex_tree->Branch("mctruth_num_exiting_protons",&m_mctruth_num_exiting_protons);
    vertex_tree->Branch("mctruth_num_exiting_photons",&m_mctruth_num_exiting_photons);
    vertex_tree->Branch("mctruth_num_exiting_neutrons",&m_mctruth_num_exiting_neutrons);
    vertex_tree->Branch("mctruth_num_exiting_pi0",&m_mctruth_num_exiting_pi0);
    vertex_tree->Branch("mctruth_num_exiting_pipm",&m_mctruth_num_exiting_pipm);
    vertex_tree->Branch("mctruth_num_exiting_delta0",&m_mctruth_num_exiting_delta0);
    vertex_tree->Branch("mctruth_num_exiting_deltapm",&m_mctruth_num_exiting_deltapm);
    vertex_tree->Branch("mctruth_num_exiting_deltapp",&m_mctruth_num_exiting_deltapp);

    vertex_tree->Branch("mctruth_leading_exiting_proton_energy",&m_mctruth_leading_exiting_proton_energy);
    vertex_tree->Branch("mctruth_is_delta_radiative",&m_mctruth_is_delta_radiative);
    vertex_tree->Branch("mctruth_delta_radiative_1g1p_or_1g1n",&m_mctruth_delta_radiative_1g1p_or_1g1n);
    vertex_tree->Branch("mctruth_delta_photon_energy",&m_mctruth_delta_photon_energy);
    vertex_tree->Branch("mctruth_delta_proton_energy",&m_mctruth_delta_proton_energy);
    vertex_tree->Branch("mctruth_delta_neutron_energy",&m_mctruth_delta_neutron_energy);
    vertex_tree->Branch("mctruth_exiting_delta0_num_daughters",&m_mctruth_exiting_delta0_num_daughters);

    vertex_tree->Branch("mctruth_exiting_photon_trackID",&m_mctruth_exiting_photon_trackID);
    vertex_tree->Branch("mctruth_exiting_photon_mother_trackID",&m_mctruth_exiting_photon_mother_trackID);
    vertex_tree->Branch("mctruth_exiting_photon_from_delta_decay",&m_mctruth_exiting_photon_from_delta_decay);
    vertex_tree->Branch("mctruth_exiting_photon_energy",&m_mctruth_exiting_photon_energy);
    vertex_tree->Branch("mctruth_exiting_photon_px",&m_mctruth_exiting_photon_px);
    vertex_tree->Branch("mctruth_exiting_photon_py",&m_mctruth_exiting_photon_py);
    vertex_tree->Branch("mctruth_exiting_photon_pz",&m_mctruth_exiting_photon_pz);

    vertex_tree->Branch("mctruth_exiting_proton_trackID",&m_mctruth_exiting_proton_trackID);
    vertex_tree->Branch("mctruth_exiting_proton_mother_trackID",&m_mctruth_exiting_proton_mother_trackID);
    vertex_tree->Branch("mctruth_exiting_proton_from_delta_decay",&m_mctruth_exiting_proton_from_delta_decay);
    vertex_tree->Branch("mctruth_exiting_proton_energy",&m_mctruth_exiting_proton_energy);
    vertex_tree->Branch("mctruth_exiting_proton_px",&m_mctruth_exiting_proton_px);
    vertex_tree->Branch("mctruth_exiting_proton_py",&m_mctruth_exiting_proton_py);
    vertex_tree->Branch("mctruth_exiting_proton_pz",&m_mctruth_exiting_proton_pz);

    vertex_tree->Branch("mctruth_exiting_neutron_trackID",&m_mctruth_exiting_neutron_trackID);
    vertex_tree->Branch("mctruth_exiting_neutron_mother_trackID",&m_mctruth_exiting_neutron_mother_trackID);
    vertex_tree->Branch("mctruth_exiting_neutron_from_delta_decay",&m_mctruth_exiting_neutron_from_delta_decay);
    vertex_tree->Branch("mctruth_exiting_neutron_energy",&m_mctruth_exiting_neutron_energy);
    vertex_tree->Branch("mctruth_exiting_neutron_px",&m_mctruth_exiting_neutron_px);
    vertex_tree->Branch("mctruth_exiting_neutron_py",&m_mctruth_exiting_neutron_py);
    vertex_tree->Branch("mctruth_exiting_neutron_pz",&m_mctruth_exiting_neutron_pz);


    vertex_tree->Branch("mctruth_is_reconstructable_1g1p",&m_mctruth_is_reconstructable_1g1p);

    vertex_tree->Branch("mctruth_is_reconstructable_1g1p",&m_mctruth_is_reconstructable_1g1p);
    vertex_tree->Branch("mctruth_num_reconstructable_protons",&m_mctruth_num_reconstructable_protons);

    vertex_tree->Branch("mctruth_pi0_leading_photon_energy",&m_mctruth_pi0_leading_photon_energy);
    vertex_tree->Branch("mctruth_pi0_leading_photon_mom",&m_mctruth_pi0_leading_photon_mom);
    vertex_tree->Branch("mctruth_pi0_leading_photon_start",&m_mctruth_pi0_leading_photon_start);
    vertex_tree->Branch("mctruth_pi0_leading_photon_end",&m_mctruth_pi0_leading_photon_end);
    vertex_tree->Branch("mctruth_pi0_leading_photon_exiting_TPC",&m_mctruth_pi0_leading_photon_exiting_TPC);
    vertex_tree->Branch("mctruth_pi0_subleading_photon_energy",&m_mctruth_pi0_subleading_photon_energy);
    vertex_tree->Branch("mctruth_pi0_subleading_photon_mom",&m_mctruth_pi0_subleading_photon_mom);
    vertex_tree->Branch("mctruth_pi0_subleading_photon_end_process",&m_mctruth_pi0_subleading_photon_end_process);
    vertex_tree->Branch("mctruth_pi0_subleading_photon_start",&m_mctruth_pi0_subleading_photon_start);
    vertex_tree->Branch("mctruth_pi0_subleading_photon_end",&m_mctruth_pi0_subleading_photon_end);
    vertex_tree->Branch("mctruth_pi0_subleading_photon_exiting_TPC",&m_mctruth_pi0_subleading_photon_exiting_TPC);


    vertex_tree->Branch("mctruth_exiting_pi0_E",&m_mctruth_exiting_pi0_E);
    vertex_tree->Branch("mctruth_exiting_pi0_mom",&m_mctruth_exiting_pi0_mom);
    vertex_tree->Branch("mctruth_exiting_pi0_px",&m_mctruth_exiting_pi0_px);
    vertex_tree->Branch("mctruth_exiting_pi0_py",&m_mctruth_exiting_pi0_py);
    vertex_tree->Branch("mctruth_exiting_pi0_pz",&m_mctruth_exiting_pi0_pz);
  }

  void ResizeMCTruths(size_t size){
    m_mctruth_daughters_pdg.resize(size);
    m_mctruth_daughters_E.resize(size);
    m_mctruth_daughters_status_code.resize(size);
    m_mctruth_daughters_trackID.resize(size);
    m_mctruth_daughters_mother_trackID.resize(size);
    m_mctruth_daughters_px.resize(size);
    m_mctruth_daughters_py.resize(size);
    m_mctruth_daughters_pz.resize(size);
    m_mctruth_daughters_startx.resize(size);
    m_mctruth_daughters_starty.resize(size);
    m_mctruth_daughters_startz.resize(size);
    m_mctruth_daughters_time.resize(size);
    m_mctruth_daughters_endx.resize(size);
    m_mctruth_daughters_endy.resize(size);
    m_mctruth_daughters_endz.resize(size);
    m_mctruth_daughters_endtime.resize(size);
    m_mctruth_daughters_end_process.resize(size);
    m_mctruth_daughters_process.resize(size);
  }

  //analyze_EventWeight.h
  void ClearEventWeightBranches(){
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

  void CreateEventWeightBranches(){
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
  }

  //analyze_Geant4.h
  void ClearGeant4Branches(var_geant4 m_collection_geant4){

    m_collection_geant4.m_geant4_pdg.clear();
    m_collection_geant4.m_geant4_trackid.clear();
    m_collection_geant4.m_geant4_mother.clear();
    m_collection_geant4.m_geant4_statuscode.clear();
    m_collection_geant4.m_geant4_E.clear();
    m_collection_geant4.m_geant4_mass.clear();
    m_collection_geant4.m_geant4_px.clear();
    m_collection_geant4.m_geant4_py.clear();
    m_collection_geant4.m_geant4_pz.clear();
    m_collection_geant4.m_geant4_dx.clear();
    m_collection_geant4.m_geant4_dy.clear();
    m_collection_geant4.m_geant4_dz.clear();

    m_collection_geant4.m_geant4_vx.clear();
    m_collection_geant4.m_geant4_vy.clear();
    m_collection_geant4.m_geant4_vz.clear();
    m_collection_geant4.m_geant4_process.clear();
    m_collection_geant4.m_geant4_end_process.clear();

    m_collection_geant4.m_geant4_costheta.clear();
  }

  void CreateGeant4Branches(var_geant4 m_collection_geant4){
    geant4_tree->Branch("geant4_pdg",&m_collection_geant4.m_geant4_pdg);
    geant4_tree->Branch("geant4_trackid",&m_collection_geant4.m_geant4_trackid);
    geant4_tree->Branch("geant4_mother",&m_collection_geant4.m_geant4_mother);
    geant4_tree->Branch("geant4_statuscode",&m_collection_geant4.m_geant4_statuscode);
    geant4_tree->Branch("geant4_E",&m_collection_geant4.m_geant4_E);
    geant4_tree->Branch("geant4_mass",&m_collection_geant4.m_geant4_mass);
    geant4_tree->Branch("geant4_px", &m_collection_geant4.m_geant4_px);
    geant4_tree->Branch("geant4_py", &m_collection_geant4.m_geant4_py);
    geant4_tree->Branch("geant4_pz", &m_collection_geant4.m_geant4_pz);

    geant4_tree->Branch("geant4_dx", &m_collection_geant4.m_geant4_dx);
    geant4_tree->Branch("geant4_dy", &m_collection_geant4.m_geant4_dy);
    geant4_tree->Branch("geant4_dz", &m_collection_geant4.m_geant4_dz);

    geant4_tree->Branch("geant4_vx", &m_collection_geant4.m_geant4_vx);
    geant4_tree->Branch("geant4_vy", &m_collection_geant4.m_geant4_vy);
    geant4_tree->Branch("geant4_vz", &m_collection_geant4.m_geant4_vz);
    geant4_tree->Branch("geant4_costheta",&m_collection_geant4.m_geant4_costheta);

    geant4_tree->Branch("geant4_end_process", &m_collection_geant4.m_geant4_end_process);
    geant4_tree->Branch("geant4_process", &m_collection_geant4.m_geant4_process);
  }

  //analyze_Slice.h
  void ClearSlices(){
    m_reco_slice_num = 0;
    m_reco_slice_nuscore.clear();
    m_matched_signal_shower_overlay_fraction.clear();
    //std::vector<double> m_matched_signal_shower_conversion_length;
    m_matched_signal_shower_true_E.clear();
    m_matched_signal_shower_nuscore.clear();
    m_matched_signal_shower_sliceId.clear();
    m_matched_signal_shower_is_clearcosmic.clear();
    m_matched_signal_shower_num = 0;
    m_matched_signal_shower_is_nuslice.clear();
    m_matched_signal_shower_tracks_in_slice.clear();
    m_matched_signal_shower_showers_in_slice.clear();

    m_reco_slice_num_pfps.clear();
    m_reco_slice_num_showers.clear();
    m_reco_slice_num_tracks.clear();


    m_matched_signal_track_true_E.clear();
    m_matched_signal_track_nuscore.clear();
    m_matched_signal_track_sliceId.clear();
    m_matched_signal_track_is_clearcosmic.clear();
    m_matched_signal_track_is_nuslice.clear();
    m_matched_signal_track_tracks_in_slice.clear();
    m_matched_signal_track_showers_in_slice.clear();


    m_matched_signal_track_num = 0;  


    //int m_matched_signal_total_num_slices;

    m_reco_1g1p_is_same_slice = false;
    m_reco_1g1p_is_multiple_slices = false;
    m_reco_1g1p_is_nuslice = false;
    m_reco_1g0p_is_nuslice = false;
    m_reco_1g1p_nuscore = -999;
    m_reco_1g0p_nuscore = -999;
    m_is_matched_1g1p = false;
    m_is_matched_1g0p = false;
    m_no_matched_showers = false;
    m_multiple_matched_showers = false;
    m_multiple_matched_tracks = false;


    /*  m_reco_slice_shower_num_matched_signal = -999;
      m_reco_slice_track_num_matched_signal = -999;
      m_reco_slice_shower_matched_sliceId.clear();
      m_reco_slice_track_matched_sliceId.clear();
      m_reco_slice_shower_matched_energy.clear();
      m_reco_slice_track_matched_energy.clear();
      m_reco_slice_shower_matched_conversion.clear();
      m_reco_slice_shower_matched_overlay_frac.clear();
      */  
  }


  void CreateSliceBranches(){
    vertex_tree->Branch("reco_slice_nuscore",&m_reco_slice_nuscore);
    vertex_tree->Branch("reco_slice_num",&m_reco_slice_num);
    vertex_tree->Branch("reco_slice_shower_num_matched_signal",& m_reco_slice_shower_num_matched_signal);
    vertex_tree->Branch("reco_slice_track_num_matched_signal",& m_reco_slice_track_num_matched_signal);

    ncdelta_slice_tree->Branch("matched_signal_shower_overlay_fraction", &m_matched_signal_shower_overlay_fraction);
    //std::vector<double> m_matched_signal_shower_conversion_length;
    ncdelta_slice_tree->Branch("matched_signal_shower_true_E", &m_matched_signal_shower_true_E);
    ncdelta_slice_tree->Branch("matched_signal_shower_nuscore", &m_matched_signal_shower_nuscore);
    ncdelta_slice_tree->Branch("matched_signal_shower_sliceId", &m_matched_signal_shower_sliceId);
    ncdelta_slice_tree->Branch("matched_signal_shower_is_clearcosmic", &m_matched_signal_shower_is_clearcosmic);
    ncdelta_slice_tree->Branch("matched_signal_shower_num", &m_matched_signal_shower_num);
    ncdelta_slice_tree->Branch("matched_signal_shower_is_nuslice", &m_matched_signal_shower_is_nuslice);
    ncdelta_slice_tree->Branch("matched_signal_shower_tracks_in_slice", &m_matched_signal_shower_tracks_in_slice);
    ncdelta_slice_tree->Branch("matched_signal_shower_showers_in_slice", &m_matched_signal_shower_showers_in_slice);

    ncdelta_slice_tree->Branch("reco_slice_num_pfps", & m_reco_slice_num_pfps);
    ncdelta_slice_tree->Branch("reco_slice_num_showers", & m_reco_slice_num_showers);
    ncdelta_slice_tree->Branch("reco_slice_num_tracks", & m_reco_slice_num_tracks);

    // ncdelta_slice_tree->Branch("matched_signal_track_overlay_fraction", &m_matched_signal_track_overlay_fraction);
    ncdelta_slice_tree->Branch("matched_signal_track_true_E", &m_matched_signal_track_true_E);
    ncdelta_slice_tree->Branch("matched_signal_track_nuscore", &m_matched_signal_track_nuscore);
    ncdelta_slice_tree->Branch("matched_signal_track_sliceId", &m_matched_signal_track_sliceId);
    ncdelta_slice_tree->Branch("matched_signal_track_is_clearcosmic", &m_matched_signal_track_is_clearcosmic);
    ncdelta_slice_tree->Branch("matched_signal_track_num", &m_matched_signal_track_num);
    ncdelta_slice_tree->Branch("matched_signal_track_is_nuslice", &m_matched_signal_track_is_nuslice);
    ncdelta_slice_tree->Branch("matched_signal_track_tracks_in_slice", &m_matched_signal_track_tracks_in_slice);
    ncdelta_slice_tree->Branch("matched_signal_track_showers_in_slice", &m_matched_signal_track_showers_in_slice);

    //int m_matched_signal_total_num_slices;
    ncdelta_slice_tree->Branch("reco_1g1p_is_same_slice",&m_reco_1g1p_is_same_slice);
    ncdelta_slice_tree->Branch("reco_1g1p_is_nuslice",&m_reco_1g1p_is_nuslice);
    ncdelta_slice_tree->Branch("reco_1g1p_is_multiple_slices",&m_reco_1g1p_is_multiple_slices);
    ncdelta_slice_tree->Branch("reco_1g1p_nuscore",&m_reco_1g1p_nuscore);
    ncdelta_slice_tree->Branch("is_matched_1g1p",&m_is_matched_1g1p);

    ncdelta_slice_tree->Branch("reco_1g0p_nuscore",&m_reco_1g0p_nuscore);
    ncdelta_slice_tree->Branch("reco_1g0p_is_nuslice",&m_reco_1g0p_is_nuslice);
    ncdelta_slice_tree->Branch("is_matched_1g0p",&m_is_matched_1g0p);

    ncdelta_slice_tree->Branch("no_matched_showers",& m_no_matched_showers);
    ncdelta_slice_tree->Branch("multiple_matched_showers",& m_multiple_matched_showers);
    ncdelta_slice_tree->Branch("multiple_matched_tracks",& m_multiple_matched_tracks);
  }

  void ResizeSlices(size_t size){
    m_reco_slice_nuscore.resize(size,-999);
    m_reco_slice_num_pfps.resize(size,0);
    m_reco_slice_num_showers.resize(size,0);
    m_reco_slice_num_tracks.resize(size,0);
  }

  void Output_EventMeta( art::Event &evt){
    
    //Some event based properties
        m_subrun_counts++;
        m_number_of_events++;
        m_run_number  = evt.run();
        m_subrun_number = evt.subRun();
        m_event_number  = evt.id().event();
  }
  
  //set the vertex for now;
  void Output_PFParticleInfo( std::vector<PandoraPFParticle> PPFPs){

    int pfp_size = PPFPs.size();

    for(int index = 0; index < pfp_size; index++){

      PandoraPFParticle* temp_p = &PPFPs[index];
      if(!(pfp_w_bestnuID == temp_p->get_SliceID() && temp_p->get_IsNeutrino()) ) continue;
      m_vertex_pos_x = temp_p->get_Vertex_pos()[0];
      m_vertex_pos_y = temp_p->get_Vertex_pos()[1];
      m_vertex_pos_z = temp_p->get_Vertex_pos()[2];
//      std::cout<<"Best NuScore is found, define the vertice as: ("<<temp_p->get_Vertex_pos()[0]<<","<<temp_p->get_Vertex_pos()[1]<<","<<temp_p->get_Vertex_pos()[2]<<")"<<std::endl;
    
      std::vector<double> tmp = {m_vertex_pos_x, m_vertex_pos_y, m_vertex_pos_z};
      m_reco_vertex_in_SCB = distToSCB(m_reco_vertex_dist_to_SCB,tmp);
      m_reco_vertex_dist_to_active_TPC =  distToTPCActive(tmp);
      m_reco_vertex_dist_to_CPA =  distToCPA(tmp);

      if(temp_p->get_IsNeutrino() ){ 
        m_reco_slice_num++;
        m_reco_slice_nuscore.push_back(temp_p->get_NuScore());

      }
    }

    //resize slice variables size;
    //    this->ResizeSlices(m_reco_slice_num); 
  }
}
