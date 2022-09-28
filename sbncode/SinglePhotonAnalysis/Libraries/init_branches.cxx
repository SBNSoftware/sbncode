#include "sbncode/SinglePhotonAnalysis/Libraries/init_branches.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/fiducial_volume.h"

namespace single_photon
{
  //Process of initialize branches:
  //ClearBranches, ResizeBranches, CreateBranches

  void ClearMeta(var_all& vars){
    //------------ Event related Variables -------------
    vars.m_event_number = -99;
    vars.m_subrun_number = -99;
    vars.m_run_number = -99;
    vars.m_test_matched_hits = 0;

    vars.m_pot_per_event = 0;
    vars.m_pot_per_subrun = vars.m_subrun_pot;
    vars.m_number_of_events_in_subrun = 0;

    vars.m_genie_spline_weight = 1.0;

    //------------ Vertex related Variables -------------
    vars.m_reco_vertex_size = 0;
    vars.m_vertex_pos_x=-99999;
    vars.m_vertex_pos_y=-99999;
    vars.m_vertex_pos_z=-99999;
    vars.m_vertex_pos_tick=-9999;
    vars.m_vertex_pos_wire_p0=-9999;
    vars.m_vertex_pos_wire_p1=-9999;
    vars.m_vertex_pos_wire_p2=-9999;
    vars.m_reco_vertex_in_SCB = -9999;
    vars.m_reco_vertex_dist_to_SCB = -9999;
    vars.m_reco_vertex_dist_to_active_TPC= -9999;
    vars.m_reco_vertex_dist_to_CPA= -9999;

//    vars.m_reco_vertex_to_nearest_dead_wire_plane0=-99999;
//    vars.m_reco_vertex_to_nearest_dead_wire_plane1=-99999;
//    vars.m_reco_vertex_to_nearest_dead_wire_plane2=-99999;

    vars.m_reco_slice_objects = 0;
  }


  void CreateMetaBranches(var_all& vars){

    //true_eventweight_tree
//    vars.true_eventweight_tree->Branch("mcweight", "std::map<std::string, std::vector<double>>",&vars.fmcweight);

    //run_subrun_tree
    vars.run_subrun_tree->Branch("run",&vars.m_run,"run/I");
    vars.run_subrun_tree->Branch("subrun",&vars.m_subrun,"subrun/I");
    vars.run_subrun_tree->Branch("subrun_pot",&vars.m_subrun_pot,"subrun_pot/D");
    vars.run_subrun_tree->Branch("subrun_counts",&vars.m_subrun_counts,"subrun_counts/I");

    //pot_tree
    vars.pot_tree->Branch("number_of_events",&vars.m_number_of_events,"number_of_events/I");
    vars.pot_tree->Branch("number_of_vertices",&vars.m_number_of_vertices,"number_of_vertices/I");
    vars.pot_tree->Branch("POT",&vars.m_pot_count,"POT/D");

    //vars.vertex_tree -- part of it
    // --------------------- Event Related variables ------------
    vars.vertex_tree->Branch("run_number", &vars.m_run_number, "run_number/I");
    vars.vertex_tree->Branch("subrun_number", &vars.m_subrun_number, "subrun_number/I");
    vars.vertex_tree->Branch("event_number", &vars.m_event_number, "event_number/I");

    vars.vertex_tree->Branch("pot_per_event",&vars.m_pot_per_event,"pot_per_event/D");
    vars.vertex_tree->Branch("pot_per_subrun",&vars.m_pot_per_subrun,"pot_per_subrun/D");
    vars.vertex_tree->Branch("number_of_events_in_subrun",&vars.m_number_of_events_in_subrun,"number_of_events_in_subrun/D");


    vars.vertex_tree->Branch("genie_spline_weight", &vars.m_genie_spline_weight, "genie_spline_weight/D");
    vars.vertex_tree->Branch("genie_CV_tune_weight", &vars.m_genie_CV_tune_weight, "genie_CV_tune_weight/D");

    vars.vertex_tree->Branch("photonu_weight_low", &vars.m_photonu_weight_low, "photonu_weight_low/D");
    vars.vertex_tree->Branch("photonu_weight_high", &vars.m_photonu_weight_high, "photonu_weight_high/D");

    vars.vertex_tree->Branch("test_matched_hits", &vars.m_test_matched_hits, "test_matched_hits/I");

    // --------------------- Vertex Related variables ------------
    vars.vertex_tree->Branch("reco_vertex_size", &vars.m_reco_vertex_size);
    vars.vertex_tree->Branch("reco_vertex_x", &vars.m_vertex_pos_x);
    vars.vertex_tree->Branch("reco_vertex_y", &vars.m_vertex_pos_y);
    vars.vertex_tree->Branch("reco_vertex_z", &vars.m_vertex_pos_z);
    vars.vertex_tree->Branch("reco_vertex_in_SCB", &vars.m_reco_vertex_in_SCB);
    vars.vertex_tree->Branch("reco_vertex_dist_to_SCB",&vars.m_reco_vertex_dist_to_SCB);
    vars.vertex_tree->Branch("reco_vertex_dist_to_active_TPC",&vars.m_reco_vertex_dist_to_active_TPC);
    vars.vertex_tree->Branch("reco_vertex_dist_to_CPA",&vars.m_reco_vertex_dist_to_CPA);
//  vars.vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane0",&vars.m_reco_vertex_to_nearest_dead_wire_plane0);
//  vars.vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane1",&vars.m_reco_vertex_to_nearest_dead_wire_plane1);
//  vars.vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane2",&vars.m_reco_vertex_to_nearest_dead_wire_plane2);

    vars.vertex_tree->Branch("reco_slice_objects", &vars.m_reco_slice_objects, "reco_slice_objects/I");

    vars.vertex_tree->Branch("vars.m_flash_optfltr_pe_beam",&vars.m_flash_optfltr_pe_beam);
    vars.vertex_tree->Branch("vars.m_flash_optfltr_pe_veto",&vars.m_flash_optfltr_pe_veto);
    vars.vertex_tree->Branch("vars.m_flash_optfltr_pe_veto_tot",&vars.m_flash_optfltr_pe_veto_tot);
    vars.vertex_tree->Branch("vars.m_flash_optfltr_pe_beavars.m_tot",&vars.m_flash_optfltr_pe_beam_tot);
  }

  //isolation.h
  void ClearIsolation(var_all& vars){
    vars.m_isolation_min_dist_trk_shr.clear();
    vars.m_isolation_min_dist_trk_unassoc.clear();

    vars.m_isolation_num_shr_hits_win_1cm_trk.clear();
    vars.m_isolation_num_shr_hits_win_2cm_trk.clear();
    vars.m_isolation_num_shr_hits_win_5cm_trk.clear();
    vars.m_isolation_num_shr_hits_win_10cm_trk.clear();

    vars.m_isolation_num_unassoc_hits_win_1cm_trk.clear();
    vars.m_isolation_num_unassoc_hits_win_2cm_trk.clear();
    vars.m_isolation_num_unassoc_hits_win_5cm_trk.clear();
    vars.m_isolation_num_unassoc_hits_win_10cm_trk.clear();

    vars.m_isolation_nearest_shr_hit_to_trk_wire.clear();
    vars.m_isolation_nearest_shr_hit_to_trk_time.clear();

    vars.m_isolation_nearest_unassoc_hit_to_trk_wire.clear();
    vars.m_isolation_nearest_unassoc_hit_to_trk_time.clear();
  }

  void CreateIsolationBranches(var_all& vars){
    vars.vertex_tree->Branch("isolation_min_dist_trk_shr", &vars.m_isolation_min_dist_trk_shr);
    vars.vertex_tree->Branch("isolation_min_dist_trk_unassoc", &vars.m_isolation_min_dist_trk_unassoc);

    vars.vertex_tree->Branch("isolation_num_shr_hits_win_1cm_trk",  &vars.m_isolation_num_shr_hits_win_1cm_trk);
    vars.vertex_tree->Branch("isolation_num_shr_hits_win_2cm_trk",  &vars.m_isolation_num_shr_hits_win_2cm_trk);
    vars.vertex_tree->Branch("isolation_num_shr_hits_win_5cm_trk",  &vars.m_isolation_num_shr_hits_win_5cm_trk);
    vars.vertex_tree->Branch("isolation_num_shr_hits_win_10cm_trk", &vars.m_isolation_num_shr_hits_win_10cm_trk);

    vars.vertex_tree->Branch("isolation_num_unassoc_hits_win_1cm_trk",  &vars.m_isolation_num_unassoc_hits_win_1cm_trk);
    vars.vertex_tree->Branch("isolation_num_unassoc_hits_win_2cm_trk",  &vars.m_isolation_num_unassoc_hits_win_2cm_trk);
    vars.vertex_tree->Branch("isolation_num_unassoc_hits_win_5cm_trk",  &vars.m_isolation_num_unassoc_hits_win_5cm_trk);
    vars.vertex_tree->Branch("isolation_num_unassoc_hits_win_10cm_trk", &vars.m_isolation_num_unassoc_hits_win_10cm_trk);


    vars.vertex_tree->Branch("isolation_nearest_shr_hit_to_trk_wire", &vars.m_isolation_nearest_shr_hit_to_trk_wire);
    vars.vertex_tree->Branch("isolation_nearest_shr_hit_to_trk_time", &vars.m_isolation_nearest_shr_hit_to_trk_time);

    vars.vertex_tree->Branch("isolation_nearest_unassoc_hit_to_trk_wire", &vars.m_isolation_nearest_unassoc_hit_to_trk_wire);
    vars.vertex_tree->Branch("isolation_nearest_unassoc_hit_to_trk_time", &vars.m_isolation_nearest_unassoc_hit_to_trk_time);

  }

  //second_shower_search.h
  void ClearSecondShowers(var_all& vars){
    vars.m_sss_num_unassociated_hits=0;
    vars.m_sss_num_unassociated_hits_below_threshold=0;
    vars.m_sss_num_associated_hits=0;

    vars.m_sss_num_candidates = 0;

    vars.m_sss_candidate_in_nu_slice.clear();
    vars.m_sss_candidate_num_hits.clear();
    vars.m_sss_candidate_num_wires.clear();
    vars.m_sss_candidate_num_ticks.clear();
    vars.m_sss_candidate_plane.clear();
    vars.m_sss_candidate_PCA.clear();
    vars.m_sss_candidate_mean_ADC.clear();
    vars.m_sss_candidate_ADC_RMS.clear();
    vars.m_sss_candidate_impact_parameter.clear();
    vars.m_sss_candidate_fit_slope.clear();
    vars.m_sss_candidate_veto_score.clear();
    vars.m_sss_candidate_fit_constant.clear();
    vars.m_sss_candidate_mean_tick.clear();
    vars.m_sss_candidate_max_tick.clear();
    vars.m_sss_candidate_min_tick.clear();
    vars.m_sss_candidate_min_wire.clear();
    vars.m_sss_candidate_max_wire.clear();
    vars.m_sss_candidate_mean_wire.clear();
    vars.m_sss_candidate_min_dist.clear();
    vars.m_sss_candidate_wire_tick_based_length.clear();
    vars.m_sss_candidate_energy.clear();
    vars.m_sss_candidate_angle_to_shower.clear();
    vars.m_sss_candidate_closest_neighbour.clear();
    vars.m_sss_candidate_matched.clear();
    vars.m_sss_candidate_matched_energy_fraction_best_plane.clear();
    vars.m_sss_candidate_pdg.clear();
    vars.m_sss_candidate_parent_pdg.clear();
    vars.m_sss_candidate_trackid.clear();
    vars.m_sss_candidate_true_energy.clear();
    vars.m_sss_candidate_overlay_fraction.clear();
    vars.m_sss_candidate_remerge.clear();
  }

  void ClearSecondShowers3D(var_all& vars){

    vars.m_sss3d_num_showers = 0;
    vars.m_sss3d_shower_start_x.clear();
    vars.m_sss3d_shower_start_y.clear();
    vars.m_sss3d_shower_start_z.clear();
    vars.m_sss3d_shower_dir_x.clear();
    vars.m_sss3d_shower_dir_y.clear();
    vars.m_sss3d_shower_dir_z.clear();
    vars.m_sss3d_shower_length.clear();
    vars.m_sss3d_shower_conversion_dist.clear();
    vars.m_sss3d_shower_invariant_mass.clear();
    vars.m_sss3d_shower_implied_invariant_mass.clear();
    vars.m_sss3d_shower_impact_parameter.clear();
    vars.m_sss3d_shower_energy_max.clear();
    vars.m_sss3d_shower_score.clear();
    vars.m_sss3d_slice_nu.clear();
    vars.m_sss3d_slice_clear_cosmic.clear();
    vars.m_sss3d_shower_ioc_ratio.clear();
  }

  void ClearStubs(var_all& vars){
    vars.m_trackstub_num_unassociated_hits = 0; 
    vars.m_trackstub_unassociated_hits_below_threshold = 0; 
    vars.m_trackstub_associated_hits=0; 
    vars.m_trackstub_num_candidates=0; 
    vars.m_trackstub_candidate_in_nu_slice.clear();
    vars.m_trackstub_candidate_num_hits.clear();
    vars.m_trackstub_candidate_num_wires.clear(); 
    vars.m_trackstub_candidate_num_ticks.clear();
    vars.m_trackstub_candidate_plane.clear(); 
    vars.m_trackstub_candidate_PCA.clear();
    vars.m_trackstub_candidate_mean_ADC.clear();
    vars.m_trackstub_candidate_ADC_RMS.clear();
    vars.m_trackstub_candidate_veto_score.clear();
    vars.m_trackstub_candidate_mean_tick.clear();
    vars.m_trackstub_candidate_max_tick.clear();
    vars.m_trackstub_candidate_min_tick.clear();
    vars.m_trackstub_candidate_min_wire.clear();
    vars.m_trackstub_candidate_max_wire.clear();
    vars.m_trackstub_candidate_mean_wire.clear();
    vars.m_trackstub_candidate_min_dist.clear();  
    vars.m_trackstub_candidate_min_impact_parameter_to_shower.clear(); 
    vars.m_trackstub_candidate_min_conversion_dist_to_shower_start.clear();  
    vars.m_trackstub_candidate_min_ioc_to_shower_start.clear();        
    vars.m_trackstub_candidate_ioc_based_length.clear();         
    vars.m_trackstub_candidate_wire_tick_based_length.clear();           
    vars.m_trackstub_candidate_mean_ADC_first_half.clear();              
    vars.m_trackstub_candidate_mean_ADC_second_half.clear();
    vars.m_trackstub_candidate_mean_ADC_first_to_second_ratio.clear(); 
    vars.m_trackstub_candidate_track_angle_wrt_shower_direction.clear();   
    vars.m_trackstub_candidate_linear_fit_chi2.clear();          
    vars.m_trackstub_candidate_energy.clear();
    vars.m_trackstub_candidate_remerge.clear(); 
    vars.m_trackstub_candidate_matched.clear(); 
    vars.m_trackstub_candidate_matched_energy_fraction_best_plane.clear(); 
    vars.m_trackstub_candidate_pdg.clear();   
    vars.m_trackstub_candidate_parent_pdg.clear();
    vars.m_trackstub_candidate_trackid.clear(); 
    vars.m_trackstub_candidate_true_energy.clear();
    vars.m_trackstub_candidate_overlay_fraction.clear(); 

    vars.m_trackstub_num_candidate_groups = 0;                
    vars.m_grouped_trackstub_candidate_indices.clear(); 
    vars.m_trackstub_candidate_group_timeoverlap_fraction.clear();   
  }

  void CreateSecondShowerBranches(var_all& vars){
    vars.vertex_tree->Branch("sss_num_unassociated_hits",&vars.m_sss_num_unassociated_hits,"sss_num_unassociated_hits/I");
    vars.vertex_tree->Branch("sss_num_unassociated_hits_below_threshold",&vars.m_sss_num_unassociated_hits_below_threshold,"sss_num_unassociated_hits_below_threshold/I");
    vars.vertex_tree->Branch("sss_num_associated_hits",&vars.m_sss_num_associated_hits,"sss_num_associated_hits/I");

    vars.vertex_tree->Branch("sss_num_candidates",&vars.m_sss_num_candidates,"sss_num_candidates/I");
    vars.vertex_tree->Branch("sss_candidate_veto_score",&vars.m_sss_candidate_veto_score);
    vars.vertex_tree->Branch("sss_candidate_in_nu_slice", &vars.m_sss_candidate_in_nu_slice);
    vars.vertex_tree->Branch("sss_candidate_num_hits",&vars.m_sss_candidate_num_hits);
    vars.vertex_tree->Branch("sss_candidate_num_wires",&vars.m_sss_candidate_num_wires);
    vars.vertex_tree->Branch("sss_candidate_num_ticks",&vars.m_sss_candidate_num_ticks);
    vars.vertex_tree->Branch("sss_candidate_plane",&vars.m_sss_candidate_plane);
    vars.vertex_tree->Branch("sss_candidate_PCA",&vars.m_sss_candidate_PCA);
    vars.vertex_tree->Branch("sss_candidate_mean_ADC",&vars.m_sss_candidate_mean_ADC);
    vars.vertex_tree->Branch("sss_candidate_ADC_RMS", &vars.m_sss_candidate_ADC_RMS);
    vars.vertex_tree->Branch("sss_candidate_impact_parameter",&vars.m_sss_candidate_impact_parameter); 
    vars.vertex_tree->Branch("sss_candidate_fit_slope",&vars.m_sss_candidate_fit_slope);
    vars.vertex_tree->Branch("sss_candidate_fit_constant",&vars.m_sss_candidate_fit_constant);
    vars.vertex_tree->Branch("sss_candidate_mean_tick",&vars.m_sss_candidate_mean_tick);
    vars.vertex_tree->Branch("sss_candidate_max_tick",&vars.m_sss_candidate_max_tick);
    vars.vertex_tree->Branch("sss_candidate_min_tick",&vars.m_sss_candidate_min_tick);
    vars.vertex_tree->Branch("sss_candidate_mean_wire",&vars.m_sss_candidate_mean_wire);
    vars.vertex_tree->Branch("sss_candidate_max_wire",&vars.m_sss_candidate_max_wire);
    vars.vertex_tree->Branch("sss_candidate_min_wire",&vars.m_sss_candidate_min_wire);
    vars.vertex_tree->Branch("sss_candidate_min_dist",&vars.m_sss_candidate_min_dist);
    vars.vertex_tree->Branch("sss_candidate_wire_tick_based_length", &vars.m_sss_candidate_wire_tick_based_length);
    vars.vertex_tree->Branch("sss_candidate_energy",&vars.m_sss_candidate_energy);
    vars.vertex_tree->Branch("sss_candidate_angle_to_shower",&vars.m_sss_candidate_angle_to_shower);
    vars.vertex_tree->Branch("sss_candidate_closest_neighbour",&vars.m_sss_candidate_closest_neighbour);
    vars.vertex_tree->Branch("sss_candidate_remerge",&vars.m_sss_candidate_remerge);

    vars.vertex_tree->Branch("sss_candidate_matched",&vars.m_sss_candidate_matched);
    vars.vertex_tree->Branch("sss_candidate_pdg",&vars.m_sss_candidate_pdg);
    vars.vertex_tree->Branch("sss_candidate_parent_pdg",&vars.m_sss_candidate_parent_pdg);
    vars.vertex_tree->Branch("sss_candidate_trackid",&vars.m_sss_candidate_trackid);
    vars.vertex_tree->Branch("sss_candidate_true_energy", &vars.m_sss_candidate_true_energy);
    vars.vertex_tree->Branch("sss_candidate_overlay_fraction",&vars.m_sss_candidate_overlay_fraction);
    vars.vertex_tree->Branch("sss_candidate_matched_energy_fraction_best_plane", &vars.m_sss_candidate_matched_energy_fraction_best_plane);


    vars.vertex_tree->Branch("sss3d_ioc_ranked_en",&vars.m_sss3d_ioc_ranked_en);
    vars.vertex_tree->Branch("sss3d_ioc_ranked_conv",&vars.m_sss3d_ioc_ranked_conv);
    vars.vertex_tree->Branch("sss3d_ioc_ranked_invar",&vars.m_sss3d_ioc_ranked_invar);
    vars.vertex_tree->Branch("sss3d_ioc_ranked_implied_invar",&vars.m_sss3d_ioc_ranked_implied_invar);
    vars.vertex_tree->Branch("sss3d_ioc_ranked_ioc",&vars.m_sss3d_ioc_ranked_ioc);
    vars.vertex_tree->Branch("sss3d_ioc_ranked_opang",&vars.m_sss3d_ioc_ranked_opang);
    vars.vertex_tree->Branch("sss3d_ioc_ranked_implied_opang",&vars.m_sss3d_ioc_ranked_implied_opang);
    vars.vertex_tree->Branch("sss3d_ioc_ranked_id",&vars.m_sss3d_ioc_ranked_id);

    vars.vertex_tree->Branch("sss3d_invar_ranked_en",&vars.m_sss3d_invar_ranked_en);
    vars.vertex_tree->Branch("sss3d_invar_ranked_conv",&vars.m_sss3d_invar_ranked_conv);
    vars.vertex_tree->Branch("sss3d_invar_ranked_invar",&vars.m_sss3d_invar_ranked_invar);
    vars.vertex_tree->Branch("sss3d_invar_ranked_implied_invar",&vars.m_sss3d_invar_ranked_implied_invar);
    vars.vertex_tree->Branch("sss3d_invar_ranked_ioc",&vars.m_sss3d_invar_ranked_ioc);
    vars.vertex_tree->Branch("sss3d_invar_ranked_opang",&vars.m_sss3d_invar_ranked_opang);
    vars.vertex_tree->Branch("sss3d_invar_ranked_implied_opang",&vars.m_sss3d_invar_ranked_implied_opang);
    vars.vertex_tree->Branch("sss3d_invar_ranked_id",&vars.m_sss3d_invar_ranked_id);


    vars.vertex_tree->Branch("sss2d_ioc_ranked_en",&vars.m_sss2d_ioc_ranked_en);
    vars.vertex_tree->Branch("sss2d_ioc_ranked_conv",&vars.m_sss2d_ioc_ranked_conv);
    vars.vertex_tree->Branch("sss2d_ioc_ranked_ioc",&vars.m_sss2d_ioc_ranked_ioc);
    vars.vertex_tree->Branch("sss2d_ioc_ranked_pca",&vars.m_sss2d_ioc_ranked_pca);
    vars.vertex_tree->Branch("sss2d_ioc_ranked_invar",&vars.m_sss2d_ioc_ranked_invar);
    vars.vertex_tree->Branch("sss2d_ioc_ranked_angle_to_shower",&vars.m_sss2d_ioc_ranked_angle_to_shower);
    vars.vertex_tree->Branch("sss2d_ioc_ranked_num_planes",&vars.m_sss2d_ioc_ranked_num_planes);

    vars.vertex_tree->Branch("sss2d_invar_ranked_en",&vars.m_sss2d_invar_ranked_en);
    vars.vertex_tree->Branch("sss2d_invar_ranked_conv",&vars.m_sss2d_invar_ranked_conv);
    vars.vertex_tree->Branch("sss2d_invar_ranked_ioc",&vars.m_sss2d_invar_ranked_ioc);
    vars.vertex_tree->Branch("sss2d_invar_ranked_pca",&vars.m_sss2d_invar_ranked_pca);
    vars.vertex_tree->Branch("sss2d_invar_ranked_invar",&vars.m_sss2d_invar_ranked_invar);
    vars.vertex_tree->Branch("sss2d_invar_ranked_angle_to_shower",&vars.m_sss2d_invar_ranked_angle_to_shower);
    vars.vertex_tree->Branch("sss2d_invar_ranked_num_planes",&vars.m_sss2d_invar_ranked_num_planes);

    vars.vertex_tree->Branch("sss2d_conv_ranked_en",&vars.m_sss2d_conv_ranked_en);
    vars.vertex_tree->Branch("sss2d_conv_ranked_conv",&vars.m_sss2d_conv_ranked_conv);
    vars.vertex_tree->Branch("sss2d_conv_ranked_ioc",&vars.m_sss2d_conv_ranked_ioc);
    vars.vertex_tree->Branch("sss2d_conv_ranked_pca",&vars.m_sss2d_conv_ranked_pca);
    vars.vertex_tree->Branch("sss2d_conv_ranked_invar",&vars.m_sss2d_conv_ranked_invar);
    vars.vertex_tree->Branch("sss2d_conv_ranked_angle_to_shower",&vars.m_sss2d_conv_ranked_angle_to_shower);
    vars.vertex_tree->Branch("sss2d_conv_ranked_num_planes",&vars.m_sss2d_conv_ranked_num_planes);

  }

  void CreateSecondShowerBranches3D(var_all& vars){
    vars.vertex_tree->Branch("sss3d_num_showers",&vars.m_sss3d_num_showers,"sss3d_num_showers/I");

    vars.vertex_tree->Branch("sss3d_shower_start_x",&vars.m_sss3d_shower_start_x);
    vars.vertex_tree->Branch("sss3d_shower_start_y",&vars.m_sss3d_shower_start_y);
    vars.vertex_tree->Branch("sss3d_shower_start_z",&vars.m_sss3d_shower_start_z);
    vars.vertex_tree->Branch("sss3d_shower_dir_x",&vars.m_sss3d_shower_dir_x);
    vars.vertex_tree->Branch("sss3d_shower_dir_y",&vars.m_sss3d_shower_dir_y);
    vars.vertex_tree->Branch("sss3d_shower_dir_z",&vars.m_sss3d_shower_dir_z);

    vars.vertex_tree->Branch("sss3d_shower_length",&vars.m_sss3d_shower_length);
    vars.vertex_tree->Branch("sss3d_shower_conversion_dist",&vars.m_sss3d_shower_conversion_dist);
    vars.vertex_tree->Branch("sss3d_shower_invariant_mass",&vars.m_sss3d_shower_invariant_mass);
    vars.vertex_tree->Branch("sss3d_shower_implied_invariant_mass",&vars.m_sss3d_shower_implied_invariant_mass);
    vars.vertex_tree->Branch("sss3d_shower_impact_parameter",&vars.m_sss3d_shower_impact_parameter);
    vars.vertex_tree->Branch("sss3d_shower_ioc_ratio",&vars.m_sss3d_shower_ioc_ratio);
    vars.vertex_tree->Branch("sss3d_shower_energy_max",&vars.m_sss3d_shower_energy_max);
    vars.vertex_tree->Branch("sss3d_shower_score",&vars.m_sss3d_shower_score);
    vars.vertex_tree->Branch("sss3d_slice_nu",&vars.m_sss3d_slice_nu);
    vars.vertex_tree->Branch("sss3d_slice_clear_cosmic",&vars.m_sss3d_slice_clear_cosmic);
  }

  void CreateStubBranches(var_all& vars){

    vars.vertex_tree->Branch("trackstub_num_unassociated_hits",&vars.m_trackstub_num_unassociated_hits,"trackstub_num_unassociated_hits/I");
    vars.vertex_tree->Branch("trackstub_unassociated_hits_below_threshold",&vars.m_trackstub_unassociated_hits_below_threshold,"trackstub_unassociated_hits_below_threshold/I");
    vars.vertex_tree->Branch("trackstub_associated_hits",&vars.m_trackstub_associated_hits,"trackstub_associated_hits/I");
    vars.vertex_tree->Branch("trackstub_num_candidates", &vars.m_trackstub_num_candidates, "trackstub_num_candidates/I");
    vars.vertex_tree->Branch("trackstub_candidate_in_nu_slice", &vars.m_trackstub_candidate_in_nu_slice);
    vars.vertex_tree->Branch("trackstub_candidate_num_hits", &vars.m_trackstub_candidate_num_hits);
    vars.vertex_tree->Branch("trackstub_candidate_num_wires", &vars.m_trackstub_candidate_num_wires);
    vars.vertex_tree->Branch("trackstub_candidate_num_ticks", &vars.m_trackstub_candidate_num_ticks);
    vars.vertex_tree->Branch("trackstub_candidate_plane", &vars.m_trackstub_candidate_plane);
    vars.vertex_tree->Branch("trackstub_candidate_PCA", &vars.m_trackstub_candidate_PCA);
    vars.vertex_tree->Branch("trackstub_candidate_mean_ADC", &vars.m_trackstub_candidate_mean_ADC);
    vars.vertex_tree->Branch("trackstub_candidate_ADC_RMS", &vars.m_trackstub_candidate_ADC_RMS);
    vars.vertex_tree->Branch("trackstub_candidate_veto_score", &vars.m_trackstub_candidate_veto_score);
    vars.vertex_tree->Branch("trackstub_candidate_mean_tick", &vars.m_trackstub_candidate_mean_tick);
    vars.vertex_tree->Branch("trackstub_candidate_max_tick", &vars.m_trackstub_candidate_max_tick);
    vars.vertex_tree->Branch("trackstub_candidate_min_tick", &vars.m_trackstub_candidate_min_tick);
    vars.vertex_tree->Branch("trackstub_candidate_min_wire", &vars.m_trackstub_candidate_min_wire);
    vars.vertex_tree->Branch("trackstub_candidate_max_wire", &vars.m_trackstub_candidate_max_wire);
    vars.vertex_tree->Branch("trackstub_candidate_mean_wire", &vars.m_trackstub_candidate_mean_wire);
    vars.vertex_tree->Branch("trackstub_candidate_min_dist", &vars.m_trackstub_candidate_min_dist);
    vars.vertex_tree->Branch("trackstub_candidate_min_impact_parameter_to_shower", &vars.m_trackstub_candidate_min_impact_parameter_to_shower);
    vars.vertex_tree->Branch("trackstub_candidate_min_conversion_dist_to_shower_start", &vars.m_trackstub_candidate_min_conversion_dist_to_shower_start);
    vars.vertex_tree->Branch("trackstub_candidate_min_ioc_to_shower_start", &vars.m_trackstub_candidate_min_ioc_to_shower_start);
    vars.vertex_tree->Branch("trackstub_candidate_ioc_based_length", &vars.m_trackstub_candidate_ioc_based_length);
    vars.vertex_tree->Branch("trackstub_candidate_wire_tick_based_length", &vars.m_trackstub_candidate_wire_tick_based_length);
    vars.vertex_tree->Branch("trackstub_candidate_mean_ADC_first_half", &vars.m_trackstub_candidate_mean_ADC_first_half);
    vars.vertex_tree->Branch("trackstub_candidate_mean_ADC_second_half", &vars.m_trackstub_candidate_mean_ADC_second_half);
    vars.vertex_tree->Branch("trackstub_candidate_mean_ADC_first_to_second_ratio", &vars.m_trackstub_candidate_mean_ADC_first_to_second_ratio);
    vars.vertex_tree->Branch("trackstub_candidate_track_angle_wrt_shower_direction", &vars.m_trackstub_candidate_track_angle_wrt_shower_direction);
    vars.vertex_tree->Branch("trackstub_candidate_linear_fit_chi2", &vars.m_trackstub_candidate_linear_fit_chi2);
    vars.vertex_tree->Branch("trackstub_candidate_energy", &vars.m_trackstub_candidate_energy);
    vars.vertex_tree->Branch("trackstub_candidate_remerge", &vars.m_trackstub_candidate_remerge);
    vars.vertex_tree->Branch("trackstub_candidate_matched", &vars.m_trackstub_candidate_matched);
    vars.vertex_tree->Branch("trackstub_candidate_matched_energy_fraction_best_plane", &vars.m_trackstub_candidate_matched_energy_fraction_best_plane);
    vars.vertex_tree->Branch("trackstub_candidate_pdg", &vars.m_trackstub_candidate_pdg);
    vars.vertex_tree->Branch("trackstub_candidate_parent_pdg", &vars.m_trackstub_candidate_parent_pdg);
    vars.vertex_tree->Branch("trackstub_candidate_trackid", &vars.m_trackstub_candidate_trackid);
    vars.vertex_tree->Branch("trackstub_candidate_true_energy", &vars.m_trackstub_candidate_true_energy);
    vars.vertex_tree->Branch("trackstub_candidate_overlay_fraction", &vars.m_trackstub_candidate_overlay_fraction);


    vars.vertex_tree->Branch("trackstub_num_candidate_groups", &vars.m_trackstub_num_candidate_groups, "trackstub_num_candidate_groups/I");
    vars.vertex_tree->Branch("grouped_trackstub_candidate_indices", &vars.m_grouped_trackstub_candidate_indices);
    vars.vertex_tree->Branch("trackstub_candidate_group_timeoverlap_fraction", &vars.m_trackstub_candidate_group_timeoverlap_fraction);

  }

//  void ResizeSecondShowers(size_t size){}

  //analyze_OpFlashes.h
  void ClearFlashes(var_all& vars){
    vars.m_reco_num_flashes =0;
    vars.m_reco_num_flashes_in_beamgate =0;
    vars.m_reco_flash_total_pe.clear();
    vars.m_reco_flash_time.clear();
    vars.m_reco_flash_time_width.clear();
    vars.m_reco_flash_abs_time.clear();
    vars.m_reco_flash_frame.clear();
    vars.m_reco_flash_ycenter.clear();
    vars.m_reco_flash_ywidth.clear();
    vars.m_reco_flash_zcenter.clear();
    vars.m_reco_flash_zwidth.clear();
    vars.m_reco_flash_total_pe_in_beamgate.clear();
    vars.m_reco_flash_time_in_beamgate.clear();
    vars.m_reco_flash_ycenter_in_beamgate.clear();
    vars.m_reco_flash_zcenter_in_beamgate.clear();
    vars.m_CRT_veto_nhits = -999;
    vars.m_CRT_veto_hit_PE.clear();
    vars.m_CRT_min_hit_time = -999;
    vars.m_CRT_min_hit_PE = -999;
    vars.m_CRT_min_hit_x = -999;
    vars.m_CRT_min_hit_y = -999;
    vars.m_CRT_min_hit_z = -999;
    vars.m_CRT_hits_time.clear();
    vars.m_CRT_hits_PE.clear();
    vars.m_CRT_hits_x.clear(); 
    vars.m_CRT_hits_y.clear();
    vars.m_CRT_hits_z.clear();
    vars.m_CRT_dt = -999;

  }

  void CreateFlashBranches(var_all& vars){
//    vars.vertex_tree->Branch("beamgate_flash_start",&vars.m_beamgate_flash_start,"beamgate_flash_start/D");
//    vars.vertex_tree->Branch("beamgate_flash_end",&vars.m_beamgate_flash_end,"beamgate_flash_end/D");
    vars.vertex_tree->Branch("reco_num_flashes",&vars.m_reco_num_flashes,"reco_num_flashes/I");
    vars.vertex_tree->Branch("reco_num_flashes_in_beamgate",&vars.m_reco_num_flashes_in_beamgate,"reco_num_flashes_in_beamgate/I");
    vars.vertex_tree->Branch("reco_flash_total_pe", &vars.m_reco_flash_total_pe);
    vars.vertex_tree->Branch("reco_flash_time", &vars.m_reco_flash_time);
    vars.vertex_tree->Branch("reco_flash_time_width",&vars.m_reco_flash_time_width);
    vars.vertex_tree->Branch("reco_flash_abs_time",&vars.m_reco_flash_abs_time);
    vars.vertex_tree->Branch("reco_flash_frame",&vars.m_reco_flash_frame);
    vars.vertex_tree->Branch("reco_flash_ycenter",&vars.m_reco_flash_ycenter);
    vars.vertex_tree->Branch("reco_flash_ywidth",&vars.m_reco_flash_ywidth);
    vars.vertex_tree->Branch("reco_flash_zcenter",&vars.m_reco_flash_zcenter);
    vars.vertex_tree->Branch("reco_flash_zwidth",&vars.m_reco_flash_zwidth);
    vars.vertex_tree->Branch("reco_flash_total_pe_in_beamgate", &vars.m_reco_flash_total_pe_in_beamgate);
    vars.vertex_tree->Branch("reco_flash_time_in_beamgate", &vars.m_reco_flash_time_in_beamgate);
    vars.vertex_tree->Branch("reco_flash_ycenter_in_beamgate",&vars.m_reco_flash_ycenter_in_beamgate);
    vars.vertex_tree->Branch("reco_flash_zcenter_in_beamgate",&vars.m_reco_flash_zcenter_in_beamgate);

    vars.vertex_tree->Branch("CRT_veto_nhits",&vars.m_CRT_veto_nhits,"CRT_veto_nhits/I");
    vars.vertex_tree->Branch("CRT_veto_hit_PE",&vars.m_CRT_veto_hit_PE);
    vars.vertex_tree->Branch("CRT_dt",& vars.m_CRT_dt," CRT_dt/D");
    vars.vertex_tree->Branch("CRT_min_hit_time",&vars.m_CRT_min_hit_time,"CRT_min_hit_time/D");
    vars.vertex_tree->Branch("CRT_min_hit_PE",&vars.m_CRT_min_hit_PE,"CRT_min_hit_PE/D");
    vars.vertex_tree->Branch("CRT_min_hit_x",&vars.m_CRT_min_hit_x,"CRT_min_hit_x/D");
    vars.vertex_tree->Branch("CRT_min_hit_y",&vars.m_CRT_min_hit_y,"CRT_min_hit_y/D");
    vars.vertex_tree->Branch("CRT_min_hit_z",&vars.m_CRT_min_hit_z,"CRT_min_hit_z/D");
    vars.vertex_tree->Branch("CRT_hits_time",&vars.m_CRT_hits_time);
    vars.vertex_tree->Branch("CRT_hits_PE",&vars.m_CRT_hits_PE);
    vars.vertex_tree->Branch("CRT_hits_x",&vars.m_CRT_hits_x);
    vars.vertex_tree->Branch("CRT_hits_y",&vars.m_CRT_hits_y);
    vars.vertex_tree->Branch("CRT_hits_z",&vars.m_CRT_hits_z);
  }

  void ResizeFlashes(size_t size, var_all& vars){
    vars.m_reco_flash_total_pe.resize(size);
    vars.m_reco_flash_time.resize(size);
    vars.m_reco_flash_time_width.resize(size);
    vars.m_reco_flash_abs_time.resize(size);
    vars.m_reco_flash_frame.resize(size);
    vars.m_reco_flash_ycenter.resize(size);
    vars.m_reco_flash_ywidth.resize(size);
    vars.m_reco_flash_zcenter.resize(size);
    vars.m_reco_flash_zwidth.resize(size);
    vars.m_reco_flash_total_pe_in_beamgate.resize(size);
    vars.m_reco_flash_time_in_beamgate.resize(size);
    vars.m_reco_flash_ycenter_in_beamgate.resize(size);
    vars.m_reco_flash_zcenter_in_beamgate.resize(size);
    vars.m_CRT_veto_hit_PE.resize(size);

    vars.m_CRT_hits_time.resize(size);
    vars.m_CRT_hits_PE.resize(size);
    vars.m_CRT_hits_x.resize(size); 
    vars.m_CRT_hits_y.resize(size);
    vars.m_CRT_hits_z.resize(size);
  }

  //analyze_Tracks.h
  void ClearTracks(var_all& vars){
    vars.m_reco_asso_tracks=0;
    vars.m_reco_track_length.clear();
    vars.m_reco_track_num_daughters.clear();
    vars.m_reco_track_daughter_trackscore.clear();
    vars.m_reco_track_dirx.clear();
    vars.m_reco_track_diry.clear();
    vars.m_reco_track_dirz.clear();
    vars.m_reco_track_startx.clear();
    vars.m_reco_track_starty.clear();
    vars.m_reco_track_startz.clear();
    vars.m_reco_track_endx.clear();
    vars.m_reco_track_endy.clear();
    vars.m_reco_track_endz.clear();
    vars.m_reco_track_end_dist_to_active_TPC.clear();
    vars.m_reco_track_start_dist_to_active_TPC.clear();
    vars.m_reco_track_end_dist_to_CPA.clear();
    vars.m_reco_track_start_dist_to_CPA.clear();
    vars.m_reco_track_end_dist_to_SCB.clear();
    vars.m_reco_track_start_dist_to_SCB.clear();
    vars.m_reco_track_end_in_SCB.clear();
    vars.m_reco_track_start_in_SCB.clear();

    vars.m_reco_track_theta_yz.clear();
    vars.m_reco_track_phi_yx.clear();

    vars.m_reco_track_calo_energy_plane0.clear();
    vars.m_reco_track_calo_energy_plane1.clear();
    vars.m_reco_track_calo_energy_plane2.clear();
    vars.m_reco_track_calo_energy_max.clear();

    vars.m_reco_track_num_trajpoints.clear();
    vars.m_reco_track_num_spacepoints.clear();
    vars.m_reco_track_proton_kinetic_energy.clear();
    vars.m_reco_track_ordered_energy_index.clear();
    vars.m_reco_track_ordered_displacement_index.clear();

    vars.m_reco_track_spacepoint_principal0.clear();
    vars.m_reco_track_spacepoint_principal1.clear();
    vars.m_reco_track_spacepoint_principal2.clear();

    vars.m_reco_track_spacepoint_chi.clear();
    vars.m_reco_track_spacepoint_max_dist.clear();

    vars.m_reco_track_best_calo_plane.clear();

    vars.m_reco_track_mean_dEdx_best_plane.clear();
    vars.m_reco_track_mean_dEdx_start_half_best_plane.clear();
    vars.m_reco_track_mean_dEdx_end_half_best_plane.clear();
    vars.m_reco_track_good_calo_best_plane.clear();
    vars.m_reco_track_trunc_dEdx_best_plane.clear();
    vars.m_reco_track_mean_trunc_dEdx_best_plane.clear();
    vars.m_reco_track_mean_trunc_dEdx_start_half_best_plane.clear();
    vars.m_reco_track_mean_trunc_dEdx_end_half_best_plane.clear();
    vars.m_reco_track_trunc_PIDA_best_plane.clear();
    vars.m_reco_track_resrange_best_plane.clear();
    vars.m_reco_track_dEdx_best_plane.clear();


    vars.m_reco_track_mean_dEdx_p0.clear();
    vars.m_reco_track_mean_dEdx_start_half_p0.clear();
    vars.m_reco_track_mean_dEdx_end_half_p0.clear();
    vars.m_reco_track_good_calo_p0.clear();
    vars.m_reco_track_trunc_dEdx_p0.clear();
    vars.m_reco_track_mean_trunc_dEdx_p0.clear();
    vars.m_reco_track_mean_trunc_dEdx_start_half_p0.clear();
    vars.m_reco_track_mean_trunc_dEdx_end_half_p0.clear();
    vars.m_reco_track_trunc_PIDA_p0.clear();
    vars.m_reco_track_resrange_p0.clear();
    vars.m_reco_track_dEdx_p0.clear();

    vars.m_reco_track_mean_dEdx_p1.clear();
    vars.m_reco_track_mean_dEdx_start_half_p1.clear();
    vars.m_reco_track_mean_dEdx_end_half_p1.clear();
    vars.m_reco_track_good_calo_p1.clear();
    vars.m_reco_track_trunc_dEdx_p1.clear();
    vars.m_reco_track_mean_trunc_dEdx_p1.clear();
    vars.m_reco_track_mean_trunc_dEdx_start_half_p1.clear();
    vars.m_reco_track_mean_trunc_dEdx_end_half_p1.clear();
    vars.m_reco_track_trunc_PIDA_p1.clear();
    vars.m_reco_track_resrange_p1.clear();
    vars.m_reco_track_dEdx_p1.clear();

    vars.m_reco_track_mean_dEdx_p2.clear();
    vars.m_reco_track_mean_dEdx_start_half_p2.clear();
    vars.m_reco_track_mean_dEdx_end_half_p2.clear();
    vars.m_reco_track_good_calo_p2.clear();
    vars.m_reco_track_trunc_dEdx_p2.clear();
    vars.m_reco_track_mean_trunc_dEdx_p2.clear();
    vars.m_reco_track_mean_trunc_dEdx_start_half_p2.clear();
    vars.m_reco_track_mean_trunc_dEdx_end_half_p2.clear();
    vars.m_reco_track_trunc_PIDA_p2.clear();
    vars.m_reco_track_resrange_p2.clear();
    vars.m_reco_track_dEdx_p2.clear();

    vars.m_reco_track_num_calo_hits_p1.clear();
    vars.m_reco_track_num_calo_hits_p0.clear();
    vars.m_reco_track_num_calo_hits_p2.clear();

    vars.m_sim_track_matched.clear();
    vars.m_sim_track_overlay_fraction.clear();
    vars.m_sim_track_energy.clear();
    vars.m_sim_track_kinetic_energy.clear();
    vars.m_sim_track_mass.clear();
    vars.m_sim_track_pdg.clear();
    vars.m_sim_track_origin.clear();
    vars.m_sim_track_parent_pdg.clear();
    vars.m_sim_track_process.clear();
    vars.m_sim_track_startx.clear();
    vars.m_sim_track_starty.clear();
    vars.m_sim_track_startz.clear();
    vars.m_sim_track_endx.clear();
    vars.m_sim_track_endy.clear();
    vars.m_sim_track_endz.clear();
    vars.m_sim_track_length.clear();

    vars.m_sim_track_px.clear();
    vars.m_sim_track_py.clear();
    vars.m_sim_track_pz.clear();
    vars.m_sim_track_trackID.clear();

    // PID
    vars.m_reco_track_pid_bragg_likelihood_mu_plane0.clear();
    vars.m_reco_track_pid_bragg_likelihood_mu_plane1.clear();
    vars.m_reco_track_pid_bragg_likelihood_mu_plane2.clear();
    vars.m_reco_track_pid_bragg_likelihood_p_plane0.clear();
    vars.m_reco_track_pid_bragg_likelihood_p_plane1.clear();
    vars.m_reco_track_pid_bragg_likelihood_p_plane2.clear();
    vars.m_reco_track_pid_bragg_likelihood_mip_plane0.clear();
    vars.m_reco_track_pid_bragg_likelihood_mip_plane1.clear();
    vars.m_reco_track_pid_bragg_likelihood_mip_plane2.clear();
    vars.m_reco_track_pid_chi2_mu_plane0.clear();
    vars.m_reco_track_pid_chi2_mu_plane1.clear();
    vars.m_reco_track_pid_chi2_mu_plane2.clear();
    vars.m_reco_track_pid_chi2_p_plane0.clear();
    vars.m_reco_track_pid_chi2_p_plane1.clear();
    vars.m_reco_track_pid_chi2_p_plane2.clear();
    vars.m_reco_track_pid_pida_plane0.clear();
    vars.m_reco_track_pid_pida_plane1.clear();
    vars.m_reco_track_pid_pida_plane2.clear();
    vars.m_reco_track_pid_three_plane_proton_pid.clear();

//    vars.m_reco_track_end_to_nearest_dead_wire_plane0.clear();
//    vars.m_reco_track_end_to_nearest_dead_wire_plane1.clear();
//    vars.m_reco_track_end_to_nearest_dead_wire_plane2.clear();

    vars.m_reco_track_sliceId.clear();
    vars.m_reco_track_nuscore.clear();
    vars.m_reco_track_isclearcosmic.clear();
    vars.m_reco_track_trackscore.clear();
    vars.m_reco_track_pfparticle_pdg.clear();
    vars.m_reco_track_is_nuslice.clear();

    vars.m_sim_track_sliceId.clear();
    vars.m_sim_track_nuscore.clear();
    vars.m_sim_track_isclearcosmic.clear();
  }

  void CreateTrackBranches(var_all& vars){
    vars.vertex_tree->Branch("reco_asso_tracks",&vars.m_reco_asso_tracks,"reco_asso_tracks/I");
    vars.vertex_tree->Branch("reco_track_num_daughters",&vars.m_reco_track_num_daughters);
    vars.vertex_tree->Branch("reco_track_daughter_trackscore",&vars.m_reco_track_daughter_trackscore);
    vars.vertex_tree->Branch("reco_track_displacement", &vars.m_reco_track_length);
    vars.vertex_tree->Branch("reco_track_dirx", &vars.m_reco_track_dirx);
    vars.vertex_tree->Branch("reco_track_diry", &vars.m_reco_track_diry);
    vars.vertex_tree->Branch("reco_track_dirz", &vars.m_reco_track_dirz);
    vars.vertex_tree->Branch("reco_track_startx", &vars.m_reco_track_startx);
    vars.vertex_tree->Branch("reco_track_starty", &vars.m_reco_track_starty);
    vars.vertex_tree->Branch("reco_track_startz", &vars.m_reco_track_startz);

    vars.vertex_tree->Branch("reco_track_endx", &vars.m_reco_track_endx);
    vars.vertex_tree->Branch("reco_track_endy", &vars.m_reco_track_endy);
    vars.vertex_tree->Branch("reco_track_endz", &vars.m_reco_track_endz);
    vars.vertex_tree->Branch("reco_track_end_dist_to_active_TPC", &vars.m_reco_track_end_dist_to_active_TPC);
    vars.vertex_tree->Branch("reco_track_start_dist_to_active_TPC", &vars.m_reco_track_start_dist_to_active_TPC);
    vars.vertex_tree->Branch("reco_track_end_dist_to_CPA", &vars.m_reco_track_end_dist_to_CPA);
    vars.vertex_tree->Branch("reco_track_start_dist_to_CPA", &vars.m_reco_track_start_dist_to_CPA);
    vars.vertex_tree->Branch("reco_track_end_dist_to_SCB", &vars.m_reco_track_end_dist_to_SCB);
    vars.vertex_tree->Branch("reco_track_start_dist_to_SCB", &vars.m_reco_track_start_dist_to_SCB);
    vars.vertex_tree->Branch("reco_track_end_in_SCB", &vars.m_reco_track_end_in_SCB);
    vars.vertex_tree->Branch("reco_track_start_in_SCB", &vars.m_reco_track_start_in_SCB);


    vars.vertex_tree->Branch("reco_track_theta_yz", &vars.m_reco_track_theta_yz);
    vars.vertex_tree->Branch("reco_track_phi_yx", &vars.m_reco_track_phi_yx);

    vars.vertex_tree->Branch("reco_track_calo_energy_plane0", &vars.m_reco_track_calo_energy_plane0);
    vars.vertex_tree->Branch("reco_track_calo_energy_plane1", &vars.m_reco_track_calo_energy_plane1);
    vars.vertex_tree->Branch("reco_track_calo_energy_plane2", &vars.m_reco_track_calo_energy_plane2);
    vars.vertex_tree->Branch("reco_track_calo_energy_max", &vars.m_reco_track_calo_energy_max);

    vars.vertex_tree->Branch("reco_track_num_trajpoints", &vars.m_reco_track_num_trajpoints);
    vars.vertex_tree->Branch("reco_track_num_spacepoints", &vars.m_reco_track_num_spacepoints);
    vars.vertex_tree->Branch("reco_track_proton_kinetic_energy", &vars.m_reco_track_proton_kinetic_energy);
    vars.vertex_tree->Branch("reco_track_ordered_energy_index", &vars.m_reco_track_ordered_energy_index);
    vars.vertex_tree->Branch("reco_track_ordered_displacement_index", &vars.m_reco_track_ordered_displacement_index);
    vars.vertex_tree->Branch("i_trk", &vars.m_reco_track_ordered_displacement_index);

    vars.vertex_tree->Branch("reco_track_spacepoint_principal0",&vars.m_reco_track_spacepoint_principal0);
    vars.vertex_tree->Branch("reco_track_spacepoint_principal1",&vars.m_reco_track_spacepoint_principal1);
    vars.vertex_tree->Branch("reco_track_spacepoint_principal2",&vars.m_reco_track_spacepoint_principal2);

    vars.vertex_tree->Branch("reco_track_spacepoint_chi",&vars.m_reco_track_spacepoint_chi);
    vars.vertex_tree->Branch("reco_track_spacepoint_max_dist",&vars.m_reco_track_spacepoint_max_dist);

    vars.vertex_tree->Branch("reco_track_best_calo_plane",&vars.m_reco_track_best_calo_plane);

    vars.vertex_tree->Branch("reco_track_mean_dEdx_best_plane",&vars.m_reco_track_mean_dEdx_best_plane);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_plane0",&vars.m_reco_track_mean_dEdx_p0);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_plane1",&vars.m_reco_track_mean_dEdx_p1);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_plane2",&vars.m_reco_track_mean_dEdx_p2);

    vars.vertex_tree->Branch("reco_track_mean_dEdx_start_half_best_plane",&vars.m_reco_track_mean_dEdx_end_half_best_plane);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_start_half_plane0",&vars.m_reco_track_mean_dEdx_end_half_p0);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_start_half_plane1",&vars.m_reco_track_mean_dEdx_end_half_p1);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_start_half_plane2",&vars.m_reco_track_mean_dEdx_end_half_p2);

    vars.vertex_tree->Branch("reco_track_mean_dEdx_end_half_best_plane",&vars.m_reco_track_mean_dEdx_start_half_best_plane);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_end_half_plane0",&vars.m_reco_track_mean_dEdx_start_half_p0);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_end_half_plane1",&vars.m_reco_track_mean_dEdx_start_half_p1);
    vars.vertex_tree->Branch("reco_track_mean_dEdx_end_half_plane2",&vars.m_reco_track_mean_dEdx_start_half_p2);

    vars.vertex_tree->Branch("reco_track_good_calo_best_plane",&vars.m_reco_track_good_calo_best_plane);
    vars.vertex_tree->Branch("reco_track_good_calo_plane0",&vars.m_reco_track_good_calo_p0);
    vars.vertex_tree->Branch("reco_track_good_calo_plane1",&vars.m_reco_track_good_calo_p1);
    vars.vertex_tree->Branch("reco_track_good_calo_plane2",&vars.m_reco_track_good_calo_p2);

    vars.vertex_tree->Branch("reco_track_trunc_dEdx_best_plane",&vars.m_reco_track_trunc_dEdx_best_plane);
    vars.vertex_tree->Branch("reco_track_trunc_dEdx_plane0",&vars.m_reco_track_trunc_dEdx_p0);
    vars.vertex_tree->Branch("reco_track_trunc_dEdx_plane1",&vars.m_reco_track_trunc_dEdx_p1);
    vars.vertex_tree->Branch("reco_track_trunc_dEdx_plane2",&vars.m_reco_track_trunc_dEdx_p2);

    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_best_plane",&vars.m_reco_track_mean_trunc_dEdx_best_plane);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_plane0",&vars.m_reco_track_mean_trunc_dEdx_p0);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_plane1",&vars.m_reco_track_mean_trunc_dEdx_p1);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_plane2",&vars.m_reco_track_mean_trunc_dEdx_p2);

    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half_best_plane",&vars.m_reco_track_mean_trunc_dEdx_end_half_best_plane);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half_plane0",&vars.m_reco_track_mean_trunc_dEdx_end_half_p0);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half_plane1",&vars.m_reco_track_mean_trunc_dEdx_end_half_p1);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_start_half_plane2",&vars.m_reco_track_mean_trunc_dEdx_end_half_p2);

    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half_best_plane",&vars.m_reco_track_mean_trunc_dEdx_start_half_best_plane);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half_plane0",&vars.m_reco_track_mean_trunc_dEdx_start_half_p0);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half_plane1",&vars.m_reco_track_mean_trunc_dEdx_start_half_p1);
    vars.vertex_tree->Branch("reco_track_mean_trunc_dEdx_end_half_plane2",&vars.m_reco_track_mean_trunc_dEdx_start_half_p2);

    vars.vertex_tree->Branch("reco_track_trunc_PIDA_best_plane",&vars.m_reco_track_trunc_PIDA_best_plane);
    vars.vertex_tree->Branch("reco_track_trunc_PIDA_plane0",&vars.m_reco_track_trunc_PIDA_p0);
    vars.vertex_tree->Branch("reco_track_trunc_PIDA_plane1",&vars.m_reco_track_trunc_PIDA_p1);
    vars.vertex_tree->Branch("reco_track_trunc_PIDA_plane2",&vars.m_reco_track_trunc_PIDA_p2);

    vars.vertex_tree->Branch("reco_track_resrange_best_plane",&vars.m_reco_track_resrange_best_plane);
    vars.vertex_tree->Branch("reco_track_resrange_plane0",&vars.m_reco_track_resrange_p0);
    vars.vertex_tree->Branch("reco_track_resrange_plane1",&vars.m_reco_track_resrange_p1);
    vars.vertex_tree->Branch("reco_track_resrange_plane2",&vars.m_reco_track_resrange_p2);

    vars.vertex_tree->Branch("reco_track_dEdx_best_plane",&vars.m_reco_track_dEdx_best_plane);
    vars.vertex_tree->Branch("reco_track_dEdx_plane0",&vars.m_reco_track_dEdx_p0);
    vars.vertex_tree->Branch("reco_track_dEdx_plane1",&vars.m_reco_track_dEdx_p1);
    vars.vertex_tree->Branch("reco_track_dEdx_plane2",&vars.m_reco_track_dEdx_p2);

    vars.vertex_tree->Branch("reco_track_num_calo_hits_plane0",&vars.m_reco_track_num_calo_hits_p0);
    vars.vertex_tree->Branch("reco_track_num_calo_hits_plane1",&vars.m_reco_track_num_calo_hits_p1);
    vars.vertex_tree->Branch("reco_track_num_calo_hits_plane2",&vars.m_reco_track_num_calo_hits_p2);



    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_mu_plane0",&vars.m_reco_track_pid_bragg_likelihood_mu_plane0);
    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_mu_plane1",&vars.m_reco_track_pid_bragg_likelihood_mu_plane1);
    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_mu_plane2",&vars.m_reco_track_pid_bragg_likelihood_mu_plane2);
    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_p_plane0",&vars.m_reco_track_pid_bragg_likelihood_p_plane0);
    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_p_plane1",&vars.m_reco_track_pid_bragg_likelihood_p_plane1);
    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_p_plane2",&vars.m_reco_track_pid_bragg_likelihood_p_plane2);
    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_mip_plane0",&vars.m_reco_track_pid_bragg_likelihood_mip_plane0);
    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_mip_plane1",&vars.m_reco_track_pid_bragg_likelihood_mip_plane1);
    vars.vertex_tree->Branch("reco_track_pid_bragg_likelihood_mip_plane2",&vars.m_reco_track_pid_bragg_likelihood_mip_plane2);
    vars.vertex_tree->Branch("reco_track_pid_chi2_mu_plane0",&vars.m_reco_track_pid_chi2_mu_plane0);
    vars.vertex_tree->Branch("reco_track_pid_chi2_mu_plane1",&vars.m_reco_track_pid_chi2_mu_plane1);
    vars.vertex_tree->Branch("reco_track_pid_chi2_mu_plane2",&vars.m_reco_track_pid_chi2_mu_plane2);
    vars.vertex_tree->Branch("reco_track_pid_chi2_p_plane0",&vars.m_reco_track_pid_chi2_p_plane0);
    vars.vertex_tree->Branch("reco_track_pid_chi2_p_plane1",&vars.m_reco_track_pid_chi2_p_plane1);
    vars.vertex_tree->Branch("reco_track_pid_chi2_p_plane2",&vars.m_reco_track_pid_chi2_p_plane2);
    vars.vertex_tree->Branch("reco_track_pid_pida_plane0",&vars.m_reco_track_pid_pida_plane0);
    vars.vertex_tree->Branch("reco_track_pid_pida_plane1",&vars.m_reco_track_pid_pida_plane1);
    vars.vertex_tree->Branch("reco_track_pid_pida_plane2",&vars.m_reco_track_pid_pida_plane2);
    vars.vertex_tree->Branch("reco_track_pid_three_plane_proton_pid",&vars.m_reco_track_pid_three_plane_proton_pid);

//    vars.vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane0",&vars.m_reco_track_end_to_nearest_dead_wire_plane0);
//    vars.vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane1",&vars.m_reco_track_end_to_nearest_dead_wire_plane1);
//    vars.vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane2",&vars.m_reco_track_end_to_nearest_dead_wire_plane2);

    vars.vertex_tree->Branch("reco_track_sliceId",& vars.m_reco_track_sliceId);
    vars.vertex_tree->Branch("reco_track_nuscore",& vars.m_reco_track_nuscore);
    vars.vertex_tree->Branch("reco_track_isclearcosmic",& vars.m_reco_track_isclearcosmic);
    vars.vertex_tree->Branch("reco_track_trackscore",& vars.m_reco_track_trackscore);
    vars.vertex_tree->Branch("reco_track_pfparticle_pdg",& vars.m_reco_track_pfparticle_pdg);
    vars.vertex_tree->Branch("reco_track_is_nuslice",& vars.m_reco_track_is_nuslice);

    vars.vertex_tree->Branch("sim_track_matched",&vars.m_sim_track_matched);
    vars.vertex_tree->Branch("sim_track_overlay_fraction",&vars.m_sim_track_overlay_fraction);
    vars.vertex_tree->Branch("sim_track_energy",&vars.m_sim_track_energy);
    vars.vertex_tree->Branch("sim_track_kinetic_energy",&vars.m_sim_track_kinetic_energy);
    vars.vertex_tree->Branch("sim_track_mass",&vars.m_sim_track_mass);
    vars.vertex_tree->Branch("sim_track_pdg",&vars.m_sim_track_pdg);
    vars.vertex_tree->Branch("sim_track_parent_pdg",&vars.m_sim_track_parent_pdg);
    vars.vertex_tree->Branch("sim_track_origin",&vars.m_sim_track_origin);
    vars.vertex_tree->Branch("sim_track_process",&vars.m_sim_track_process);
    vars.vertex_tree->Branch("sim_track_startx",&vars.m_sim_track_startx);
    vars.vertex_tree->Branch("sim_track_starty",&vars.m_sim_track_starty);
    vars.vertex_tree->Branch("sim_track_startz",&vars.m_sim_track_startz);
    vars.vertex_tree->Branch("sim_track_px",&vars.m_sim_track_px);
    vars.vertex_tree->Branch("sim_track_py",&vars.m_sim_track_py);
    vars.vertex_tree->Branch("sim_track_pz",&vars.m_sim_track_pz);
    vars.vertex_tree->Branch("sim_track_endx",&vars.m_sim_track_endx);
    vars.vertex_tree->Branch("sim_track_endy",&vars.m_sim_track_endy);
    vars.vertex_tree->Branch("sim_track_endz",&vars.m_sim_track_endz);
    vars.vertex_tree->Branch("sim_track_length",&vars.m_sim_track_length);

    vars.vertex_tree->Branch("sim_track_trackID",&vars.m_sim_track_trackID);

    vars.vertex_tree->Branch("sim_track_sliceId",& vars.m_sim_track_sliceId);
    vars.vertex_tree->Branch("sim_track_nuscore",& vars.m_sim_track_nuscore);
    vars.vertex_tree->Branch("sim_track_isclearcosmic",& vars.m_sim_track_isclearcosmic);
  }

  void ResizeTracks(size_t size, var_all& vars){
    vars.m_reco_track_length.resize(size);
    vars.m_reco_track_dirx.resize(size);
    vars.m_reco_track_num_daughters.resize(size);
    vars.m_reco_track_daughter_trackscore.resize(size);

    vars.m_reco_track_diry.resize(size);
    vars.m_reco_track_dirz.resize(size);
    vars.m_reco_track_endx.resize(size);
    vars.m_reco_track_endy.resize(size);
    vars.m_reco_track_endz.resize(size);
    vars.m_reco_track_end_dist_to_active_TPC.resize(size);
    vars.m_reco_track_start_dist_to_active_TPC.resize(size);
    vars.m_reco_track_end_dist_to_CPA.resize(size);
    vars.m_reco_track_start_dist_to_CPA.resize(size);
    vars.m_reco_track_end_dist_to_SCB.resize(size);
    vars.m_reco_track_start_dist_to_SCB.resize(size);
    vars.m_reco_track_end_in_SCB.resize(size);
    vars.m_reco_track_start_in_SCB.resize(size);

    vars.m_reco_track_calo_energy_plane0.resize(size);
    vars.m_reco_track_calo_energy_plane1.resize(size);
    vars.m_reco_track_calo_energy_plane2.resize(size);
    vars.m_reco_track_calo_energy_max.resize(size);



    vars.m_reco_track_startx.resize(size);
    vars.m_reco_track_starty.resize(size);
    vars.m_reco_track_startz.resize(size);
    vars.m_reco_track_num_trajpoints.resize(size);
    vars.m_reco_track_num_spacepoints.resize(size);
    vars.m_reco_track_proton_kinetic_energy.resize(size);
    vars.m_reco_track_ordered_energy_index.resize(size);
    vars.m_reco_track_ordered_displacement_index.resize(size);


    vars.m_reco_track_spacepoint_principal0.resize(size);
    vars.m_reco_track_spacepoint_principal1.resize(size);
    vars.m_reco_track_spacepoint_principal2.resize(size);

    vars.m_reco_track_spacepoint_chi.resize(size);
    vars.m_reco_track_spacepoint_max_dist.resize(size);

    vars.m_reco_track_theta_yz.resize(size);
    vars.m_reco_track_phi_yx.resize(size);

    vars.m_reco_track_best_calo_plane.resize(size);

    vars.m_reco_track_mean_dEdx_best_plane.resize(size);
    vars.m_reco_track_mean_dEdx_start_half_best_plane.resize(size);
    vars.m_reco_track_mean_dEdx_end_half_best_plane.resize(size);
    vars.m_reco_track_good_calo_best_plane.resize(size);
    vars.m_reco_track_trunc_dEdx_best_plane.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_best_plane.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_start_half_best_plane.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_end_half_best_plane.resize(size);
    vars.m_reco_track_trunc_PIDA_best_plane.resize(size);
    vars.m_reco_track_resrange_best_plane.resize(size);
    vars.m_reco_track_dEdx_best_plane.resize(size);


    vars.m_reco_track_mean_dEdx_p0.resize(size);
    vars.m_reco_track_mean_dEdx_start_half_p0.resize(size);
    vars.m_reco_track_mean_dEdx_end_half_p0.resize(size);
    vars.m_reco_track_good_calo_p0.resize(size);
    vars.m_reco_track_trunc_dEdx_p0.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_p0.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_start_half_p0.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_end_half_p0.resize(size);
    vars.m_reco_track_trunc_PIDA_p0.resize(size);
    vars.m_reco_track_resrange_p0.resize(size);
    vars.m_reco_track_dEdx_p0.resize(size);

    vars.m_reco_track_mean_dEdx_p1.resize(size);
    vars.m_reco_track_mean_dEdx_start_half_p1.resize(size);
    vars.m_reco_track_mean_dEdx_end_half_p1.resize(size);
    vars.m_reco_track_good_calo_p1.resize(size);
    vars.m_reco_track_trunc_dEdx_p1.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_p1.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_start_half_p1.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_end_half_p1.resize(size);
    vars.m_reco_track_trunc_PIDA_p1.resize(size);
    vars.m_reco_track_resrange_p1.resize(size);
    vars.m_reco_track_dEdx_p1.resize(size);

    vars.m_reco_track_mean_dEdx_p2.resize(size);
    vars.m_reco_track_mean_dEdx_start_half_p2.resize(size);
    vars.m_reco_track_mean_dEdx_end_half_p2.resize(size);
    vars.m_reco_track_good_calo_p2.resize(size);
    vars.m_reco_track_trunc_dEdx_p2.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_p2.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_start_half_p2.resize(size);
    vars.m_reco_track_mean_trunc_dEdx_end_half_p2.resize(size);
    vars.m_reco_track_trunc_PIDA_p2.resize(size);
    vars.m_reco_track_resrange_p2.resize(size);
    vars.m_reco_track_dEdx_p2.resize(size);

    vars.m_reco_track_num_calo_hits_p1.resize(size);
    vars.m_reco_track_num_calo_hits_p0.resize(size);
    vars.m_reco_track_num_calo_hits_p2.resize(size);



    vars.m_sim_track_matched.resize(size);
    vars.m_sim_track_energy.resize(size);
    vars.m_sim_track_mass.resize(size);
    vars.m_sim_track_kinetic_energy.resize(size);
    vars.m_sim_track_pdg.resize(size);
    vars.m_sim_track_parent_pdg.resize(size);
    vars.m_sim_track_origin.resize(size);
    vars.m_sim_track_process.resize(size);
    vars.m_sim_track_startx.resize(size);
    vars.m_sim_track_starty.resize(size);
    vars.m_sim_track_startz.resize(size);
    vars.m_sim_track_endx.resize(size);
    vars.m_sim_track_endy.resize(size);
    vars.m_sim_track_endz.resize(size);
    vars.m_sim_track_length.resize(size);

    vars.m_sim_track_px.resize(size);
    vars.m_sim_track_py.resize(size);
    vars.m_sim_track_pz.resize(size);
    vars.m_sim_track_trackID.resize(size);
    vars.m_sim_track_overlay_fraction.resize(size);

    vars.m_reco_track_pid_bragg_likelihood_mu_plane0.resize(size);
    vars.m_reco_track_pid_bragg_likelihood_mu_plane1.resize(size);
    vars.m_reco_track_pid_bragg_likelihood_mu_plane2.resize(size);
    vars.m_reco_track_pid_bragg_likelihood_p_plane0.resize(size);
    vars.m_reco_track_pid_bragg_likelihood_p_plane1.resize(size);
    vars.m_reco_track_pid_bragg_likelihood_p_plane2.resize(size);
    vars.m_reco_track_pid_bragg_likelihood_mip_plane0.resize(size);
    vars.m_reco_track_pid_bragg_likelihood_mip_plane1.resize(size);
    vars.m_reco_track_pid_bragg_likelihood_mip_plane2.resize(size);
    vars.m_reco_track_pid_chi2_mu_plane0.resize(size);
    vars.m_reco_track_pid_chi2_mu_plane1.resize(size);
    vars.m_reco_track_pid_chi2_mu_plane2.resize(size);
    vars.m_reco_track_pid_chi2_p_plane0.resize(size);
    vars.m_reco_track_pid_chi2_p_plane1.resize(size);
    vars.m_reco_track_pid_chi2_p_plane2.resize(size);
    vars.m_reco_track_pid_pida_plane0.resize(size);
    vars.m_reco_track_pid_pida_plane1.resize(size);
    vars.m_reco_track_pid_pida_plane2.resize(size);
    vars.m_reco_track_pid_three_plane_proton_pid.resize(size);

//    vars.m_reco_track_end_to_nearest_dead_wire_plane0.resize(size);
//    vars.m_reco_track_end_to_nearest_dead_wire_plane1.resize(size);
//    vars.m_reco_track_end_to_nearest_dead_wire_plane2.resize(size);

    vars.m_reco_track_sliceId.resize(size);
    vars.m_reco_track_nuscore.resize(size);
    vars.m_reco_track_isclearcosmic.resize(size);
    vars.m_reco_track_trackscore.resize(size);
    vars.m_reco_track_pfparticle_pdg.resize(size);
    vars.m_reco_track_is_nuslice.resize(size);

    vars.m_sim_track_sliceId.resize(size);
    vars.m_sim_track_nuscore.resize(size);
    vars.m_sim_track_isclearcosmic.resize(size);
  }

  //analyze_Shower.h
  void ClearShowers(var_all& vars){
    vars.m_reco_asso_showers=0;
    vars.m_reco_shower_num_daughters.clear();
    vars.m_reco_shower_daughter_trackscore.clear();

    vars.m_reco_shower3d_exists.clear();

    vars.m_reco_shower3d_startx.clear();
    vars.m_reco_shower3d_starty.clear();
    vars.m_reco_shower3d_startz.clear();
    vars.m_reco_shower3d_dirx.clear();
    vars.m_reco_shower3d_diry.clear();
    vars.m_reco_shower3d_dirz.clear();
    vars.m_reco_shower3d_theta_yz.clear();
    vars.m_reco_shower3d_phi_yx.clear();
    vars.m_reco_shower3d_conversion_distance.clear();
    vars.m_reco_shower3d_impact_parameter.clear();
    vars.m_reco_shower3d_implied_dirx.clear();
    vars.m_reco_shower3d_implied_diry.clear();
    vars.m_reco_shower3d_implied_dirz.clear();
    vars.m_reco_shower3d_openingangle.clear();
    vars.m_reco_shower3d_length.clear();

    vars.m_reco_shower3d_energy_plane0.clear();
    vars.m_reco_shower3d_energy_plane1.clear();
    vars.m_reco_shower3d_energy_plane2.clear();
    vars.m_reco_shower3d_dEdx_plane0.clear();
    vars.m_reco_shower3d_dEdx_plane1.clear();
    vars.m_reco_shower3d_dEdx_plane2.clear();


    vars.m_reco_shower_startx.clear();
    vars.m_reco_shower_starty.clear();
    vars.m_reco_shower_start_dist_to_active_TPC.clear();
    vars.m_reco_shower_start_dist_to_CPA.clear();
    vars.m_reco_shower_start_dist_to_SCB.clear();
    vars.m_reco_shower_start_in_SCB.clear();
    vars.m_reco_shower_end_dist_to_active_TPC.clear();
    vars.m_reco_shower_end_dist_to_SCB.clear();

    vars.m_reco_shower_dirx.clear();
    vars.m_reco_shower_diry.clear();
    vars.m_reco_shower_dirz.clear();
    vars.m_reco_shower_theta_yz.clear();
    vars.m_reco_shower_phi_yx.clear();
    vars.m_reco_shower_conversion_distance.clear();
    vars.m_reco_shower_impact_parameter.clear();
    vars.m_reco_shower_implied_dirx.clear();
    vars.m_reco_shower_implied_diry.clear();
    vars.m_reco_shower_implied_dirz.clear();
    vars.m_reco_shower_openingangle.clear();
    vars.m_reco_shower_length.clear();
    vars.m_reco_shower_delaunay_num_triangles_plane0.clear();
    vars.m_reco_shower_delaunay_num_triangles_plane1.clear();
    vars.m_reco_shower_delaunay_num_triangles_plane2.clear();
    vars.m_reco_shower_num_hits_plane0.clear();
    vars.m_reco_shower_num_hits_plane1.clear();
    vars.m_reco_shower_num_hits_plane2.clear();
    vars.m_reco_shower_delaunay_area_plane0.clear();
    vars.m_reco_shower_delaunay_area_plane1.clear();
    vars.m_reco_shower_delaunay_area_plane2.clear();

    vars.m_reco_shower_kalman_exists.clear();
    vars.m_reco_shower_kalman_median_dEdx_plane0.clear();
    vars.m_reco_shower_kalman_median_dEdx_plane1.clear();
    vars.m_reco_shower_kalman_median_dEdx_plane2.clear();
    vars.m_reco_shower_kalman_median_dEdx_allplane.clear();
    vars.m_reco_shower_kalman_mean_dEdx_plane0.clear();
    vars.m_reco_shower_kalman_mean_dEdx_plane1.clear();
    vars.m_reco_shower_kalman_mean_dEdx_plane2.clear();

    vars.m_sim_shower_energy.clear();
    vars.m_sim_shower_matched.clear();
    vars.m_sim_shower_kinetic_energy.clear();
    vars.m_sim_shower_mass.clear();
    vars.m_sim_shower_pdg.clear();
    vars.m_sim_shower_trackID.clear();
    vars.m_sim_shower_parent_pdg.clear();
    vars.m_sim_shower_parent_trackID.clear();
    vars.m_sim_shower_origin.clear();
    vars.m_sim_shower_process.clear();
    vars.m_sim_shower_end_process.clear();
    vars.m_sim_shower_start_x.clear();
    vars.m_sim_shower_start_y.clear();
    vars.m_sim_shower_start_z.clear();
    vars.m_sim_shower_vertex_x.clear();
    vars.m_sim_shower_vertex_y.clear();
    vars.m_sim_shower_vertex_z.clear();
    vars.m_sim_shower_is_true_shower.clear();
    vars.m_sim_shower_best_matched_plane.clear();
    vars.m_sim_shower_matched_energy_fraction_plane0.clear();
    vars.m_sim_shower_matched_energy_fraction_plane1.clear();
    vars.m_sim_shower_matched_energy_fraction_plane2.clear();
    vars.m_sim_shower_overlay_fraction.clear();
    vars.m_sim_shower_px.clear();
    vars.m_sim_shower_py.clear();
    vars.m_sim_shower_pz.clear();
    vars.m_sim_shower_sliceId.clear();
    vars.m_sim_shower_nuscore.clear();
    vars.m_sim_shower_isclearcosmic.clear();
    vars.m_sim_shower_is_nuslice.clear();



    vars.m_reco_shower_ordered_energy_index.clear();
    vars.m_reco_shower_energy_max.clear();
    vars.m_reco_shower_energy_plane0.clear();
    vars.m_reco_shower_energy_plane1.clear();
    vars.m_reco_shower_energy_plane2.clear();

    vars.m_reco_shower_reclustered_energy_plane0.clear();
    vars.m_reco_shower_reclustered_energy_plane1.clear();
    vars.m_reco_shower_reclustered_energy_plane2.clear();
    vars.m_reco_shower_reclustered_energy_max.clear();

    vars.m_reco_shower_plane0_nhits.clear();
    vars.m_reco_shower_plane1_nhits.clear();
    vars.m_reco_shower_plane2_nhits.clear();
    vars.m_reco_shower_plane0_meanRMS.clear();
    vars.m_reco_shower_plane1_meanRMS.clear();
    vars.m_reco_shower_plane2_meanRMS.clear();

    vars.m_reco_shower_hit_tick.clear();
    vars.m_reco_shower_hit_wire.clear();
    vars.m_reco_shower_hit_plane.clear();
    vars.m_reco_shower_spacepoint_x.clear();
    vars.m_reco_shower_spacepoint_y.clear();
    vars.m_reco_shower_spacepoint_z.clear();


    vars.m_reco_shower_dQdx_plane0.clear();
    vars.m_reco_shower_dQdx_plane2.clear();
    vars.m_reco_shower_dQdx_plane2.clear();
    vars.m_reco_shower_dEdx_plane0.clear();
    vars.m_reco_shower_dEdx_plane1.clear();
    vars.m_reco_shower_dEdx_plane2.clear();
    vars.m_reco_shower_dEdx_plane0_median.clear();
    vars.m_reco_shower_dEdx_plane1_median.clear();
    vars.m_reco_shower_dEdx_plane2_median.clear();

    vars.m_reco_shower_angle_wrt_wires_plane0.clear();
    vars.m_reco_shower_angle_wrt_wires_plane1.clear();
    vars.m_reco_shower_angle_wrt_wires_plane2.clear();

    vars.m_reco_shower_dEdx_amalgamated.clear();
    vars.m_reco_shower_dEdx_amalgamated_nhits.clear();


    vars.m_reco_shower_dQdx_plane0_median.clear();
    vars.m_reco_shower_dQdx_plane1_median.clear();
    vars.m_reco_shower_dQdx_plane2_median.clear();

    vars.m_reco_shower_dEdx_plane0_mean.clear();
    vars.m_reco_shower_dEdx_plane1_mean.clear();
    vars.m_reco_shower_dEdx_plane2_mean.clear();  
    vars.m_reco_shower_dEdx_plane0_max.clear();
    vars.m_reco_shower_dEdx_plane1_max.clear();
    vars.m_reco_shower_dEdx_plane2_max.clear();  
    vars.m_reco_shower_dEdx_plane0_min.clear();
    vars.m_reco_shower_dEdx_plane1_min.clear();
    vars.m_reco_shower_dEdx_plane2_min.clear();  

    vars.m_reco_shower_dEdx_plane0_nhits.clear();
    vars.m_reco_shower_dEdx_plane1_nhits.clear();
    vars.m_reco_shower_dEdx_plane2_nhits.clear();  

//    vars.m_reco_shower_start_to_nearest_dead_wire_plane0.clear();
//    vars.m_reco_shower_start_to_nearest_dead_wire_plane1.clear();
//    vars.m_reco_shower_start_to_nearest_dead_wire_plane2.clear();

    vars.m_reco_shower_flash_shortest_distz.clear();
    vars.m_reco_shower_flash_shortest_index_z.clear();
    vars.m_reco_shower_flash_shortest_disty.clear();
    vars.m_reco_shower_flash_shortest_index_y.clear();

    vars.m_reco_shower_flash_shortest_distyz.clear();
    vars.m_reco_shower_flash_shortest_index_yz.clear();

    vars.m_reco_shower_sliceId.clear();
    vars.m_reco_shower_nuscore.clear();
    vars.m_reco_shower_isclearcosmic.clear();
    vars.m_reco_shower_is_nuslice.clear();
    vars.m_reco_shower_trackscore.clear();
    vars.m_reco_shower_pfparticle_pdg.clear();

  }

  void CreateShowerBranches(var_all& vars){
    vars.vertex_tree->Branch("reco_asso_showers",&vars.m_reco_asso_showers,"reco_asso_showers/I");
    vars.vertex_tree->Branch("reco_shower_num_daughters",&vars.m_reco_shower_num_daughters);
    vars.vertex_tree->Branch("reco_shower_daughter_trackscore",&vars.m_reco_shower_daughter_trackscore);

    vars.vertex_tree->Branch("reco_shower_length", &vars.m_reco_shower_length);
    vars.vertex_tree->Branch("reco_shower_opening_angle", &vars.m_reco_shower_openingangle);
    vars.vertex_tree->Branch("reco_shower_dirx", &vars.m_reco_shower_dirx);
    vars.vertex_tree->Branch("reco_shower_diry", &vars.m_reco_shower_diry);
    vars.vertex_tree->Branch("reco_shower_dirz", &vars.m_reco_shower_dirz);
    vars.vertex_tree->Branch("reco_shower_startx", &vars.m_reco_shower_startx);
    vars.vertex_tree->Branch("reco_shower_starty", &vars.m_reco_shower_starty);
    vars.vertex_tree->Branch("reco_shower_startz", &vars.m_reco_shower_startz);
    vars.vertex_tree->Branch("reco_shower_start_dist_to_active_TPC", &vars.m_reco_shower_start_dist_to_active_TPC);
    vars.vertex_tree->Branch("reco_shower_start_dist_to_CPA", &vars.m_reco_shower_start_dist_to_CPA);
    vars.vertex_tree->Branch("reco_shower_start_dist_to_SCB",  &vars.m_reco_shower_start_dist_to_SCB);
    vars.vertex_tree->Branch("reco_shower_start_in_SCB",   &vars.m_reco_shower_start_in_SCB);
    vars.vertex_tree->Branch("reco_shower_end_dist_to_active_TPC", &vars.m_reco_shower_end_dist_to_active_TPC);
    vars.vertex_tree->Branch("reco_shower_end_dist_to_SCB",  &vars.m_reco_shower_end_dist_to_SCB);


    vars.vertex_tree->Branch("reco_shower_theta_yz",&vars.m_reco_shower_theta_yz);
    vars.vertex_tree->Branch("reco_shower_phi_yx",&vars.m_reco_shower_phi_yx);
    vars.vertex_tree->Branch("reco_shower_conversion_distance",& vars.m_reco_shower_conversion_distance);
    vars.vertex_tree->Branch("reco_shower_impact_parameter",& vars.m_reco_shower_impact_parameter);
    vars.vertex_tree->Branch("reco_shower_implied_dirx", &vars.m_reco_shower_implied_dirx);
    vars.vertex_tree->Branch("reco_shower_implied_diry", &vars.m_reco_shower_implied_diry);
    vars.vertex_tree->Branch("reco_shower_implied_dirz", &vars.m_reco_shower_implied_dirz);

    vars.vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane0",&vars.m_reco_shower_delaunay_num_triangles_plane0);
    vars.vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane1",&vars.m_reco_shower_delaunay_num_triangles_plane1);
    vars.vertex_tree->Branch("reco_shower_delaunay_num_triangles_plane2",&vars.m_reco_shower_delaunay_num_triangles_plane2);
    vars.vertex_tree->Branch("reco_shower_num_hits_plane0",&vars.m_reco_shower_num_hits_plane0);
    vars.vertex_tree->Branch("reco_shower_num_hits_plane1",&vars.m_reco_shower_num_hits_plane1);
    vars.vertex_tree->Branch("reco_shower_num_hits_plane2",&vars.m_reco_shower_num_hits_plane2);
    vars.vertex_tree->Branch("reco_shower_delaunay_area_plane0",&vars.m_reco_shower_delaunay_area_plane0);
    vars.vertex_tree->Branch("reco_shower_delaunay_area_plane1",&vars.m_reco_shower_delaunay_area_plane1);
    vars.vertex_tree->Branch("reco_shower_delaunay_area_plane2",&vars.m_reco_shower_delaunay_area_plane2);
    //the calorimetry info
    vars.vertex_tree->Branch("reco_shower_energy_max",&vars.m_reco_shower_energy_max);
    vars.vertex_tree->Branch("reco_shower_energy_plane0",&vars.m_reco_shower_energy_plane0);
    vars.vertex_tree->Branch("reco_shower_energy_plane1",&vars.m_reco_shower_energy_plane1);
    vars.vertex_tree->Branch("reco_shower_energy_plane2",&vars.m_reco_shower_energy_plane2);
    vars.vertex_tree->Branch("reco_shower_plane0_nhits",&vars.m_reco_shower_plane0_nhits);
    vars.vertex_tree->Branch("reco_shower_plane1_nhits",&vars.m_reco_shower_plane1_nhits);
    vars.vertex_tree->Branch("reco_shower_plane2_nhits",&vars.m_reco_shower_plane2_nhits);
    vars.vertex_tree->Branch("reco_shower_plane0_meanRMS",&vars.m_reco_shower_plane0_meanRMS);
    vars.vertex_tree->Branch("reco_shower_plane1_meanRMS",&vars.m_reco_shower_plane1_meanRMS);
    vars.vertex_tree->Branch("reco_shower_plane2_meanRMS",&vars.m_reco_shower_plane2_meanRMS);

    vars.vertex_tree->Branch("reco_shower_reclustered_energy_plane0",&vars.m_reco_shower_reclustered_energy_plane0);
    vars.vertex_tree->Branch("reco_shower_reclustered_energy_plane1",&vars.m_reco_shower_reclustered_energy_plane1);
    vars.vertex_tree->Branch("reco_shower_reclustered_energy_plane2",&vars.m_reco_shower_reclustered_energy_plane2);
    vars.vertex_tree->Branch("reco_shower_reclustered_energy_max",&vars.m_reco_shower_reclustered_energy_max);

    vars.vertex_tree->Branch("reco_shower_hit_tick",&vars.m_reco_shower_hit_tick);
    vars.vertex_tree->Branch("reco_shower_hit_wire",&vars.m_reco_shower_hit_wire);
    vars.vertex_tree->Branch("reco_shower_hit_plane",&vars.m_reco_shower_hit_plane);

    vars.vertex_tree->Branch("reco_shower_spacepoint_x",&vars.m_reco_shower_spacepoint_x);
    vars.vertex_tree->Branch("reco_shower_spacepoint_y",&vars.m_reco_shower_spacepoint_y);
    vars.vertex_tree->Branch("reco_shower_spacepoint_z",&vars.m_reco_shower_spacepoint_z);

    vars.vertex_tree->Branch("reco_shower_ordered_energy_index",&vars.m_reco_shower_ordered_energy_index);
    vars.vertex_tree->Branch("i_shr",&vars.m_reco_shower_ordered_energy_index);
    vars.vertex_tree->Branch("reco_shower_dQdx_plane0",&vars.m_reco_shower_dQdx_plane0);
    vars.vertex_tree->Branch("reco_shower_dQdx_plane1",&vars.m_reco_shower_dQdx_plane1);
    vars.vertex_tree->Branch("reco_shower_dQdx_plane2",&vars.m_reco_shower_dQdx_plane2);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane0",&vars.m_reco_shower_dEdx_plane0);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane1",&vars.m_reco_shower_dEdx_plane1);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane2",&vars.m_reco_shower_dEdx_plane2);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane0_median",&vars.m_reco_shower_dEdx_plane0_median);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane1_median",&vars.m_reco_shower_dEdx_plane1_median);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane2_median",&vars.m_reco_shower_dEdx_plane2_median);

    vars.vertex_tree->Branch("reco_shower_angle_wrt_wires_plane0",& vars.m_reco_shower_angle_wrt_wires_plane0);
    vars.vertex_tree->Branch("reco_shower_angle_wrt_wires_plane1",& vars.m_reco_shower_angle_wrt_wires_plane1);
    vars.vertex_tree->Branch("reco_shower_angle_wrt_wires_plane2",& vars.m_reco_shower_angle_wrt_wires_plane2);

    vars.vertex_tree->Branch("reco_shower_dEdx_amalgamated",&vars.m_reco_shower_dEdx_amalgamated);
    vars.vertex_tree->Branch("reco_shower_dEdx_amalgamated_nhits",&vars.m_reco_shower_dEdx_amalgamated_nhits);


    vars.vertex_tree->Branch("reco_shower_dQdx_plane0_median",&vars.m_reco_shower_dQdx_plane0_median);
    vars.vertex_tree->Branch("reco_shower_dQdx_plane1_median",&vars.m_reco_shower_dQdx_plane1_median);
    vars.vertex_tree->Branch("reco_shower_dQdx_plane2_median",&vars.m_reco_shower_dQdx_plane2_median);

    vars.vertex_tree->Branch("reco_shower_dEdx_plane0_mean",&vars.m_reco_shower_dEdx_plane0_mean);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane1_mean",&vars.m_reco_shower_dEdx_plane1_mean);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane2_mean",&vars.m_reco_shower_dEdx_plane2_mean);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane0_max",&vars.m_reco_shower_dEdx_plane0_max);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane1_max",&vars.m_reco_shower_dEdx_plane1_max);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane2_max",&vars.m_reco_shower_dEdx_plane2_max);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane0_min",&vars.m_reco_shower_dEdx_plane0_min);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane1_min",&vars.m_reco_shower_dEdx_plane1_min);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane2_min",&vars.m_reco_shower_dEdx_plane2_min);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane0_nhits",&vars.m_reco_shower_dEdx_plane0_nhits);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane1_nhits",&vars.m_reco_shower_dEdx_plane1_nhits);
    vars.vertex_tree->Branch("reco_shower_dEdx_plane2_nhits",&vars.m_reco_shower_dEdx_plane2_nhits);

//    vars.vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane0",&vars.m_reco_shower_start_to_nearest_dead_wire_plane0);
//    vars.vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane1",&vars.m_reco_shower_start_to_nearest_dead_wire_plane1);
//    vars.vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane2",&vars.m_reco_shower_start_to_nearest_dead_wire_plane2);

    vars.vertex_tree->Branch("reco_shower_flash_shortest_distz",&vars.m_reco_shower_flash_shortest_distz);
    vars.vertex_tree->Branch("reco_shower_flash_shortest_disty",&vars.m_reco_shower_flash_shortest_disty);
    vars.vertex_tree->Branch("reco_shower_flash_shortest_distyz",&vars.m_reco_shower_flash_shortest_distyz);
    vars.vertex_tree->Branch("reco_shower_flash_shortest_index_z",&vars.m_reco_shower_flash_shortest_index_z);
    vars.vertex_tree->Branch("reco_shower_flash_shortest_index_y",&vars.m_reco_shower_flash_shortest_index_y);
    vars.vertex_tree->Branch("reco_shower_flash_shortest_index_yz",&vars.m_reco_shower_flash_shortest_index_yz);

    vars.vertex_tree->Branch("reco_shower_sliceId",& vars.m_reco_shower_sliceId);
    vars.vertex_tree->Branch("reco_shower_nuscore",& vars.m_reco_shower_nuscore);
    vars.vertex_tree->Branch("reco_shower_isclearcosmic",& vars.m_reco_shower_isclearcosmic);
    vars.vertex_tree->Branch("reco_shower_is_nuslice", & vars.m_reco_shower_is_nuslice);
    vars.vertex_tree->Branch("reco_shower_trackscore", & vars.m_reco_shower_trackscore);
    vars.vertex_tree->Branch("reco_shower_pfparticle_pdg", & vars.m_reco_shower_pfparticle_pdg);


    vars.vertex_tree->Branch("reco_shower3d_exists", &vars.m_reco_shower3d_exists);
    vars.vertex_tree->Branch("reco_shower3d_length", &vars.m_reco_shower3d_length);
    vars.vertex_tree->Branch("reco_shower3d_opening_angle", &vars.m_reco_shower3d_openingangle);
    vars.vertex_tree->Branch("reco_shower3d_dirx", &vars.m_reco_shower3d_dirx);
    vars.vertex_tree->Branch("reco_shower3d_diry", &vars.m_reco_shower3d_diry);
    vars.vertex_tree->Branch("reco_shower3d_dirz", &vars.m_reco_shower3d_dirz);
    vars.vertex_tree->Branch("reco_shower3d_startx", &vars.m_reco_shower3d_startx);
    vars.vertex_tree->Branch("reco_shower3d_starty", &vars.m_reco_shower3d_starty);
    vars.vertex_tree->Branch("reco_shower3d_startz", &vars.m_reco_shower3d_startz);
    vars.vertex_tree->Branch("reco_shower3d_theta_yz",&vars.m_reco_shower3d_theta_yz);
    vars.vertex_tree->Branch("reco_shower3d_phi_yx",&vars.m_reco_shower3d_phi_yx);
    vars.vertex_tree->Branch("reco_shower3d_conversion_distance",& vars.m_reco_shower3d_conversion_distance);
    vars.vertex_tree->Branch("reco_shower3d_impact_parameter",& vars.m_reco_shower3d_impact_parameter);
    vars.vertex_tree->Branch("reco_shower3d_implied_dirx", &vars.m_reco_shower3d_implied_dirx);
    vars.vertex_tree->Branch("reco_shower3d_implied_diry", &vars.m_reco_shower3d_implied_diry);
    vars.vertex_tree->Branch("reco_shower3d_implied_dirz", &vars.m_reco_shower3d_implied_dirz);

    vars.vertex_tree->Branch("reco_shower3d_energy_plane0", &vars.m_reco_shower3d_energy_plane0);
    vars.vertex_tree->Branch("reco_shower3d_energy_plane1", &vars.m_reco_shower3d_energy_plane1);
    vars.vertex_tree->Branch("reco_shower3d_energy_plane2", &vars.m_reco_shower3d_energy_plane2);
    vars.vertex_tree->Branch("reco_shower3d_dEdx_plane0", &vars.m_reco_shower3d_dEdx_plane0);
    vars.vertex_tree->Branch("reco_shower3d_dEdx_plane1", &vars.m_reco_shower3d_dEdx_plane1);
    vars.vertex_tree->Branch("reco_shower3d_dEdx_plane2", &vars.m_reco_shower3d_dEdx_plane2);

    vars.vertex_tree->Branch("reco_shower_kalman_exists",&vars.m_reco_shower_kalman_exists);
    vars.vertex_tree->Branch("reco_shower_kalman_dEdx_plane0_median",&vars.m_reco_shower_kalman_median_dEdx_plane0);
    vars.vertex_tree->Branch("reco_shower_kalman_dEdx_plane1_median",&vars.m_reco_shower_kalman_median_dEdx_plane1);
    vars.vertex_tree->Branch("reco_shower_kalman_dEdx_plane2_median",&vars.m_reco_shower_kalman_median_dEdx_plane2);
    vars.vertex_tree->Branch("reco_shower_kalman_dEdx_allplane_median",&vars.m_reco_shower_kalman_median_dEdx_allplane);

    vars.vertex_tree->Branch("reco_shower_kalman_dEdx_plane0_mean",&vars.m_reco_shower_kalman_mean_dEdx_plane0);
    vars.vertex_tree->Branch("reco_shower_kalman_dEdx_plane1_mean",&vars.m_reco_shower_kalman_mean_dEdx_plane1);
    vars.vertex_tree->Branch("reco_shower_kalman_dEdx_plane2_mean",&vars.m_reco_shower_kalman_mean_dEdx_plane2);


    vars.vertex_tree->Branch("sim_shower_matched",&vars.m_sim_shower_matched);
    vars.vertex_tree->Branch("sim_shower_energy",&vars.m_sim_shower_energy);
    vars.vertex_tree->Branch("sim_shower_kinetic_energy",&vars.m_sim_shower_kinetic_energy);
    vars.vertex_tree->Branch("sim_shower_mass",&vars.m_sim_shower_mass);
    vars.vertex_tree->Branch("sim_shower_pdg",&vars.m_sim_shower_pdg);
    vars.vertex_tree->Branch("sim_shower_trackID",&vars.m_sim_shower_trackID);
    vars.vertex_tree->Branch("sim_shower_parent_pdg",&vars.m_sim_shower_parent_pdg);
    vars.vertex_tree->Branch("sim_shower_parent_trackID",&vars.m_sim_shower_parent_trackID);
    vars.vertex_tree->Branch("sim_shower_origin",&vars.m_sim_shower_origin);
    vars.vertex_tree->Branch("sim_shower_process",&vars.m_sim_shower_process);
    vars.vertex_tree->Branch("sim_shower_end_process",&vars.m_sim_shower_end_process);
    vars.vertex_tree->Branch("sim_shower_start_x",&vars.m_sim_shower_start_x);
    vars.vertex_tree->Branch("sim_shower_start_y",&vars.m_sim_shower_start_y);
    vars.vertex_tree->Branch("sim_shower_start_z",&vars.m_sim_shower_start_z);
    vars.vertex_tree->Branch("sim_shower_vertex_x",&vars.m_sim_shower_vertex_x);
    vars.vertex_tree->Branch("sim_shower_vertex_y",&vars.m_sim_shower_vertex_y);
    vars.vertex_tree->Branch("sim_shower_vertex_z",&vars.m_sim_shower_vertex_z);
    vars.vertex_tree->Branch("sim_shower_px",&vars.m_sim_shower_px);
    vars.vertex_tree->Branch("sim_shower_py",&vars.m_sim_shower_py);
    vars.vertex_tree->Branch("sim_shower_pz",&vars.m_sim_shower_pz);

    vars.vertex_tree->Branch("sim_shower_is_true_shower",&vars.m_sim_shower_is_true_shower);
    vars.vertex_tree->Branch("sim_shower_best_matched_plane",&vars.m_sim_shower_best_matched_plane);
    vars.vertex_tree->Branch("sim_shower_matched_energy_fraction_plane0",&vars.m_sim_shower_matched_energy_fraction_plane0);
    vars.vertex_tree->Branch("sim_shower_matched_energy_fraction_plane1",&vars.m_sim_shower_matched_energy_fraction_plane1);
    vars.vertex_tree->Branch("sim_shower_matched_energy_fraction_plane2",&vars.m_sim_shower_matched_energy_fraction_plane2);
    vars.vertex_tree->Branch("sim_shower_overlay_fraction",&vars.m_sim_shower_overlay_fraction);
    vars.vertex_tree->Branch("sim_shower_sliceId", & vars.m_sim_shower_sliceId);
    vars.vertex_tree->Branch("sim_shower_nuscore", & vars.m_sim_shower_nuscore);
    vars.vertex_tree->Branch("sim_shower_isclearcosmic", & vars.m_sim_shower_isclearcosmic);
    vars.vertex_tree->Branch("sim_shower_is_nusclice", & vars.m_sim_shower_is_nuslice);
  }

  void ResizeShowers(size_t size, var_all& vars){
    vars.m_reco_shower_num_daughters.resize(size);
    vars.m_reco_shower_daughter_trackscore.resize(size);

    vars.m_reco_shower_kalman_exists.resize(size);
    vars.m_reco_shower_kalman_median_dEdx_plane0.resize(size);
    vars.m_reco_shower_kalman_median_dEdx_plane1.resize(size);
    vars.m_reco_shower_kalman_median_dEdx_plane2.resize(size);
    vars.m_reco_shower_kalman_median_dEdx_allplane.resize(size);
    vars.m_reco_shower_kalman_mean_dEdx_plane0.resize(size);
    vars.m_reco_shower_kalman_mean_dEdx_plane1.resize(size);
    vars.m_reco_shower_kalman_mean_dEdx_plane2.resize(size);

    vars.m_reco_shower_reclustered_energy_plane0.resize(size);
    vars.m_reco_shower_reclustered_energy_plane1.resize(size);
    vars.m_reco_shower_reclustered_energy_plane2.resize(size);
    vars.m_reco_shower_reclustered_energy_max.resize(size);


    vars.m_reco_shower3d_exists.resize(size);
    vars.m_reco_shower3d_startx.resize(size);
    vars.m_reco_shower3d_starty.resize(size);
    vars.m_reco_shower3d_startz.resize(size);
    vars.m_reco_shower3d_dirx.resize(size);
    vars.m_reco_shower3d_diry.resize(size);
    vars.m_reco_shower3d_dirz.resize(size);
    vars.m_reco_shower3d_theta_yz.resize(size);
    vars.m_reco_shower3d_phi_yx.resize(size);
    vars.m_reco_shower3d_conversion_distance.resize(size);
    vars.m_reco_shower3d_openingangle.resize(size);
    vars.m_reco_shower3d_length.resize(size);
    vars.m_reco_shower3d_impact_parameter.resize(size);
    vars.m_reco_shower3d_implied_dirx.resize(size);
    vars.m_reco_shower3d_implied_diry.resize(size);
    vars.m_reco_shower3d_implied_dirz.resize(size);
    vars.m_reco_shower3d_energy_plane0.resize(size);
    vars.m_reco_shower3d_energy_plane1.resize(size);
    vars.m_reco_shower3d_energy_plane2.resize(size);
    vars.m_reco_shower3d_dEdx_plane0.resize(size);
    vars.m_reco_shower3d_dEdx_plane1.resize(size);
    vars.m_reco_shower3d_dEdx_plane2.resize(size);

    vars.m_reco_shower_start_dist_to_active_TPC.resize(size);
    vars.m_reco_shower_start_dist_to_CPA.resize(size);
    vars.m_reco_shower_start_dist_to_SCB.resize(size);
    vars.m_reco_shower_start_in_SCB.resize(size);

    vars.m_reco_shower_end_dist_to_active_TPC.resize(size);
    vars.m_reco_shower_end_dist_to_SCB.resize(size);


    vars.m_reco_shower_startx.resize(size);
    vars.m_reco_shower_starty.resize(size);
    vars.m_reco_shower_startz.resize(size);
    vars.m_reco_shower_dirx.resize(size);
    vars.m_reco_shower_diry.resize(size);
    vars.m_reco_shower_dirz.resize(size);
    vars.m_reco_shower_theta_yz.resize(size);
    vars.m_reco_shower_phi_yx.resize(size);
    vars.m_reco_shower_conversion_distance.resize(size);
    vars.m_reco_shower_openingangle.resize(size);
    vars.m_reco_shower_length.resize(size);
    vars.m_reco_shower_impact_parameter.resize(size);
    vars.m_reco_shower_implied_dirx.resize(size);
    vars.m_reco_shower_implied_diry.resize(size);
    vars.m_reco_shower_implied_dirz.resize(size);
    vars.m_reco_shower_delaunay_num_triangles_plane0.resize(size);
    vars.m_reco_shower_delaunay_num_triangles_plane1.resize(size);
    vars.m_reco_shower_delaunay_num_triangles_plane2.resize(size);
    vars.m_reco_shower_num_hits_plane0.resize(size);
    vars.m_reco_shower_num_hits_plane1.resize(size);
    vars.m_reco_shower_num_hits_plane2.resize(size);
    vars.m_reco_shower_delaunay_area_plane0.resize(size);
    vars.m_reco_shower_delaunay_area_plane1.resize(size);
    vars.m_reco_shower_delaunay_area_plane2.resize(size);

    vars.m_reco_shower_energy_max.resize(size);
    vars.m_reco_shower_energy_plane0.resize(size);
    vars.m_reco_shower_energy_plane1.resize(size);
    vars.m_reco_shower_energy_plane2.resize(size);

    vars.m_reco_shower_plane0_nhits.resize(size);
    vars.m_reco_shower_plane1_nhits.resize(size);
    vars.m_reco_shower_plane2_nhits.resize(size);

    vars.m_reco_shower_plane0_meanRMS.resize(size);
    vars.m_reco_shower_plane1_meanRMS.resize(size);
    vars.m_reco_shower_plane2_meanRMS.resize(size);



    vars.m_reco_shower_ordered_energy_index.resize(size);
    vars.m_reco_shower_dQdx_plane0.resize(size);
    vars.m_reco_shower_dQdx_plane1.resize(size);
    vars.m_reco_shower_dQdx_plane2.resize(size);
    vars.m_reco_shower_dEdx_plane0.resize(size);
    vars.m_reco_shower_dEdx_plane1.resize(size);
    vars.m_reco_shower_dEdx_plane2.resize(size);
    vars.m_reco_shower_dEdx_plane0_median.resize(size);
    vars.m_reco_shower_dEdx_plane1_median.resize(size);
    vars.m_reco_shower_dEdx_plane2_median.resize(size);

    vars.m_reco_shower_angle_wrt_wires_plane0.resize(size);
    vars.m_reco_shower_angle_wrt_wires_plane1.resize(size);
    vars.m_reco_shower_angle_wrt_wires_plane2.resize(size);

    vars.m_reco_shower_dEdx_amalgamated.resize(size);
    vars.m_reco_shower_dEdx_amalgamated_nhits.resize(size);

    vars.m_reco_shower_dQdx_plane0_median.resize(size);
    vars.m_reco_shower_dQdx_plane1_median.resize(size);
    vars.m_reco_shower_dQdx_plane2_median.resize(size);

    vars.m_reco_shower_dEdx_plane0_min.resize(size);
    vars.m_reco_shower_dEdx_plane1_min.resize(size);
    vars.m_reco_shower_dEdx_plane2_min.resize(size);
    vars.m_reco_shower_dEdx_plane0_max.resize(size);
    vars.m_reco_shower_dEdx_plane1_max.resize(size);
    vars.m_reco_shower_dEdx_plane2_max.resize(size);
    vars.m_reco_shower_dEdx_plane0_mean.resize(size);
    vars.m_reco_shower_dEdx_plane1_mean.resize(size);
    vars.m_reco_shower_dEdx_plane2_mean.resize(size);




    vars.m_reco_shower_dEdx_plane0_nhits.resize(size);
    vars.m_reco_shower_dEdx_plane1_nhits.resize(size);
    vars.m_reco_shower_dEdx_plane2_nhits.resize(size);

//    vars.m_reco_shower_start_to_nearest_dead_wire_plane0.resize(size);
//    vars.m_reco_shower_start_to_nearest_dead_wire_plane1.resize(size);
//    vars.m_reco_shower_start_to_nearest_dead_wire_plane2.resize(size);

    vars.m_reco_shower_flash_shortest_distz.resize(size);
    vars.m_reco_shower_flash_shortest_index_z.resize(size);
    vars.m_reco_shower_flash_shortest_disty.resize(size);
    vars.m_reco_shower_flash_shortest_index_y.resize(size);

    vars.m_reco_shower_flash_shortest_distyz.resize(size);
    vars.m_reco_shower_flash_shortest_index_yz.resize(size);

    vars.m_reco_shower_sliceId.resize(size);
    vars.m_reco_shower_nuscore.resize(size);
    vars.m_reco_shower_isclearcosmic.resize(size);
    vars.m_reco_shower_is_nuslice.resize(size);
    vars.m_reco_shower_trackscore.resize(size);
    vars.m_reco_shower_pfparticle_pdg.resize(size);


    vars.m_sim_shower_energy.resize(size);
    vars.m_sim_shower_matched.resize(size);
    vars.m_sim_shower_kinetic_energy.resize(size);
    vars.m_sim_shower_mass.resize(size);
    vars.m_sim_shower_pdg.resize(size);
    vars.m_sim_shower_trackID.resize(size);
    vars.m_sim_shower_parent_pdg.resize(size);
    vars.m_sim_shower_parent_trackID.resize(size);
    vars.m_sim_shower_origin.resize(size);
    vars.m_sim_shower_process.resize(size);
    vars.m_sim_shower_end_process.resize(size);
    vars.m_sim_shower_start_x.resize(size);
    vars.m_sim_shower_start_y.resize(size);
    vars.m_sim_shower_start_z.resize(size);
    vars.m_sim_shower_vertex_x.resize(size);
    vars.m_sim_shower_vertex_y.resize(size);
    vars.m_sim_shower_vertex_z.resize(size);
    vars.m_sim_shower_is_true_shower.resize(size);
    vars.m_sim_shower_best_matched_plane.resize(size);
    vars.m_sim_shower_matched_energy_fraction_plane0.resize(size);
    vars.m_sim_shower_matched_energy_fraction_plane1.resize(size);
    vars.m_sim_shower_matched_energy_fraction_plane2.resize(size);
    vars.m_sim_shower_overlay_fraction.resize(size);
    vars.m_sim_shower_px.resize(size);
    vars.m_sim_shower_py.resize(size);
    vars.m_sim_shower_pz.resize(size);
    vars.m_sim_shower_sliceId.resize(size);
    vars.m_sim_shower_nuscore.resize(size);
    vars.m_sim_shower_isclearcosmic.resize(size);
    vars.m_sim_shower_is_nuslice.resize(size);
  }

  //analyze_MCTruth.h
  void ClearMCTruths(var_all& vars){
    vars.m_mctruth_num = 0;
    vars.m_mctruth_origin = -99;
    vars.m_mctruth_mode = -99;
    vars.m_mctruth_interaction_type = -99;
    vars.m_mctruth_nu_vertex_x = -9999;
    vars.m_mctruth_nu_vertex_y = -9999;
    vars.m_mctruth_nu_vertex_z = -9999;
    vars.m_mctruth_reco_vertex_dist = -9999;
    vars.m_mctruth_ccnc = -99;
    vars.m_mctruth_qsqr = -99;
    vars.m_mctruth_nu_E = -99;
    vars.m_mctruth_nu_pdg = 0;
    vars.m_mctruth_lepton_pdg = 0;
    vars.m_mctruth_num_daughter_particles = -99;
    vars.m_mctruth_daughters_pdg.clear();
    vars.m_mctruth_daughters_E.clear();

    vars.m_mctruth_daughters_status_code.clear();
    vars.m_mctruth_daughters_trackID.clear();
    vars.m_mctruth_daughters_mother_trackID.clear();
    vars.m_mctruth_daughters_px.clear();
    vars.m_mctruth_daughters_py.clear();
    vars.m_mctruth_daughters_pz.clear();
    vars.m_mctruth_daughters_startx.clear();
    vars.m_mctruth_daughters_starty.clear();
    vars.m_mctruth_daughters_startz.clear();
    vars.m_mctruth_daughters_time.clear();
    vars.m_mctruth_daughters_endx.clear();
    vars.m_mctruth_daughters_endy.clear();
    vars.m_mctruth_daughters_endz.clear();
    vars.m_mctruth_daughters_endtime.clear();
    vars.m_mctruth_daughters_process.clear();
    vars.m_mctruth_daughters_end_process.clear();


    vars.m_mctruth_is_delta_radiative = 0;
    vars.m_mctruth_delta_radiative_1g1p_or_1g1n = -999;

    vars.m_mctruth_delta_photon_energy=-999;
    vars.m_mctruth_delta_proton_energy=-999;
    vars.m_mctruth_delta_neutron_energy=-999;

    vars.m_mctruth_num_exiting_photons =0;
    vars.m_mctruth_num_exiting_protons =0;
    vars.m_mctruth_num_exiting_pi0 =0;
    vars.m_mctruth_num_exiting_pipm =0;
    vars.m_mctruth_num_exiting_neutrons=0;
    vars.m_mctruth_num_exiting_delta0=0;
    vars.m_mctruth_num_exiting_deltapm=0;
    vars.m_mctruth_num_exiting_deltapp=0;

    vars.m_mctruth_num_reconstructable_protons = 0;

    vars.m_mctruth_is_reconstructable_1g1p = 0;
    vars.m_mctruth_is_reconstructable_1g0p = 0;

    vars.m_mctruth_leading_exiting_proton_energy = -9999;

    vars.m_mctruth_exiting_pi0_E.clear();
    vars.m_mctruth_exiting_pi0_mom.clear();
    vars.m_mctruth_exiting_pi0_px.clear();
    vars.m_mctruth_exiting_pi0_py.clear();
    vars.m_mctruth_exiting_pi0_pz.clear();

    vars.m_mctruth_pi0_leading_photon_energy = -9999;
    vars.m_mctruth_pi0_subleading_photon_energy = -9999;
    vars.m_mctruth_pi0_leading_photon_end_process = "none";
    vars.m_mctruth_pi0_subleading_photon_end_process = "none";
    vars.m_mctruth_pi0_leading_photon_end = {-9999,-9999,-9999};
    vars.m_mctruth_pi0_leading_photon_start = {-9999,-9999,-9999};
    vars.m_mctruth_pi0_subleading_photon_end = {-9999,-9999,-9999};
    vars.m_mctruth_pi0_subleading_photon_start = {-9999,-9999,-9999};
    vars.m_mctruth_pi0_leading_photon_exiting_TPC = -999;
    vars.m_mctruth_pi0_subleading_photon_exiting_TPC = -999;
    vars.m_mctruth_pi0_leading_photon_mom = {-9999,-9999,-9999};
    vars.m_mctruth_pi0_subleading_photon_mom = {-9999,-9999,-9999};

    vars.m_mctruth_exiting_delta0_num_daughters.clear();

    vars.m_mctruth_exiting_photon_mother_trackID.clear();
    vars.m_mctruth_exiting_photon_trackID.clear();
    vars.m_mctruth_exiting_photon_from_delta_decay.clear();
    vars.m_mctruth_exiting_photon_energy.clear();
    vars.m_mctruth_exiting_photon_px.clear();
    vars.m_mctruth_exiting_photon_py.clear();
    vars.m_mctruth_exiting_photon_pz.clear();

    vars.m_mctruth_exiting_proton_mother_trackID.clear();
    vars.m_mctruth_exiting_proton_trackID.clear();
    vars.m_mctruth_exiting_proton_from_delta_decay.clear();
    vars.m_mctruth_exiting_proton_energy.clear();
    vars.m_mctruth_exiting_proton_px.clear();
    vars.m_mctruth_exiting_proton_py.clear();
    vars.m_mctruth_exiting_proton_pz.clear();

    vars.m_mctruth_exiting_neutron_mother_trackID.clear();
    vars.m_mctruth_exiting_neutron_trackID.clear();
    vars.m_mctruth_exiting_neutron_from_delta_decay.clear();
    vars.m_mctruth_exiting_neutron_energy.clear();
    vars.m_mctruth_exiting_neutron_px.clear();
    vars.m_mctruth_exiting_neutron_py.clear();
    vars.m_mctruth_exiting_neutron_pz.clear();
  }

  void CreateMCTruthBranches(var_all& vars){
    vars.vertex_tree->Branch("mctruth_num",&vars.m_mctruth_num);
    vars.vertex_tree->Branch("mctruth_origin",&vars.m_mctruth_origin);
    vars.vertex_tree->Branch("mctruth_nu_pdg",&vars.m_mctruth_nu_pdg);
    vars.vertex_tree->Branch("mctruth_nu_E",&vars.m_mctruth_nu_E);

    vars.vertex_tree->Branch("mctruth_nu_vertex_x",&vars.m_mctruth_nu_vertex_x);
    vars.vertex_tree->Branch("mctruth_nu_vertex_y",&vars.m_mctruth_nu_vertex_y);
    vars.vertex_tree->Branch("mctruth_nu_vertex_z",&vars.m_mctruth_nu_vertex_z);
    vars.vertex_tree->Branch("mctruth_reco_vertex_dist",&vars.m_mctruth_reco_vertex_dist);

    vars.vertex_tree->Branch("mctruth_lepton_pdg",&vars.m_mctruth_lepton_pdg);
    vars.vertex_tree->Branch("mctruth_lepton_E",&vars.m_mctruth_lepton_E);
    vars.vertex_tree->Branch("mctruth_mode",&vars.m_mctruth_mode);
    vars.vertex_tree->Branch("mctruth_qsqr",&vars.m_mctruth_qsqr);
    vars.vertex_tree->Branch("mctruth_cc_or_nc",&vars.m_mctruth_ccnc);
    vars.vertex_tree->Branch("mctruth_interaction_type",&vars.m_mctruth_interaction_type);

    vars.vertex_tree->Branch("mctruth_num_daughter_particles",&vars.m_mctruth_num_daughter_particles);
    vars.vertex_tree->Branch("mctruth_daughters_pdg",&vars.m_mctruth_daughters_pdg);
    vars.vertex_tree->Branch("mctruth_daughters_E",&vars.m_mctruth_daughters_E);
    vars.vertex_tree->Branch("mctruth_daughters_status_code",&vars.m_mctruth_daughters_status_code);
    vars.vertex_tree->Branch("mctruth_daughters_trackID",&vars.m_mctruth_daughters_trackID);
    vars.vertex_tree->Branch("mctruth_daughters_mother_trackID",&vars.m_mctruth_daughters_mother_trackID);
    vars.vertex_tree->Branch("mctruth_daughters_px",&vars.m_mctruth_daughters_px);
    vars.vertex_tree->Branch("mctruth_daughters_py",&vars.m_mctruth_daughters_py);
    vars.vertex_tree->Branch("mctruth_daughters_pz",&vars.m_mctruth_daughters_pz);
    vars.vertex_tree->Branch("mctruth_daughters_startx",&vars.m_mctruth_daughters_startx);
    vars.vertex_tree->Branch("mctruth_daughters_starty",&vars.m_mctruth_daughters_starty);
    vars.vertex_tree->Branch("mctruth_daughters_startz",&vars.m_mctruth_daughters_startz);
    vars.vertex_tree->Branch("mctruth_daughters_time",&vars.m_mctruth_daughters_time);
    vars.vertex_tree->Branch("mctruth_daughters_endx",&vars.m_mctruth_daughters_endx);
    vars.vertex_tree->Branch("mctruth_daughters_endy",&vars.m_mctruth_daughters_endy);
    vars.vertex_tree->Branch("mctruth_daughters_endz",&vars.m_mctruth_daughters_endz);
    vars.vertex_tree->Branch("mctruth_daughters_endtime",&vars.m_mctruth_daughters_endtime);
    vars.vertex_tree->Branch("mctruth_daughters_process",&vars.m_mctruth_daughters_process);
    vars.vertex_tree->Branch("mctruth_daughters_end_process",&vars.m_mctruth_daughters_end_process);




    vars.vertex_tree->Branch("mctruth_num_exiting_protons",&vars.m_mctruth_num_exiting_protons);
    vars.vertex_tree->Branch("mctruth_num_exiting_photons",&vars.m_mctruth_num_exiting_photons);
    vars.vertex_tree->Branch("mctruth_num_exiting_neutrons",&vars.m_mctruth_num_exiting_neutrons);
    vars.vertex_tree->Branch("mctruth_num_exiting_pi0",&vars.m_mctruth_num_exiting_pi0);
    vars.vertex_tree->Branch("mctruth_num_exiting_pipm",&vars.m_mctruth_num_exiting_pipm);
    vars.vertex_tree->Branch("mctruth_num_exiting_delta0",&vars.m_mctruth_num_exiting_delta0);
    vars.vertex_tree->Branch("mctruth_num_exiting_deltapm",&vars.m_mctruth_num_exiting_deltapm);
    vars.vertex_tree->Branch("mctruth_num_exiting_deltapp",&vars.m_mctruth_num_exiting_deltapp);

    vars.vertex_tree->Branch("mctruth_leading_exiting_proton_energy",&vars.m_mctruth_leading_exiting_proton_energy);
    vars.vertex_tree->Branch("mctruth_is_delta_radiative",&vars.m_mctruth_is_delta_radiative);
    vars.vertex_tree->Branch("mctruth_delta_radiative_1g1p_or_1g1n",&vars.m_mctruth_delta_radiative_1g1p_or_1g1n);
    vars.vertex_tree->Branch("mctruth_delta_photon_energy",&vars.m_mctruth_delta_photon_energy);
    vars.vertex_tree->Branch("mctruth_delta_proton_energy",&vars.m_mctruth_delta_proton_energy);
    vars.vertex_tree->Branch("mctruth_delta_neutron_energy",&vars.m_mctruth_delta_neutron_energy);
    vars.vertex_tree->Branch("mctruth_exiting_delta0_num_daughters",&vars.m_mctruth_exiting_delta0_num_daughters);

    vars.vertex_tree->Branch("mctruth_exiting_photon_trackID",&vars.m_mctruth_exiting_photon_trackID);
    vars.vertex_tree->Branch("mctruth_exiting_photon_mother_trackID",&vars.m_mctruth_exiting_photon_mother_trackID);
    vars.vertex_tree->Branch("mctruth_exiting_photon_from_delta_decay",&vars.m_mctruth_exiting_photon_from_delta_decay);
    vars.vertex_tree->Branch("mctruth_exiting_photon_energy",&vars.m_mctruth_exiting_photon_energy);
    vars.vertex_tree->Branch("mctruth_exiting_photon_px",&vars.m_mctruth_exiting_photon_px);
    vars.vertex_tree->Branch("mctruth_exiting_photon_py",&vars.m_mctruth_exiting_photon_py);
    vars.vertex_tree->Branch("mctruth_exiting_photon_pz",&vars.m_mctruth_exiting_photon_pz);

    vars.vertex_tree->Branch("mctruth_exiting_proton_trackID",&vars.m_mctruth_exiting_proton_trackID);
    vars.vertex_tree->Branch("mctruth_exiting_proton_mother_trackID",&vars.m_mctruth_exiting_proton_mother_trackID);
    vars.vertex_tree->Branch("mctruth_exiting_proton_from_delta_decay",&vars.m_mctruth_exiting_proton_from_delta_decay);
    vars.vertex_tree->Branch("mctruth_exiting_proton_energy",&vars.m_mctruth_exiting_proton_energy);
    vars.vertex_tree->Branch("mctruth_exiting_proton_px",&vars.m_mctruth_exiting_proton_px);
    vars.vertex_tree->Branch("mctruth_exiting_proton_py",&vars.m_mctruth_exiting_proton_py);
    vars.vertex_tree->Branch("mctruth_exiting_proton_pz",&vars.m_mctruth_exiting_proton_pz);

    vars.vertex_tree->Branch("mctruth_exiting_neutron_trackID",&vars.m_mctruth_exiting_neutron_trackID);
    vars.vertex_tree->Branch("mctruth_exiting_neutron_mother_trackID",&vars.m_mctruth_exiting_neutron_mother_trackID);
    vars.vertex_tree->Branch("mctruth_exiting_neutron_from_delta_decay",&vars.m_mctruth_exiting_neutron_from_delta_decay);
    vars.vertex_tree->Branch("mctruth_exiting_neutron_energy",&vars.m_mctruth_exiting_neutron_energy);
    vars.vertex_tree->Branch("mctruth_exiting_neutron_px",&vars.m_mctruth_exiting_neutron_px);
    vars.vertex_tree->Branch("mctruth_exiting_neutron_py",&vars.m_mctruth_exiting_neutron_py);
    vars.vertex_tree->Branch("mctruth_exiting_neutron_pz",&vars.m_mctruth_exiting_neutron_pz);


    vars.vertex_tree->Branch("mctruth_is_reconstructable_1g1p",&vars.m_mctruth_is_reconstructable_1g1p);

    vars.vertex_tree->Branch("mctruth_is_reconstructable_1g1p",&vars.m_mctruth_is_reconstructable_1g1p);
    vars.vertex_tree->Branch("mctruth_num_reconstructable_protons",&vars.m_mctruth_num_reconstructable_protons);

    vars.vertex_tree->Branch("mctruth_pi0_leading_photon_energy",&vars.m_mctruth_pi0_leading_photon_energy);
    vars.vertex_tree->Branch("mctruth_pi0_leading_photon_mom",&vars.m_mctruth_pi0_leading_photon_mom);
    vars.vertex_tree->Branch("mctruth_pi0_leading_photon_start",&vars.m_mctruth_pi0_leading_photon_start);
    vars.vertex_tree->Branch("mctruth_pi0_leading_photon_end",&vars.m_mctruth_pi0_leading_photon_end);
    vars.vertex_tree->Branch("mctruth_pi0_leading_photon_exiting_TPC",&vars.m_mctruth_pi0_leading_photon_exiting_TPC);
    vars.vertex_tree->Branch("mctruth_pi0_subleading_photon_energy",&vars.m_mctruth_pi0_subleading_photon_energy);
    vars.vertex_tree->Branch("mctruth_pi0_subleading_photon_mom",&vars.m_mctruth_pi0_subleading_photon_mom);
    vars.vertex_tree->Branch("mctruth_pi0_subleading_photon_end_process",&vars.m_mctruth_pi0_subleading_photon_end_process);
    vars.vertex_tree->Branch("mctruth_pi0_subleading_photon_start",&vars.m_mctruth_pi0_subleading_photon_start);
    vars.vertex_tree->Branch("mctruth_pi0_subleading_photon_end",&vars.m_mctruth_pi0_subleading_photon_end);
    vars.vertex_tree->Branch("mctruth_pi0_subleading_photon_exiting_TPC",&vars.m_mctruth_pi0_subleading_photon_exiting_TPC);


    vars.vertex_tree->Branch("mctruth_exiting_pi0_E",&vars.m_mctruth_exiting_pi0_E);
    vars.vertex_tree->Branch("mctruth_exiting_pi0_mom",&vars.m_mctruth_exiting_pi0_mom);
    vars.vertex_tree->Branch("mctruth_exiting_pi0_px",&vars.m_mctruth_exiting_pi0_px);
    vars.vertex_tree->Branch("mctruth_exiting_pi0_py",&vars.m_mctruth_exiting_pi0_py);
    vars.vertex_tree->Branch("mctruth_exiting_pi0_pz",&vars.m_mctruth_exiting_pi0_pz);
  }

  void ResizeMCTruths(size_t size, var_all& vars){
    vars.m_mctruth_daughters_pdg.resize(size);
    vars.m_mctruth_daughters_E.resize(size);
    vars.m_mctruth_daughters_status_code.resize(size);
    vars.m_mctruth_daughters_trackID.resize(size);
    vars.m_mctruth_daughters_mother_trackID.resize(size);
    vars.m_mctruth_daughters_px.resize(size);
    vars.m_mctruth_daughters_py.resize(size);
    vars.m_mctruth_daughters_pz.resize(size);
    vars.m_mctruth_daughters_startx.resize(size);
    vars.m_mctruth_daughters_starty.resize(size);
    vars.m_mctruth_daughters_startz.resize(size);
    vars.m_mctruth_daughters_time.resize(size);
    vars.m_mctruth_daughters_endx.resize(size);
    vars.m_mctruth_daughters_endy.resize(size);
    vars.m_mctruth_daughters_endz.resize(size);
    vars.m_mctruth_daughters_endtime.resize(size);
    vars.m_mctruth_daughters_end_process.resize(size);
    vars.m_mctruth_daughters_process.resize(size);
  }

  //analyze_EventWeight.h
  void ClearEventWeightBranches(var_all& vars){
    vars.m_mcflux_nu_pos_x=-9999;
    vars.m_mcflux_nu_pos_y=-9999;
    vars.m_mcflux_nu_pos_z=-9999;
    vars.m_mcflux_nu_mom_x=-9999;
    vars.m_mcflux_nu_mom_y=-9999;
    vars.m_mcflux_nu_mom_z=-9999;
    vars.m_mcflux_nu_mom_z=-9999;
    vars.m_mcflux_nu_mom_E=-9999;
    vars.m_mcflux_ntype=0;
    vars.m_mcflux_ptype=0;
    vars.m_mcflux_nimpwt=-9999;
    vars.m_mcflux_dk2gen=-9999;
    vars.m_mcflux_nenergyn=-9999;
    vars.m_mcflux_tpx=-9999;
    vars.m_mcflux_tpy=-9999;
    vars.m_mcflux_tpz=-9999;
    vars.m_mcflux_vx=-9999;
    vars.m_mcflux_vy=-9999;
    vars.m_mcflux_vz=-9999;
    vars.m_mcflux_tptype=0;
    vars.m_mctruth_nparticles=0;
//    vars.fmcweight.clear();

	//vars.m_mctruth_particles_track_ID[];
	//vars.m_mctruth_particles_pdg_code[];
	//vars.m_mctruth_particles_mother[];
	//vars.m_mctruth_particles_status_code[];
	//vars.m_mctruth_particles_num_daughters[]; //other similar variables
	//vars.m_mctruth_particles_daughters[];
    //vars.m_mctruth_particles_Gvx.clear();
    //vars.m_mctruth_particles_Gvy.clear();
    //vars.m_mctruth_particles_Gvz.clear();
    //vars.m_mctruth_particles_Gvt.clear();
    //vars.m_mctruth_particles_px0.clear();
    //vars.m_mctruth_particles_py0.clear();
    //vars.m_mctruth_particles_pz0.clear();
    //vars.m_mctruth_particles_e0.clear();
    ////int vars.m_mctruth_particles_rescatter.clear();
    //vars.m_mctruth_particles_polx.clear();
    //vars.m_mctruth_particles_poly.clear();
    //vars.m_mctruth_particles_polz.clear();

    //int vars.m_mctruth_neutrino_CCNC;
    //int vars.m_mctruth_neutrino_mode: "vars.m_mctruth_neutrino_mode" //declared in mctruth vars
    //vars.m_mctruth_neutrino_interactionType: "vars.m_mctruth_neutrino_interactionType" 
    //int vars.m_mctruth_neutrino_target.clear();
    //int vars.m_mctruth_neutrino_nucleon.clear();
    //int vars.m_mctruth_neutrino_quark.clear();
    //vars.m_mctruth_neutrino_w.clear();
    //vars.m_mctruth_neutrino_x.clear();
    //vars.m_mctruth_neutrino_y.clear();
    //vars.m_mctruth_neutrino_QSqr: "vars.m_mctruth_neutrino_QSqr"
    vars.m_gtruth_is_sea_quark=false;
    vars.m_gtruth_tgt_pdg=0;
    vars.m_gtruth_tgt_Z = -9999;
    vars.m_gtruth_tgt_A = -9999;
    vars.m_gtruth_tgt_p4_x = -9999;
    vars.m_gtruth_tgt_p4_y = -9999;
    vars.m_gtruth_tgt_p4_z = -9999;
    vars.m_gtruth_tgt_p4_E = -9999;
    vars.m_gtruth_weight=-9999;
    vars.m_gtruth_probability=-9999;
    vars.m_gtruth_xsec=-9999;
    vars.m_gtruth_diff_xsec=-9999;
    vars.m_gtruth_gphase_space=-9999;
    vars.m_gtruth_vertex_x=-9999;
    vars.m_gtruth_vertex_y=-9999;
    vars.m_gtruth_vertex_z=-9999;
    vars.m_gtruth_vertex_T=-9999;
    vars.m_gtruth_gscatter=-9999;
    vars.m_gtruth_gint=-9999;
    vars.m_gtruth_res_num=-9999;
    vars.m_gtruth_num_piplus=-9999;
    vars.m_gtruth_num_pi0=-9999;
    vars.m_gtruth_num_piminus=-9999;
    vars.m_gtruth_num_proton=-9999;
    vars.m_gtruth_num_neutron=-9999;
    vars.m_gtruth_is_charm=false;
    vars.m_gtruth_is_strange=false;
    vars.m_gtruth_charm_hadron_pdg = -9999;
    vars.m_gtruth_strange_hadron_pdg = -9999;
    vars.m_gtruth_decay_mode = -9999;
    vars.m_gtruth_gx=-9999;
    vars.m_gtruth_gy=-9999;
    vars.m_gtruth_gy=-9999;
    vars.m_gtruth_gt=-9999;
    vars.m_gtruth_gw=-9999;
    vars.m_gtruth_gQ2=-9999;
    vars.m_gtruth_gq2=-9999;
    vars.m_gtruth_probe_pdg=0;
    vars.m_gtruth_probe_p4_x=-9999;
    vars.m_gtruth_probe_p4_y=-9999;
    vars.m_gtruth_probe_p4_z=-9999;
    vars.m_gtruth_probe_p4_E=-9999;
    vars.m_gtruth_hit_nuc_p4_x=-9999;
    vars.m_gtruth_hit_nuc_p4_y=-9999;
    vars.m_gtruth_hit_nuc_p4_z=-9999;
    vars.m_gtruth_hit_nuc_p4_E=-9999;
    vars.m_gtruth_hit_nuc_pos=-9999;
    vars.m_gtruth_fs_had_syst_p4_x=-9999;
    vars.m_gtruth_fs_had_syst_p4_y=-9999;
    vars.m_gtruth_fs_had_syst_p4_z=-9999;
    vars.m_gtruth_fs_had_syst_p4_E=-9999;
  }

  void CreateEventWeightBranches(var_all& vars){
    //-----------------run info
    vars.eventweight_tree->Branch("run", &vars.m_run_number_eventweight); 
    vars.eventweight_tree->Branch("subrun", &vars.m_subrun_number_eventweight);
    vars.eventweight_tree->Branch("event",  &vars.m_event_number_eventweight);
    //------------------mcflux
    vars.eventweight_tree->Branch("MCFlux_NuPosX",  &vars.m_mcflux_nu_pos_x );
    vars.eventweight_tree->Branch("MCFlux_NuPosY",  &vars.m_mcflux_nu_pos_y );
    vars.eventweight_tree->Branch("MCFlux_NuPosZ",  &vars.m_mcflux_nu_pos_z );
    vars.eventweight_tree->Branch("MCFlux_NuMomX",  &vars.m_mcflux_nu_mom_x);
    vars.eventweight_tree->Branch("MCFlux_NuMomY",  &vars.m_mcflux_nu_mom_y );
    vars.eventweight_tree->Branch("MCFlux_NuMomZ",  &vars.m_mcflux_nu_mom_z);
    vars.eventweight_tree->Branch("MCFlux_NuMomE",  &vars.m_mcflux_nu_mom_E);
    vars.eventweight_tree->Branch("MCFlux_ntype",  &vars.m_mcflux_ntype );
    vars.eventweight_tree->Branch("MCFlux_ptype",  &vars.m_mcflux_ptype );
    vars.eventweight_tree->Branch("MCFlux_nimpwt",  &vars.m_mcflux_nimpwt );
    vars.eventweight_tree->Branch("MCFlux_dk2gen",  &vars.m_mcflux_dk2gen );
    vars.eventweight_tree->Branch("MCFlux_nenergyn",  &vars.m_mcflux_nenergyn); 
    vars.eventweight_tree->Branch("MCFlux_tpx",  &vars.m_mcflux_tpx );
    vars.eventweight_tree->Branch("MCFlux_tpy",  &vars.m_mcflux_tpy );
    vars.eventweight_tree->Branch("MCFlux_tpz",  &vars.m_mcflux_tpz );
    vars.eventweight_tree->Branch("MCFlux_vx",  &vars.m_mcflux_vx);
    vars.eventweight_tree->Branch("MCFlux_vy",  &vars.m_mcflux_vy );
    vars.eventweight_tree->Branch("MCFlux_vz",  &vars.m_mcflux_vz );
    vars.eventweight_tree->Branch("MCFlux_tptype",  & vars.m_mcflux_tptype );
    //---------------mctruth
    vars.eventweight_tree->Branch("MCTruth_NParticles",  &vars.m_mctruth_nparticles); 
    //vars.eventweight_tree->Branch("MCTruth_particles_TrackId",  &vars.m_mctruth_particles_track_Id, "vars.m_mctruth_particles_track_Id[vars.m_mctruth_nparticles]/I" );
    vars.eventweight_tree->Branch("MCTruth_particles_TrackId",  &vars.m_mctruth_particles_track_Id, "MCTruth_particles_TrackId[MCTruth_NParticles]/I" );
    //vars.eventweight_tree->Branch("MCTruth_particles_TrackId",  &vars.m_mctruth_particles_track_Id);
    vars.eventweight_tree->Branch("MCTruth_particles_PdgCode",  &vars.m_mctruth_particles_pdg_code, "MCTruth_particles_PdgCode[MCTruth_NParticles]/I" );
    vars.eventweight_tree->Branch("MCTruth_particles_Mother",  &vars.m_mctruth_particles_mother, "MCTruth_particles_Mother[MCTruth_NParticles]/I" );
    vars.eventweight_tree->Branch("MCTruth_particles_StatusCode",  &vars.m_mctruth_particles_status_code, "MCTruth_particles_StatusCode[MCTruth_NParticles]/I");
    vars.eventweight_tree->Branch("MCTruth_particles_NumberDaughters",  &vars.m_mctruth_particles_num_daughters ,"MCTruth_particles_NumberDaughters[MCTruth_NParticles]/I" );
    vars.eventweight_tree->Branch("MCTruth_particles_Daughters",  &vars.m_mctruth_particles_daughters,"MCTruth_particles_Daughters[MCTruth_NParticles][100]" );
    vars.eventweight_tree->Branch("MCTruth_particles_Gvx",  &vars.m_mctruth_particles_Gvx,"MCTruth_particles_Gvx[MCTruth_NParticles]/D");
    vars.eventweight_tree->Branch("MCTruth_particles_Gvy",  &vars.m_mctruth_particles_Gvy,"MCTruth_particles_Gvy[MCTruth_NParticles]/D" );
    vars.eventweight_tree->Branch("MCTruth_particles_Gvz",  &vars.m_mctruth_particles_Gvz,"MCTruth_particles_Gvz[MCTruth_NParticles]/D" );
    vars.eventweight_tree->Branch("MCTruth_particles_Gvt",  &vars.m_mctruth_particles_Gvt,"MCTruth_particles_Gvt[MCTruth_NParticles]/D" );
    vars.eventweight_tree->Branch("MCTruth_particles_px0",  &vars.m_mctruth_particles_px0, "MCTruth_particles_px0[MCTruth_NParticles]/D"  );
    vars.eventweight_tree->Branch("MCTruth_particles_py0",  &vars.m_mctruth_particles_py0, "MCTruth_particles_py0[MCTruth_NParticles]/D"  );
    vars.eventweight_tree->Branch("MCTruth_particles_pz0",  &vars.m_mctruth_particles_pz0,  "MCTruth_particles_pz0[MCTruth_NParticles]/D"  );
    vars.eventweight_tree->Branch("MCTruth_particles_e0",  &vars.m_mctruth_particles_e0, "MCTruth_particles_e0[MCTruth_NParticles]/D" );
    vars.eventweight_tree->Branch("MCTruth_particles_Rescatter",  &vars.m_mctruth_particles_rescatter,"MCTruth_particles_Rescatter[MCTruth_NParticles]/I" );
    vars.eventweight_tree->Branch("MCTruth_particles_polx",  &vars.m_mctruth_particles_polx, "MCTruth_particles_polx[MCTruth_NParticles]/D" );
    vars.eventweight_tree->Branch("MCTruth_particles_poly",  &vars.m_mctruth_particles_poly, "MCTruth_particles_poly[MCTruth_NParticles]/D" );
    vars.eventweight_tree->Branch("MCTruth_particles_polz",  &vars.m_mctruth_particles_polz, "MCTruth_particles_polz[MCTruth_NParticles]/D");
    vars.eventweight_tree->Branch("MCTruth_neutrino_CCNC",  &vars.m_mctruth_neutrino_ccnc );
    vars.eventweight_tree->Branch("MCTruth_neutrino_mode",  &vars.m_mctruth_neutrino_mode );
    vars.eventweight_tree->Branch("MCTruth_neutrino_interactionType",  &vars.m_mctruth_neutrino_interaction_type );
    vars.eventweight_tree->Branch("MCTruth_neutrino_target",  &vars.m_mctruth_neutrino_target );
    vars.eventweight_tree->Branch("MCTruth_neutrino_nucleon",  &vars.m_mctruth_neutrino_nucleon );
    vars.eventweight_tree->Branch("MCTruth_neutrino_quark",  &vars.m_mctruth_neutrino_quark );
    vars.eventweight_tree->Branch("MCTruth_neutrino_W",  &vars.m_mctruth_neutrino_w );
    vars.eventweight_tree->Branch("MCTruth_neutrino_X",  &vars.m_mctruth_neutrino_x );
    vars.eventweight_tree->Branch("MCTruth_neutrino_Y",  &vars.m_mctruth_neutrino_y );
    vars.eventweight_tree->Branch("MCTruth_neutrino_QSqr",  &vars.m_mctruth_neutrino_qsqr );

    //---------------------gtruth
    vars.eventweight_tree->Branch("GTruth_IsSeaQuark",  &vars.m_gtruth_is_sea_quark );
    vars.eventweight_tree->Branch("GTruth_tgtPDG",  &vars.m_gtruth_tgt_pdg );
    vars.eventweight_tree->Branch("GTruth_tgtA",  &vars.m_gtruth_tgt_A );
    vars.eventweight_tree->Branch("GTruth_tgtZ",  &vars.m_gtruth_tgt_Z );
    vars.eventweight_tree->Branch("GTruth_TgtP4x",  &vars.m_gtruth_tgt_p4_x );
    vars.eventweight_tree->Branch("GTruth_TgtP4y",  &vars.m_gtruth_tgt_p4_y );
    vars.eventweight_tree->Branch("GTruth_TgtP4z",  &vars.m_gtruth_tgt_p4_z );
    vars.eventweight_tree->Branch("GTruth_TgtP4E",  &vars.m_gtruth_tgt_p4_E );
    vars.eventweight_tree->Branch("GTruth_weight",  &vars.m_gtruth_weight );
    vars.eventweight_tree->Branch("GTruth_probability",  &vars.m_gtruth_probability );
    vars.eventweight_tree->Branch("GTruth_Xsec",  &vars.m_gtruth_xsec );
    vars.eventweight_tree->Branch("GTruth_DiffXsec", &vars.m_gtruth_diff_xsec );
    vars.eventweight_tree->Branch("GTruth_GPhaseSpace", &vars.m_gtruth_gphase_space );
    vars.eventweight_tree->Branch("GTruth_vertexX",  &vars.m_gtruth_vertex_x );
    vars.eventweight_tree->Branch("GTruth_vertexY",  &vars.m_gtruth_vertex_y );
    vars.eventweight_tree->Branch("GTruth_vertexZ",  &vars.m_gtruth_vertex_z );
    vars.eventweight_tree->Branch("GTruth_vertexT",  &vars.m_gtruth_vertex_T );
    vars.eventweight_tree->Branch("GTruth_Gscatter", &vars.m_gtruth_gscatter );
    vars.eventweight_tree->Branch("GTruth_Gint",  &vars.m_gtruth_gint );
    vars.eventweight_tree->Branch("GTruth_ResNum", &vars.m_gtruth_res_num); 
    vars.eventweight_tree->Branch("GTruth_NumPiPlus",  &vars.m_gtruth_num_piplus); 
    vars.eventweight_tree->Branch("GTruth_NumPi0",  &vars.m_gtruth_num_pi0);
    vars.eventweight_tree->Branch("GTruth_NumPiMinus",  &vars.m_gtruth_num_piminus);  
    vars.eventweight_tree->Branch("GTruth_NumProton",  &vars.m_gtruth_num_proton );
    vars.eventweight_tree->Branch("GTruth_NumNeutron", &vars.m_gtruth_num_neutron );
    vars.eventweight_tree->Branch("GTruth_IsCharm",  &vars.m_gtruth_is_charm );
    vars.eventweight_tree->Branch("GTruth_IsStrange",  &vars.m_gtruth_is_strange );
    vars.eventweight_tree->Branch("GTruth_StrangeHadronPDG",  &vars.m_gtruth_strange_hadron_pdg );
    vars.eventweight_tree->Branch("GTruth_CharmHadronPDG",  &vars.m_gtruth_charm_hadron_pdg );
    vars.eventweight_tree->Branch("GTruth_DecayMode",&vars.m_gtruth_decay_mode);
    vars.eventweight_tree->Branch("GTruth_gX",  &vars.m_gtruth_gx );
    vars.eventweight_tree->Branch("GTruth_gY", &vars.m_gtruth_gy );
    vars.eventweight_tree->Branch("GTruth_gT", &vars.m_gtruth_gt );
    vars.eventweight_tree->Branch("GTruth_gW", &vars.m_gtruth_gw );
    vars.eventweight_tree->Branch("GTruth_gQ2", &vars.m_gtruth_gQ2 );
    vars.eventweight_tree->Branch("GTruth_gq2",  &vars.m_gtruth_gq2 );
    vars.eventweight_tree->Branch("GTruth_ProbePDG",  &vars.m_gtruth_probe_pdg );
    vars.eventweight_tree->Branch("GTruth_ProbeP4x",  &vars.m_gtruth_probe_p4_x );
    vars.eventweight_tree->Branch("GTruth_ProbeP4y",  &vars.m_gtruth_probe_p4_y );
    vars.eventweight_tree->Branch("GTruth_ProbeP4z",  &vars.m_gtruth_probe_p4_z );
    vars.eventweight_tree->Branch("GTruth_ProbeP4E",  &vars.m_gtruth_probe_p4_E );
    vars.eventweight_tree->Branch("GTruth_HitNucP4x", &vars.m_gtruth_hit_nuc_p4_x );
    vars.eventweight_tree->Branch("GTruth_HitNucP4y", &vars.m_gtruth_hit_nuc_p4_y );
    vars.eventweight_tree->Branch("GTruth_HitNucP4z", &vars.m_gtruth_hit_nuc_p4_z );
    vars.eventweight_tree->Branch("GTruth_HitNucP4E", &vars.m_gtruth_hit_nuc_p4_E );
    vars.eventweight_tree->Branch("GTruth_HitNucPos", &vars.m_gtruth_hit_nuc_pos );
    vars.eventweight_tree->Branch("GTruth_FShadSystP4x", &vars.m_gtruth_fs_had_syst_p4_x );
    vars.eventweight_tree->Branch("GTruth_FShadSystP4y", &vars.m_gtruth_fs_had_syst_p4_y );
    vars.eventweight_tree->Branch("GTruth_FShadSystP4z",  &vars.m_gtruth_fs_had_syst_p4_z );
    vars.eventweight_tree->Branch("GTruth_FShadSystP4E",  &vars.m_gtruth_fs_had_syst_p4_E );
  }

  //analyze_Geant4.h
  void ClearGeant4Branches(var_all& vars){

    vars.m_geant4_pdg.clear();
    vars.m_geant4_trackid.clear();
    vars.m_geant4_mother.clear();
    vars.m_geant4_statuscode.clear();
    vars.m_geant4_E.clear();
    vars.m_geant4_mass.clear();
    vars.m_geant4_px.clear();
    vars.m_geant4_py.clear();
    vars.m_geant4_pz.clear();
    vars.m_geant4_dx.clear();
    vars.m_geant4_dy.clear();
    vars.m_geant4_dz.clear();

    vars.m_geant4_vx.clear();
    vars.m_geant4_vy.clear();
    vars.m_geant4_vz.clear();
    vars.m_geant4_process.clear();
    vars.m_geant4_end_process.clear();

    vars.m_geant4_costheta.clear();
  }

  void CreateGeant4Branches(var_all& vars){
    vars.geant4_tree->Branch("geant4_pdg",&vars.m_geant4_pdg);
    vars.geant4_tree->Branch("geant4_trackid",&vars.m_geant4_trackid);
    vars.geant4_tree->Branch("geant4_mother",&vars.m_geant4_mother);
    vars.geant4_tree->Branch("geant4_statuscode",&vars.m_geant4_statuscode);
    vars.geant4_tree->Branch("geant4_E",&vars.m_geant4_E);
    vars.geant4_tree->Branch("geant4_mass",&vars.m_geant4_mass);
    vars.geant4_tree->Branch("geant4_px", &vars.m_geant4_px);
    vars.geant4_tree->Branch("geant4_py", &vars.m_geant4_py);
    vars.geant4_tree->Branch("geant4_pz", &vars.m_geant4_pz);

    vars.geant4_tree->Branch("geant4_dx", &vars.m_geant4_dx);
    vars.geant4_tree->Branch("geant4_dy", &vars.m_geant4_dy);
    vars.geant4_tree->Branch("geant4_dz", &vars.m_geant4_dz);

    vars.geant4_tree->Branch("geant4_vx", &vars.m_geant4_vx);
    vars.geant4_tree->Branch("geant4_vy", &vars.m_geant4_vy);
    vars.geant4_tree->Branch("geant4_vz", &vars.m_geant4_vz);
    vars.geant4_tree->Branch("geant4_costheta",&vars.m_geant4_costheta);

    vars.geant4_tree->Branch("geant4_end_process", &vars.m_geant4_end_process);
    vars.geant4_tree->Branch("geant4_process", &vars.m_geant4_process);
  }

  //analyze_Slice.h
  void ClearSlices(var_all& vars){
    vars.m_reco_slice_num = 0;
    vars.m_reco_slice_nuscore.clear();
    vars.m_matched_signal_shower_overlay_fraction.clear();
    //std::vector<double> vars.m_matched_signal_shower_conversion_length;
    vars.m_matched_signal_shower_true_E.clear();
    vars.m_matched_signal_shower_nuscore.clear();
    vars.m_matched_signal_shower_sliceId.clear();
    vars.m_matched_signal_shower_is_clearcosmic.clear();
    vars.m_matched_signal_shower_num = 0;
    vars.m_matched_signal_shower_is_nuslice.clear();
    vars.m_matched_signal_shower_tracks_in_slice.clear();
    vars.m_matched_signal_shower_showers_in_slice.clear();

    vars.m_reco_slice_num_pfps.clear();
    vars.m_reco_slice_num_showers.clear();
    vars.m_reco_slice_num_tracks.clear();


    vars.m_matched_signal_track_true_E.clear();
    vars.m_matched_signal_track_nuscore.clear();
    vars.m_matched_signal_track_sliceId.clear();
    vars.m_matched_signal_track_is_clearcosmic.clear();
    vars.m_matched_signal_track_is_nuslice.clear();
    vars.m_matched_signal_track_tracks_in_slice.clear();
    vars.m_matched_signal_track_showers_in_slice.clear();


    vars.m_matched_signal_track_num = 0;  


    //int vars.m_matched_signal_total_num_slices;

    vars.m_reco_1g1p_is_same_slice = false;
    vars.m_reco_1g1p_is_multiple_slices = false;
    vars.m_reco_1g1p_is_nuslice = false;
    vars.m_reco_1g0p_is_nuslice = false;
    vars.m_reco_1g1p_nuscore = -999;
    vars.m_reco_1g0p_nuscore = -999;
    vars.m_is_matched_1g1p = false;
    vars.m_is_matched_1g0p = false;
    vars.m_no_matched_showers = false;
    vars.m_multiple_matched_showers = false;
    vars.m_multiple_matched_tracks = false;


    /*  vars.m_reco_slice_shower_num_matched_signal = -999;
      vars.m_reco_slice_track_num_matched_signal = -999;
      vars.m_reco_slice_shower_matched_sliceId.clear();
      vars.m_reco_slice_track_matched_sliceId.clear();
      vars.m_reco_slice_shower_matched_energy.clear();
      vars.m_reco_slice_track_matched_energy.clear();
      vars.m_reco_slice_shower_matched_conversion.clear();
      vars.m_reco_slice_shower_matched_overlay_frac.clear();
      */  
  }


  void CreateSliceBranches(var_all& vars){
    vars.vertex_tree->Branch("reco_slice_nuscore",&vars.m_reco_slice_nuscore);
    vars.vertex_tree->Branch("reco_slice_num",&vars.m_reco_slice_num);
    vars.vertex_tree->Branch("reco_slice_shower_num_matched_signal",& vars.m_reco_slice_shower_num_matched_signal);
    vars.vertex_tree->Branch("reco_slice_track_num_matched_signal",& vars.m_reco_slice_track_num_matched_signal);

    vars.ncdelta_slice_tree->Branch("matched_signal_shower_overlay_fraction", &vars.m_matched_signal_shower_overlay_fraction);
    //std::vector<double> vars.m_matched_signal_shower_conversion_length;
    vars.ncdelta_slice_tree->Branch("matched_signal_shower_true_E", &vars.m_matched_signal_shower_true_E);
    vars.ncdelta_slice_tree->Branch("matched_signal_shower_nuscore", &vars.m_matched_signal_shower_nuscore);
    vars.ncdelta_slice_tree->Branch("matched_signal_shower_sliceId", &vars.m_matched_signal_shower_sliceId);
    vars.ncdelta_slice_tree->Branch("matched_signal_shower_is_clearcosmic", &vars.m_matched_signal_shower_is_clearcosmic);
    vars.ncdelta_slice_tree->Branch("matched_signal_shower_num", &vars.m_matched_signal_shower_num);
    vars.ncdelta_slice_tree->Branch("matched_signal_shower_is_nuslice", &vars.m_matched_signal_shower_is_nuslice);
    vars.ncdelta_slice_tree->Branch("matched_signal_shower_tracks_in_slice", &vars.m_matched_signal_shower_tracks_in_slice);
    vars.ncdelta_slice_tree->Branch("matched_signal_shower_showers_in_slice", &vars.m_matched_signal_shower_showers_in_slice);

    vars.ncdelta_slice_tree->Branch("reco_slice_num_pfps", & vars.m_reco_slice_num_pfps);
    vars.ncdelta_slice_tree->Branch("reco_slice_num_showers", & vars.m_reco_slice_num_showers);
    vars.ncdelta_slice_tree->Branch("reco_slice_num_tracks", & vars.m_reco_slice_num_tracks);

    // vars.ncdelta_slice_tree->Branch("matched_signal_track_overlay_fraction", &vars.m_matched_signal_track_overlay_fraction);
    vars.ncdelta_slice_tree->Branch("matched_signal_track_true_E", &vars.m_matched_signal_track_true_E);
    vars.ncdelta_slice_tree->Branch("matched_signal_track_nuscore", &vars.m_matched_signal_track_nuscore);
    vars.ncdelta_slice_tree->Branch("matched_signal_track_sliceId", &vars.m_matched_signal_track_sliceId);
    vars.ncdelta_slice_tree->Branch("matched_signal_track_is_clearcosmic", &vars.m_matched_signal_track_is_clearcosmic);
    vars.ncdelta_slice_tree->Branch("matched_signal_track_num", &vars.m_matched_signal_track_num);
    vars.ncdelta_slice_tree->Branch("matched_signal_track_is_nuslice", &vars.m_matched_signal_track_is_nuslice);
    vars.ncdelta_slice_tree->Branch("matched_signal_track_tracks_in_slice", &vars.m_matched_signal_track_tracks_in_slice);
    vars.ncdelta_slice_tree->Branch("matched_signal_track_showers_in_slice", &vars.m_matched_signal_track_showers_in_slice);

    //int vars.m_matched_signal_total_num_slices;
    vars.ncdelta_slice_tree->Branch("reco_1g1p_is_same_slice",&vars.m_reco_1g1p_is_same_slice);
    vars.ncdelta_slice_tree->Branch("reco_1g1p_is_nuslice",&vars.m_reco_1g1p_is_nuslice);
    vars.ncdelta_slice_tree->Branch("reco_1g1p_is_multiple_slices",&vars.m_reco_1g1p_is_multiple_slices);
    vars.ncdelta_slice_tree->Branch("reco_1g1p_nuscore",&vars.m_reco_1g1p_nuscore);
    vars.ncdelta_slice_tree->Branch("is_matched_1g1p",&vars.m_is_matched_1g1p);

    vars.ncdelta_slice_tree->Branch("reco_1g0p_nuscore",&vars.m_reco_1g0p_nuscore);
    vars.ncdelta_slice_tree->Branch("reco_1g0p_is_nuslice",&vars.m_reco_1g0p_is_nuslice);
    vars.ncdelta_slice_tree->Branch("is_matched_1g0p",&vars.m_is_matched_1g0p);

    vars.ncdelta_slice_tree->Branch("no_matched_showers",& vars.m_no_matched_showers);
    vars.ncdelta_slice_tree->Branch("multiple_matched_showers",& vars.m_multiple_matched_showers);
    vars.ncdelta_slice_tree->Branch("multiple_matched_tracks",& vars.m_multiple_matched_tracks);
  }

  void ResizeSlices(size_t size, var_all& vars){
    vars.m_reco_slice_nuscore.resize(size,-999);
    vars.m_reco_slice_num_pfps.resize(size,0);
    vars.m_reco_slice_num_showers.resize(size,0);
    vars.m_reco_slice_num_tracks.resize(size,0);
  }

  void Save_EventMeta( art::Event &evt, var_all& vars){
    
    //Some event based properties
        vars.m_subrun_counts++;
        vars.m_number_of_events++;
        vars.m_run_number  = evt.run();
        vars.m_subrun_number = evt.subRun();
        vars.m_event_number  = evt.id().event();
  }
  
  //set the vertex for now;
  void Save_PFParticleInfo( std::vector<PandoraPFParticle> PPFPs, var_all& vars, para_all& paras){

    int pfp_size = PPFPs.size();

    for(int index = 0; index < pfp_size; index++){

      PandoraPFParticle* temp_p = &PPFPs[index];
      if(!(vars.pfp_w_bestnuID == temp_p->get_SliceID() && temp_p->get_IsNeutrino()) ) continue;
      vars.m_vertex_pos_x = temp_p->get_Vertex_pos()[0];
      vars.m_vertex_pos_y = temp_p->get_Vertex_pos()[1];
      vars.m_vertex_pos_z = temp_p->get_Vertex_pos()[2];
//      std::cout<<"Best NuScore is found, define the vertice as: ("<<temp_p->get_Vertex_pos()[0]<<","<<temp_p->get_Vertex_pos()[1]<<","<<temp_p->get_Vertex_pos()[2]<<")"<<std::endl;
    
      std::vector<double> tmp = {vars.m_vertex_pos_x, vars.m_vertex_pos_y, vars.m_vertex_pos_z};
      vars.m_reco_vertex_in_SCB = distToSCB(vars.m_reco_vertex_dist_to_SCB,tmp, paras);
      vars.m_reco_vertex_dist_to_active_TPC =  distToTPCActive(tmp, paras);
      vars.m_reco_vertex_dist_to_CPA =  distToCPA(tmp, paras);

      if(temp_p->get_IsNeutrino() ){ 
        vars.m_reco_slice_num++;
        vars.m_reco_slice_nuscore.push_back(temp_p->get_NuScore());

      }
    }

    //resize slice variables size;
    //    this->ResizeSlices(vars.m_reco_slice_num); 
  }
}
