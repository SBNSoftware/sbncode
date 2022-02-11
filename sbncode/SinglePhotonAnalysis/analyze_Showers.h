#include "SinglePhoton_module.h"
#include "reco_truth_matching.h"

namespace single_photon
{
    void SinglePhoton::ClearShowers(){
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

        m_reco_shower_start_to_nearest_dead_wire_plane0.clear();
        m_reco_shower_start_to_nearest_dead_wire_plane1.clear();
        m_reco_shower_start_to_nearest_dead_wire_plane2.clear();

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

    void SinglePhoton::ResizeShowers(size_t size){
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

        m_reco_shower_start_to_nearest_dead_wire_plane0.resize(size);
        m_reco_shower_start_to_nearest_dead_wire_plane1.resize(size);
        m_reco_shower_start_to_nearest_dead_wire_plane2.resize(size);

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

    void SinglePhoton::CreateShowerBranches(){
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

        vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane0",&m_reco_shower_start_to_nearest_dead_wire_plane0);
        vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane1",&m_reco_shower_start_to_nearest_dead_wire_plane1);
        vertex_tree->Branch("reco_shower_start_to_nearest_dead_wire_plane2",&m_reco_shower_start_to_nearest_dead_wire_plane2);

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

    void SinglePhoton::AnalyzeShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap, std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > & pfParticleToClusterMap,std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> >  & clusterToHitMap , 
            std::map<int, double>& sliceIdToNuScoreMap,
            std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
            std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap, 
            std::map<art::Ptr<recob::PFParticle>,bool> &PFPToNuSliceMap, 
            std::map<art::Ptr<recob::PFParticle>,double> &PFPToTrackScoreMap,
            PFParticleIdMap &pfParticleMap,
            std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>>& PFPtoShowerReco3DMap
            ){

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Begininning recob::Shower analysis suite"<<std::endl;;

        m_reco_asso_showers=showers.size();
        int i_shr = 0;
        this->ResizeShowers(m_reco_asso_showers);

        for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
        {

            const art::Ptr<recob::Shower> shower = *iter;
            const art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap[shower];

            art::Ptr<recob::Shower> shower3d;
            if(PFPtoShowerReco3DMap.count(pfp)==0){
                //std::cout<<"SHOWER_3D : <<ERROR!! No shower3d in map for this pfp"<<std::endl;
                //std::cout<<"Reverting to normal recob::Shower"<<std::endl;
                m_reco_shower3d_exists[i_shr] = 0;
                shower3d = shower;
            }else{
                shower3d   = PFPtoShowerReco3DMap[pfp];
                m_reco_shower3d_exists[i_shr] = 1;
            }
            const std::vector<art::Ptr<recob::Hit>> hits =  pfParticleToHitMap[pfp];
            const std::vector<art::Ptr<recob::Cluster>> clusters = pfParticleToClusterMap[pfp];


            //int m_shrid = shower->ID(); This is an used variable, always -999
            double m_length = shower->Length();
            double m_open_angle = shower->OpenAngle();

            TVector3 shr_start = shower->ShowerStart();
            TVector3 shr_dir = shower->Direction();

            TVector3 shr3d_start = shower3d->ShowerStart();
            TVector3 shr3d_dir = shower3d->Direction();

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t On Shower: "<<i_shr<<" which has length: "<<m_length<<""<<std::endl;;

            m_reco_shower_startx[i_shr] = shr_start.X();
            m_reco_shower_starty[i_shr] = shr_start.Y();
            m_reco_shower_startz[i_shr] = shr_start.Z();

 
            std::vector<double> hstart = {m_reco_shower_startx[i_shr],m_reco_shower_starty[i_shr],m_reco_shower_startz[i_shr]};
            m_reco_shower_start_dist_to_active_TPC[i_shr] = distToTPCActive(hstart);
            m_reco_shower_start_in_SCB[i_shr] = this->distToSCB(m_reco_shower_start_dist_to_SCB[i_shr],hstart);

            m_reco_shower_dirx[i_shr] = shr_dir.X();
            m_reco_shower_diry[i_shr] = shr_dir.Y();
            m_reco_shower_dirz[i_shr] = shr_dir.Z();
            m_reco_shower_length[i_shr] = m_length;
            m_reco_shower_openingangle[i_shr] = m_open_angle;

            m_reco_shower3d_startx[i_shr] = shr3d_start.X();
            m_reco_shower3d_starty[i_shr] = shr3d_start.Y();
            m_reco_shower3d_startz[i_shr] = shr3d_start.Z();
            m_reco_shower3d_dirx[i_shr] = shr3d_dir.X();
            m_reco_shower3d_diry[i_shr] = shr3d_dir.Y();
            m_reco_shower3d_dirz[i_shr] = shr3d_dir.Z();
            m_reco_shower3d_length[i_shr] = shower3d->Length();
            m_reco_shower3d_openingangle[i_shr] = shower3d->OpenAngle();


            m_reco_shower_conversion_distance[i_shr] = sqrt( pow(shr_start.X()-m_vertex_pos_x,2)+pow(shr_start.Y()-m_vertex_pos_y,2)+ pow(shr_start.Z()-m_vertex_pos_z,2)  );
            m_reco_shower3d_conversion_distance[i_shr] = sqrt( pow(shr3d_start.X()-m_vertex_pos_x,2)+pow(shr3d_start.Y()-m_vertex_pos_y,2)+ pow(shr3d_start.Z()-m_vertex_pos_z,2)  );

            //pandroa shower
            std::vector<double> shr_ts = {shr_start.X(), shr_start.Y(), shr_start.Z()};
            std::vector<double> shr_te = {shr_start.X()-shr_dir.X(),shr_start.Y()-shr_dir.Y(),shr_start.Z()-shr_dir.Z()};
            std::vector<double> shr_tv = {m_vertex_pos_x,m_vertex_pos_y,m_vertex_pos_z};

            m_reco_shower_impact_parameter[i_shr] = dist_line_point(shr_ts,shr_te,shr_tv );
            m_reco_shower_implied_dirx[i_shr] = shr_start.X()-m_vertex_pos_x;;
            m_reco_shower_implied_diry[i_shr] = shr_start.Y()-m_vertex_pos_y;
            m_reco_shower_implied_dirz[i_shr] = shr_start.Z()-m_vertex_pos_z;

            double norm = sqrt(pow(m_reco_shower_implied_dirx[i_shr],2)+pow(m_reco_shower_implied_diry[i_shr],2)+pow(m_reco_shower_implied_dirz[i_shr],2));
            m_reco_shower_implied_dirx[i_shr] = m_reco_shower_implied_dirx[i_shr]/norm;
            m_reco_shower_implied_diry[i_shr] = m_reco_shower_implied_diry[i_shr]/norm;
            m_reco_shower_implied_dirz[i_shr] = m_reco_shower_implied_dirz[i_shr]/norm;

            //now 3D shower
            std::vector<double> shr3d_ts = {shr3d_start.X(), shr3d_start.Y(), shr3d_start.Z()};
            std::vector<double> shr3d_te = {shr3d_start.X()-shr3d_dir.X(),shr3d_start.Y()-shr3d_dir.Y(),shr3d_start.Z()-shr3d_dir.Z()};
            std::vector<double> shr3d_tv = {m_vertex_pos_x,m_vertex_pos_y,m_vertex_pos_z};

            m_reco_shower3d_impact_parameter[i_shr] = dist_line_point(shr3d_ts,shr3d_te,shr3d_tv );
            m_reco_shower3d_implied_dirx[i_shr] = shr3d_start.X()-m_vertex_pos_x;;
            m_reco_shower3d_implied_diry[i_shr] = shr3d_start.Y()-m_vertex_pos_y;
            m_reco_shower3d_implied_dirz[i_shr] = shr3d_start.Z()-m_vertex_pos_z;

            double shr3d_norm = sqrt(pow(m_reco_shower3d_implied_dirx[i_shr],2)+pow(m_reco_shower3d_implied_diry[i_shr],2)+pow(m_reco_shower3d_implied_dirz[i_shr],2));
            m_reco_shower3d_implied_dirx[i_shr] = m_reco_shower3d_implied_dirx[i_shr]/shr3d_norm;
            m_reco_shower3d_implied_diry[i_shr] = m_reco_shower3d_implied_diry[i_shr]/shr3d_norm;
            m_reco_shower3d_implied_dirz[i_shr] = m_reco_shower3d_implied_dirz[i_shr]/shr3d_norm;



            m_reco_shower_theta_yz[i_shr] = atan2(m_reco_shower_diry[i_shr],m_reco_shower_dirz[i_shr]);
            m_reco_shower_phi_yx[i_shr] = atan2(m_reco_shower_diry[i_shr],m_reco_shower_dirx[i_shr]);

            m_reco_shower3d_theta_yz[i_shr] = atan2(m_reco_shower3d_diry[i_shr],m_reco_shower3d_dirz[i_shr]);
            m_reco_shower3d_phi_yx[i_shr] = atan2(m_reco_shower3d_diry[i_shr],m_reco_shower3d_dirx[i_shr]);


            m_reco_shower_start_to_nearest_dead_wire_plane0[i_shr] = distanceToNearestDeadWire(0, m_reco_shower_starty[i_shr], m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
            m_reco_shower_start_to_nearest_dead_wire_plane1[i_shr] = distanceToNearestDeadWire(1, m_reco_shower_starty[i_shr], m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
            m_reco_shower_start_to_nearest_dead_wire_plane2[i_shr] = distanceToNearestDeadWire(2, m_reco_shower_starty[i_shr], m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
            std::vector<int> t_num(3,0);   // num of triangles on each plane
            std::vector<int> t_numhits(3,0);  // num of hits on each plane
            std::vector<double> t_area(3,0.0);

            //Right, this basically loops over all hits in all planes and for each plane forms the Delaunay triangilization of it and calculates the 2D area inscribed by the convex hull
            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Starting Delaunay Triangleization"<<std::endl;;

            //auto start = std::chrono::high_resolution_clock::now();
            this->delaunay_hit_wrapper(hits, t_numhits, t_num, t_area);

            //auto finish = std::chrono::high_resolution_clock::now();
            //auto microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
            //if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Finished Delaunay Triangleization. It took "<< microseconds.count() << "ms and found "<<t_num[0]+t_num[1]+t_num[2]<<" triangles"<<std::endl;;

            m_reco_shower_delaunay_num_triangles_plane0[i_shr] = t_num[0];
            m_reco_shower_delaunay_num_triangles_plane1[i_shr] = t_num[1];
            m_reco_shower_delaunay_num_triangles_plane2[i_shr] = t_num[2];

            m_reco_shower_delaunay_area_plane0[i_shr] = t_area[0];
            m_reco_shower_delaunay_area_plane1[i_shr] = t_area[1];
            m_reco_shower_delaunay_area_plane2[i_shr] = t_area[2];

            m_reco_shower_num_hits_plane0[i_shr] = t_numhits[0];
            m_reco_shower_num_hits_plane1[i_shr] = t_numhits[1];
            m_reco_shower_num_hits_plane2[i_shr] = t_numhits[2];
            //-------------- Calorimetry 3D --------------------


            const std::vector< double > shr3d_energy = shower3d->Energy();
            const std::vector< double > shr3d_dEdx = shower3d->dEdx();
            //const int shr3d_bestplane = shower3d->best_plane();

 //           std::cout<<"SHOWER3D_ENERGY: best plane: "<<shr3d_bestplane<<std::endl;
            //for(auto &en:shr3d_energy){
            //    std::cout<<en<<" ";
            //}
            if(shr3d_energy.size()==3){
                m_reco_shower3d_energy_plane0[i_shr] = shr3d_energy[0];
                m_reco_shower3d_energy_plane1[i_shr] = shr3d_energy[1];
                m_reco_shower3d_energy_plane2[i_shr] = shr3d_energy[2];
            }else{
                m_reco_shower3d_energy_plane0[i_shr] =-99;
                m_reco_shower3d_energy_plane1[i_shr] =-99;
                m_reco_shower3d_energy_plane2[i_shr] =-999;
            }

   //         std::cout<<std::endl<<"SHOWER3D_DEDX: "<<std::endl;
            //for(auto &dedx: shr3d_dEdx){
            //    std::cout<<dedx<<" ";
            //}
            std::cout<<std::endl;
            if(shr3d_dEdx.size()==3){
                m_reco_shower3d_dEdx_plane0[i_shr] = shr3d_dEdx[0];
                m_reco_shower3d_dEdx_plane1[i_shr] = shr3d_dEdx[1];
                m_reco_shower3d_dEdx_plane2[i_shr] = shr3d_dEdx[2];
            }else{
                m_reco_shower3d_dEdx_plane0[i_shr] =-99;
                m_reco_shower3d_dEdx_plane1[i_shr] =-99;
                m_reco_shower3d_dEdx_plane2[i_shr] =-999;
            }


            //------------- calorimetry ------------

            m_reco_shower_energy_max[i_shr] = CalcEShower(hits);
            m_reco_shower_energy_plane0[i_shr] = CalcEShowerPlane(hits, 0);
            m_reco_shower_energy_plane1[i_shr] = CalcEShowerPlane(hits, 1);
            m_reco_shower_energy_plane2[i_shr] = CalcEShowerPlane(hits, 2);

            m_reco_shower_plane0_nhits[i_shr] = getNHitsPlane(hits, 0);
            m_reco_shower_plane1_nhits[i_shr] = getNHitsPlane(hits, 1);
            m_reco_shower_plane2_nhits[i_shr] = getNHitsPlane(hits, 2);

            m_reco_shower_plane0_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 0);
            m_reco_shower_plane1_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 1);
            m_reco_shower_plane2_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 2);


            //currently only run on 1 shower events
            if(showers.size()==1){
               for(auto &h: hits){ 

                    int plane= h->View();
                    int wire = h->WireID().Wire;
                    int tick = h->PeakTime();

                    m_reco_shower_hit_tick.push_back(tick);
                    m_reco_shower_hit_plane.push_back(plane);
                    m_reco_shower_hit_wire.push_back(wire);




               }

            }



            //std::cout<<"The energy on each plane is 0: "<< m_reco_shower_energy_plane0[i_shr]<<", 1: "<< m_reco_shower_energy_plane1[i_shr]<<", 2: "<<  m_reco_shower_energy_plane2[i_shr]<<std::endl;


            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 0"<<std::endl;
            m_reco_shower_dQdx_plane0[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 0 ); 
            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 1"<<std::endl;
            m_reco_shower_dQdx_plane1[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 1 ); 
            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 2"<<std::endl; 
            m_reco_shower_dQdx_plane2[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 2 ); 

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tstarting de/dx plane 0"<<std::endl;
            m_reco_shower_dEdx_plane0[i_shr] = CalcdEdxFromdQdx(m_reco_shower_dQdx_plane0[i_shr]);
            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tstarting de/dx plane 1"<<std::endl;
            m_reco_shower_dEdx_plane1[i_shr] = CalcdEdxFromdQdx(m_reco_shower_dQdx_plane1[i_shr]);
            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tstarting de/dx plane 2"<<std::endl;
            m_reco_shower_dEdx_plane2[i_shr] = CalcdEdxFromdQdx(m_reco_shower_dQdx_plane2[i_shr]);

            m_reco_shower_dEdx_plane0_median[i_shr] = getMedian(m_reco_shower_dEdx_plane0[i_shr]);
            m_reco_shower_dEdx_plane1_median[i_shr] = getMedian(m_reco_shower_dEdx_plane1[i_shr]);
            m_reco_shower_dEdx_plane2_median[i_shr] = getMedian(m_reco_shower_dEdx_plane2[i_shr]);

            m_reco_shower_angle_wrt_wires_plane0[i_shr] = getAnglewrtWires(shr_dir,0);
            m_reco_shower_angle_wrt_wires_plane1[i_shr] = getAnglewrtWires(shr_dir,1);
            m_reco_shower_angle_wrt_wires_plane2[i_shr] = getAnglewrtWires(shr_dir,2);


            m_reco_shower_dQdx_plane0_median[i_shr] = getMedian(m_reco_shower_dQdx_plane0[i_shr]);
            m_reco_shower_dQdx_plane1_median[i_shr] = getMedian(m_reco_shower_dQdx_plane1[i_shr]);
            m_reco_shower_dQdx_plane2_median[i_shr] = getMedian(m_reco_shower_dQdx_plane2[i_shr]);


            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tPAR: "<<m_reco_shower_dEdx_plane0[i_shr].size()<<" "<<m_reco_shower_dEdx_plane1[i_shr].size()<<" "<<m_reco_shower_dEdx_plane2[i_shr].size()<<std::endl;

            m_reco_shower_dEdx_plane0_mean[i_shr] = std::accumulate(m_reco_shower_dEdx_plane0[i_shr].begin(), m_reco_shower_dEdx_plane0[i_shr].end(), 0.0)/((double)m_reco_shower_dEdx_plane0[i_shr].size()); 
            m_reco_shower_dEdx_plane1_mean[i_shr] = std::accumulate(m_reco_shower_dEdx_plane1[i_shr].begin(), m_reco_shower_dEdx_plane1[i_shr].end(), 0.0)/((double)m_reco_shower_dEdx_plane1[i_shr].size()); 
            m_reco_shower_dEdx_plane2_mean[i_shr] = std::accumulate(m_reco_shower_dEdx_plane2[i_shr].begin(), m_reco_shower_dEdx_plane2[i_shr].end(), 0.0)/((double)m_reco_shower_dEdx_plane2[i_shr].size()); 

            auto maxp0 = std::max_element(m_reco_shower_dEdx_plane0[i_shr].begin(), m_reco_shower_dEdx_plane0[i_shr].end());
            auto maxp1 = std::max_element(m_reco_shower_dEdx_plane1[i_shr].begin(), m_reco_shower_dEdx_plane1[i_shr].end());
            auto maxp2 = std::max_element(m_reco_shower_dEdx_plane2[i_shr].begin(), m_reco_shower_dEdx_plane2[i_shr].end());
            auto minp0 = std::min_element(m_reco_shower_dEdx_plane0[i_shr].begin(), m_reco_shower_dEdx_plane0[i_shr].end());
            auto minp1 = std::min_element(m_reco_shower_dEdx_plane1[i_shr].begin(), m_reco_shower_dEdx_plane1[i_shr].end());
            auto minp2 = std::min_element(m_reco_shower_dEdx_plane2[i_shr].begin(), m_reco_shower_dEdx_plane2[i_shr].end());


            if(maxp0 == m_reco_shower_dEdx_plane0[i_shr].end()){
                m_reco_shower_dEdx_plane0_max[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane0_max[i_shr] = *maxp0; 
            }

            if(maxp1 == m_reco_shower_dEdx_plane1[i_shr].end()){
                m_reco_shower_dEdx_plane1_max[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane1_max[i_shr] = *maxp1; 
            }

            if(maxp2 == m_reco_shower_dEdx_plane2[i_shr].end()){
                m_reco_shower_dEdx_plane2_max[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane2_max[i_shr] = *maxp2; 
            }


            if(minp0 == m_reco_shower_dEdx_plane0[i_shr].end()){
                m_reco_shower_dEdx_plane0_min[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane0_min[i_shr] = *minp0; 
            }

            if(minp1 == m_reco_shower_dEdx_plane1[i_shr].end()){
                m_reco_shower_dEdx_plane1_min[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane1_min[i_shr] = *minp1; 
            }

            if(minp2 == m_reco_shower_dEdx_plane2[i_shr].end()){
                m_reco_shower_dEdx_plane2_min[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane2_min[i_shr] = *minp2; 
            }


            m_reco_shower_dEdx_plane0_nhits[i_shr] = m_reco_shower_dEdx_plane0[i_shr].size();
            m_reco_shower_dEdx_plane1_nhits[i_shr] = m_reco_shower_dEdx_plane1[i_shr].size();
            m_reco_shower_dEdx_plane2_nhits[i_shr] = m_reco_shower_dEdx_plane2[i_shr].size();

            m_reco_shower_dEdx_amalgamated[i_shr] = getAmalgamateddEdx( m_reco_shower_angle_wrt_wires_plane0[i_shr],  m_reco_shower_angle_wrt_wires_plane1[i_shr],  m_reco_shower_angle_wrt_wires_plane2[i_shr], m_reco_shower_dEdx_plane0_median[i_shr], m_reco_shower_dEdx_plane1_median[i_shr], m_reco_shower_dEdx_plane2_median[i_shr],m_reco_shower_dEdx_plane0_nhits[i_shr], m_reco_shower_dEdx_plane1_nhits[i_shr], m_reco_shower_dEdx_plane2_nhits[i_shr] );

            m_reco_shower_dEdx_amalgamated_nhits[i_shr] = getAmalgamateddEdxNHits(m_reco_shower_dEdx_amalgamated[i_shr], m_reco_shower_dEdx_plane0_median[i_shr], m_reco_shower_dEdx_plane1_median[i_shr], m_reco_shower_dEdx_plane2_median[i_shr],m_reco_shower_dEdx_plane0_nhits[i_shr], m_reco_shower_dEdx_plane1_nhits[i_shr], m_reco_shower_dEdx_plane2_nhits[i_shr] );

            //-------------- Flashes : Was there a flash in the beam_time and if so was it near in Z? --------------------
            double zmin = m_reco_shower_startz[i_shr];
            double zmax = zmin + m_reco_shower_dirz[i_shr]*m_reco_shower_length[i_shr];
            if(zmin > zmax) std::swap(zmin, zmax);

            double ymin = m_reco_shower_starty[i_shr];
            double ymax = zmin + m_reco_shower_diry[i_shr]*m_reco_shower_length[i_shr];
            if(ymin > ymax) std::swap(ymin, ymax);

            //Code property of Gray Yarbrough (all rights reserved)
            //int optical_flash_in_beamgate_counter=0;
            double shortest_dist_to_flash_z=DBL_MAX;
            double shortest_dist_to_flash_y=DBL_MAX;
            double shortest_dist_to_flash_yz=DBL_MAX;
            //-999 my nonsenese int can change
            int shortest_dist_to_flash_index_z=-999;
            int shortest_dist_to_flash_index_y=-999;
            int shortest_dist_to_flash_index_yz=-999;

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tnumber of flashes: "<< m_reco_num_flashes<< ""<<std::endl;;
            for(int i_flash = 0; i_flash < m_reco_num_flashes; ++i_flash) {

                double const zcenter=m_reco_flash_zcenter[i_flash];
                if(m_is_verbose) std::cout<< "SinglePhoton::AnalyzeShowers()\t||\tflash z center:" <<m_reco_flash_zcenter[i_flash]<< ""<<std::endl;;
                double const ycenter=m_reco_flash_ycenter[i_flash];
                if(m_is_verbose) std::cout<< "SinglePhoton::AnaluzeShowers()\t||\tflash y center:" <<m_reco_flash_ycenter[i_flash]<< ""<<std::endl;;

                //z plane
                double dist_z=DBL_MAX;
                if(zcenter < zmin) {
                    dist_z = zmin - zcenter;
                }
                else if(zcenter > zmax) {
                    dist_z = zcenter - zmax;
                }
                else {
                    dist_z = 0;
                }	    
                if(dist_z < shortest_dist_to_flash_z){
                    shortest_dist_to_flash_z = dist_z;
                    shortest_dist_to_flash_index_z=i_flash;
                }


                //y plane

                double dist_y=DBL_MAX;
                if(ycenter < ymin) {
                    dist_y = ymin - ycenter;
                }
                else if(ycenter > ymax) {
                    dist_y = ycenter - ymax;
                }
                else {
                    dist_y= 0;
                }	    
                if(dist_y < shortest_dist_to_flash_y){
                    shortest_dist_to_flash_y = dist_y;
                    shortest_dist_to_flash_index_y=i_flash;
                }

                double dist_yz=DBL_MAX;
                dist_yz=std::sqrt(dist_y*dist_y+dist_z*dist_z);
                if(dist_yz<shortest_dist_to_flash_yz){
                    shortest_dist_to_flash_yz = dist_yz;
                    shortest_dist_to_flash_index_yz=i_flash;
                }

            }


            //assume setting to nonsense value
            if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_z = -2;
            m_reco_shower_flash_shortest_distz[i_shr]=shortest_dist_to_flash_z;
            if(m_is_verbose) std::cout << "the shortest distance in z plane between a flash and the shower: " << shortest_dist_to_flash_z << ""<<std::endl;;
            m_reco_shower_flash_shortest_index_z[i_shr]=shortest_dist_to_flash_index_z;
            if(m_is_verbose) std::cout << "the index closest flash to shower z_plane: " << shortest_dist_to_flash_index_z << ""<<std::endl;;

            if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_y = -2;
            m_reco_shower_flash_shortest_disty[i_shr]=shortest_dist_to_flash_y;
            if(m_is_verbose) std::cout <<"the shortest distance in y plane between a flash and the shower: " << shortest_dist_to_flash_y << ""<<std::endl;;
            m_reco_shower_flash_shortest_index_y[i_shr]=shortest_dist_to_flash_index_y;
            if(m_is_verbose) std::cout << "the index closest flash to shower y_plane: " << shortest_dist_to_flash_index_y << ""<<std::endl;;

            if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_yz = -2;
            m_reco_shower_flash_shortest_distyz[i_shr]=shortest_dist_to_flash_yz;
            if(m_is_verbose) std::cout <<"the shortest distance in yz between a flash and the shower: " << shortest_dist_to_flash_yz << ""<<std::endl;;
            m_reco_shower_flash_shortest_index_yz[i_shr]=shortest_dist_to_flash_index_yz;
            if(m_is_verbose) std::cout << "the index closest flash to shower yz: " << shortest_dist_to_flash_index_yz << ""<<std::endl;;


            //end optical flash code


            m_reco_shower_num_daughters[i_shr] = pfp->NumDaughters();  //corresponding PFParticle
            if(m_reco_shower_num_daughters[i_shr]>0){
                //currently just look at 1 daughter
                m_reco_shower_daughter_trackscore[i_shr] = PFPToTrackScoreMap[pfParticleMap[pfp->Daughters().front()]];
            }


            //------------and finally some slice info-----------------

            m_reco_shower_sliceId[i_shr] = PFPToSliceIdMap[pfp];
            m_reco_shower_nuscore[i_shr] = sliceIdToNuScoreMap[ m_reco_shower_sliceId[i_shr]] ;
            m_reco_shower_isclearcosmic[i_shr] = PFPToClearCosmicMap[pfp];
            m_reco_shower_is_nuslice[i_shr] = PFPToNuSliceMap[pfp];
            //m_reco_shower_trackscore[i_shr] = PFPToTrackScoreMap[pfp];
            //std::cout<<"m_reco_shower_is_nuslice[i_shr] = "<<m_reco_shower_is_nuslice[i_shr]<<" for shr with pfp "<<pfp->Self()<<std::endl; 

            if ( PFPToTrackScoreMap.find(pfp) != PFPToTrackScoreMap.end() ) {
                m_reco_shower_trackscore[i_shr] = PFPToTrackScoreMap[pfp];
                m_reco_shower_pfparticle_pdg[i_shr] = pfp->PdgCode();
            } else{
                m_reco_shower_trackscore[i_shr] = -999; 
                m_reco_shower_pfparticle_pdg[i_shr] = -999;
            }

            if ( m_reco_shower_sliceId[i_shr] >0) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t On Shower: "<<i_shr<<". Pfp id = "<< pfp->Self()<<". The slice id for this shower is "<< m_reco_shower_sliceId[i_shr]<<", the neutrino score for this slice is "<< m_reco_shower_nuscore[i_shr]<<", and is_nuslice = "<<  m_reco_shower_is_nuslice[i_shr]<<". The track score is : "<< m_reco_shower_trackscore[i_shr]<<std::endl;





            i_shr++;
        }


        //Lets sort and order the showers
        m_reco_shower_ordered_energy_index = sort_indexes(m_reco_shower_energy_max);




        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Finished."<<std::endl;;
    }

    void SinglePhoton::AnalyzeKalmanShowers( const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> &showerToPFParticleMap,                         std::map<art::Ptr<recob::PFParticle>,art::Ptr<recob::Track>> &  pfParticlesToShowerKalmanMap, std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>>>&  kalmanTrackToCaloMap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap){

        std::cout<<"Singlephoton::AnalyzeKalmanShowerrs\t||\tStarting to Analyze Showers ("<<showers.size()<<") via Kalman "<<std::endl;

        std::vector<double> gains = {0,0,0};
        for(int plane =0; plane < 3; plane++){
            if (m_is_data == false &&  m_is_overlayed == false){
                gains[plane] = m_gain_mc[plane] ;
            } if (m_is_data == true ||  m_is_overlayed == true){
                gains[plane] = m_gain_data[plane] ;
            }
        }


        int i_shr=0;
        for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
        {
            const art::Ptr<recob::Shower> shower = *iter;
            const art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap[shower];
            std::vector<art::Ptr<recob::Hit>> hitz = pfParticleToHitMap[pfp];

            if( pfParticlesToShowerKalmanMap.count(pfp) == 0 ){
                std::cout<<"Singlephoton::AnalyzeKalmanShowerrs\t||\t Warning, no match for a Kalman track for this PFP."<<std::endl;
                continue;
            }
            const art::Ptr<recob::Track> kalman = pfParticlesToShowerKalmanMap[pfp];

            if(kalmanTrackToCaloMap.count(kalman)==0){
                std::cout<<"Singlephoton::AnalyzeKalmanShowerrs\t||\t Warning, no match for a Calo for this Kalman track."<<std::endl;
                m_reco_shower_kalman_exists[i_shr]=0;
                continue;
            }
            
            m_reco_shower_kalman_exists[i_shr]=1;
            const std::vector<art::Ptr<anab::Calorimetry>> calo = kalmanTrackToCaloMap[kalman];

            if(calo.size()!=3){
                std::cout<<"Singlephoton::AnalyzeKalmanShowerrs\t||\tERROR!! ERROR!!! anab::Calorimetery vector is not of length 3!!! "<<std::endl;
                continue;
            }
            //std::cout<<"index 0 is "<<calo[0]->PlaneID()<<", 1 is "<<calo[1]->PlaneID()<<", 2 is "<<calo[2]->PlaneID()<<std::endl;

            double res_range_lim = m_length_dqdx_box; //4cm 

            for(size_t p=0; p<calo.size();p++){
    
                int plane = calo[p]->PlaneID().Plane;
                if(plane<0 || plane > 3) continue; // Guanqun: plane == 3 is allowed??

                std::vector<double> t_dEdx; //in XX cm only  (4 for now)
                std::vector<double> t_res;
            

                for(size_t ix=0; ix<calo[p]->ResidualRange().size(); ++ix){

                    double rr = calo[p]->ResidualRange().back() - calo[p]->ResidualRange()[ix]; 
                    if(rr <= res_range_lim){
			// Guanqun: why is here a gains[plane], it's not converting ADC to E?
                        t_dEdx.push_back(gains[plane]*m_work_function*calo[p]->dQdx()[ix]*1e-6 /m_recombination_factor);
                        //t_dQdx.push_back(*calo[p]->dQdx()[x]);
                        t_res.push_back(rr);
                    }
                }
                /*std::cout<<"KAL: plane "<<calo[p]->PlaneID()<<" En: "<<std::endl;
                for(auto &xx: t_res) std::cout<<xx<<" ";
                std::cout<<std::endl;
                //for(auto &xx: t_dEdx) std::cout<<xx<<" ";
                std::cout<<std::endl<<" dQdx: "<<std::endl;
                //for(auto &xx: t_dQdx) std::cout<<xx<<" ";
                std::cout<<std::endl;
                */

                double tmean = NAN;
                double tmedian = NAN;


                if(t_dEdx.size()>0) tmedian = this->getMedian(t_dEdx);
                if(t_dEdx.size()>0) tmean = std::accumulate(t_dEdx.begin(), t_dEdx.end(), 0)/((double)t_dEdx.size());

                switch(plane){
                    case 0:
                            m_reco_shower_kalman_mean_dEdx_plane0[i_shr] = tmean ;
                            m_reco_shower_kalman_median_dEdx_plane0[i_shr] = tmedian ;
                            break;
                    case 1:
                            m_reco_shower_kalman_mean_dEdx_plane1[i_shr] = tmean;
                            m_reco_shower_kalman_median_dEdx_plane1[i_shr] = tmedian;
                            break;
                    case 2:
                            m_reco_shower_kalman_mean_dEdx_plane2[i_shr] = tmean;
                            m_reco_shower_kalman_median_dEdx_plane2[i_shr] = tmedian;
                            break;
                    default:
                            break;
                }


                const std::vector< anab::Point_t >  kal_pts = calo[p]->XYZ();   
                double circle = 1.0;//in cm
                std::vector<double> pts_within;
                std::vector<double> pts_x;

                for(size_t ix=0; ix<kal_pts.size(); ++ix){
                    //std::cout<<"KAL: "<<kal_pts[ix].X()<<" "<<kal_pts[ix].Y()<<" "<<kal_pts[ix].Z()<<std::endl;
                    pts_within.push_back(0);   
                    pts_x.push_back(calo[p]->ResidualRange().back()-calo[p]->ResidualRange()[ix]);

                    double wire = (double)calcWire(kal_pts[ix].Y(), kal_pts[ix].Z(), plane, m_TPC, m_Cryostat, *geom);
                    double time = calcTime(kal_pts[ix].X(), plane, m_TPC,m_Cryostat, *theDetector);

                    //loop over all hits  
                    for(auto &hit: hitz){
                        if(plane != hit->View())continue;
                        double this_w = (double)hit->WireID().Wire;
                        double this_t = (double)hit->PeakTime();
                        double dist = sqrt(pow(wire*0.3-this_w*0.3,2)+pow(time/25.0-this_t/25.0,2));
                        if(dist<=circle) pts_within.back()++;
                    }
                    //std::cout<<"KAL "<<ix<<" "<<pts_within.back()<<" "<<calo[p]->ResidualRange().back()-calo[p]->ResidualRange()[ix]<<std::endl;
                }
                if(false && pts_x.size()>2){
                TCanvas *c = new TCanvas();
                c->cd();
               
                TGraph *g  = new TGraph(pts_x.size(), &pts_x[0], &pts_within[0]);
                g->SetLineColor(kRed);
                g->SetLineWidth(2);
                g->Draw("alp");
                g->SetTitle(("kal_"+std::to_string(plane)+"_"+std::to_string(i_shr)+"_"+std::to_string(m_event_number) +".pdf").c_str());
                c->SaveAs(("kal_"+std::to_string(plane)+"_"+std::to_string(i_shr)+"_"+std::to_string(m_event_number) +".pdf").c_str(),"pdf");
                }
            }


            // some kalman averaging
            double tmp_kal_2 = m_reco_shower_kalman_median_dEdx_plane2[i_shr];
            double tmp_kal_1 = m_reco_shower_kalman_median_dEdx_plane1[i_shr];
            double tmp_kal_0 = m_reco_shower_kalman_median_dEdx_plane0[i_shr];
            double wei_0 = fabs(cos(m_reco_shower_angle_wrt_wires_plane0[i_shr]));
            double wei_1 = fabs(cos(m_reco_shower_angle_wrt_wires_plane1[i_shr]));
            double wei_2 = 20.0*fabs(cos(m_reco_shower_angle_wrt_wires_plane2[i_shr]));
           
            double thresh = 0.01;

            if(tmp_kal_2!=tmp_kal_2 || tmp_kal_2< thresh){
                tmp_kal_2 = 0;
                wei_2 = 0.0;
            }
            if(tmp_kal_1!=tmp_kal_1 || tmp_kal_1 < thresh){
                tmp_kal_1 = 0;
                wei_1 = 0.0;
            }
            if(tmp_kal_0!=tmp_kal_0 || tmp_kal_0 < thresh){
                tmp_kal_0 = 0;
                wei_0 = 0.0;
            }
            double kal_norm = wei_0+wei_1+wei_2;

            if(kal_norm!=0.0){
                m_reco_shower_kalman_median_dEdx_allplane[i_shr] = (tmp_kal_2*wei_2+tmp_kal_1*wei_1+tmp_kal_0*wei_0)/(kal_norm);
            }else{
                m_reco_shower_kalman_median_dEdx_allplane[i_shr] = NAN;
            }



            


        i_shr++;
        }
        return;
    }


    //-----------------------------------------------------------------------------------------------------------------------------------------
    void SinglePhoton::RecoMCShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  
            std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap, 
            std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,
            std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,
            std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector
            ){
        //OBSOLETE OBSOLETE


        if(m_is_verbose) std::cout<<"SinglePhoton::RecoMCShowers()\t||\t Begininning recob::Shower Reco-MC suite"<<std::endl;;

        int i_shr = 0;
        for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
        {
            const art::Ptr<recob::Shower> shower = *iter;
            m_sim_shower_matched[i_shr] = 0;
            if(showerToMCParticleMap.count(shower) > 0){



                const art::Ptr<simb::MCParticle> mcparticle = showerToMCParticleMap[shower];
                const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruthMap[mcparticle];

                std::vector<double> corrected(3);
                this->spacecharge_correction(mcparticle, corrected);
                m_sim_shower_matched[i_shr] = 1;
                m_sim_shower_energy[i_shr] = mcparticle->E();
                m_sim_shower_mass[i_shr] = mcparticle->Mass();
                m_sim_shower_kinetic_energy[i_shr] = mcparticle->E()-mcparticle->Mass();
                m_sim_shower_pdg[i_shr] = mcparticle->PdgCode();
                m_sim_shower_process[i_shr] = mcparticle->Process();
                m_sim_shower_start_x[i_shr] = corrected[0];
                m_sim_shower_start_y[i_shr] = corrected[1];
                m_sim_shower_start_z[i_shr] =corrected[2];

                m_sim_shower_origin[i_shr] = mctruth->Origin();
                //so this might be broken still due to mcparticle. must chcek
                if(mcparticle->Mother()>=(int)mcParticleVector.size()){
                    m_sim_shower_parent_pdg[i_shr] = -999;
                }else{
                    m_sim_shower_parent_pdg[i_shr] = mcParticleVector[mcparticle->Mother()]->PdgCode();
                }

                //OK is this photon matched to a delta?

                /*
                   std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t -- gamma ("<<mcparticle->TrackId()<<"| pdg: "<<mcparticle->PdgCode()<<") of status_code "<<mcparticle->StatusCode()<<std::endl;

                   art::Ptr<simb::MCParticle> nth_mother = crap_map[mcparticle->Mother()];
                   std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t ---- with mother "<<nth_mother->PdgCode()<<" ("<<nth_mother->TrackId()<<") status_code "<<nth_mother->StatusCode()<<std::endl;
                   int n_generation = 2;

                   while(nth_mother->StatusCode() != 0){

                   nth_mother = crap_map[nth_mother->Mother()]; 
                   std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t ---- and "<<n_generation<<"-mother "<<nth_mother->PdgCode()<<" ("<<nth_mother->TrackId()<<") and status_code "<<nth_mother->StatusCode()<<std::endl;
                   if( is_delta_map.count(nth_mother->PdgCode())>0 && nth_mother->StatusCode()==3){
                   std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t ------ Defintely From a Delta Decay! : "<<is_delta_map[nth_mother->PdgCode()]<<std::endl;

                   }
                   n_generation++;
                   }
                   */



            }
            i_shr++;
        }

    }


    void SinglePhoton::CollectCalo(const art::Event &evt, const art::Ptr<recob::Shower> &shower)
    {

    }

    double SinglePhoton::CalcEShowerPlane(const std::vector<art::Ptr<recob::Hit>>& hits, int this_plane){    
        double energy = 0.;

        //for each hit in the shower
        for (auto &thishitptr : hits){
            //check the plane
            int plane= thishitptr->View();

            //skip invalid planes     	
            if (plane != this_plane )	continue;

            //calc the energy of the hit
            double E = QtoEConversionHit(thishitptr, plane);	

            //add the energy to the plane
            energy += E;
        }//for each hit

        return energy;

    }

    double SinglePhoton::CalcEShower(const std::vector<art::Ptr<recob::Hit>> &hits){    
        double energy[3] = {0., 0., 0.};

        //std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t Looking at shower with "<<hits.size() <<" hits on all planes"<<std::endl;

        //for each hit in the shower
        for (auto &thishitptr : hits){
            //check the plane
            int plane= thishitptr->View();

            //skip invalid planes     	
            if (plane > 2 || plane < 0)	continue;

            //calc the energy of the hit
            double E = QtoEConversionHit(thishitptr, plane);    

            //add the energy to the plane
            energy[plane] += E;
        }//for each hiti

        //find the max energy on a single plane
        double max = energy[0];
        for (double en: energy){
            if( en > max){
                max = en;
            }
        }
        // std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The energy on each plane for this shower is "<<energy[0]<<", "<<energy[1]<<", "<<energy[2]<<std::endl;

        //return the highest energy on any of the planes
        return max;

    }

    double SinglePhoton::GetQHit(art::Ptr<recob::Hit> thishitptr, int plane){
        double gain;
        //choose gain based on whether data/mc and by plane
        if (m_is_data == false &&  m_is_overlayed == false){
            gain = m_gain_mc[plane] ;
            //if (m_is_verbose) std::cout<<"the gain for mc on plane "<<plane<<" is "<<gain<<std::endl;
        } if (m_is_data == true ||  m_is_overlayed == true){
            gain = m_gain_data[plane] ;
            //if (m_is_verbose) std::cout<<"the gain for data on plane "<<plane<<" is "<<gain<<std::endl;

        }

        double Q = thishitptr->Integral()*gain;
        return Q;
    }

    double SinglePhoton::QtoEConversionHit(art::Ptr<recob::Hit> thishitptr, int plane){
        return QtoEConversion(GetQHit(thishitptr, plane));

    }

    double SinglePhoton::QtoEConversion(double Q){
        //return the energy value converted to MeV (the factor of 1e-6)
        double E = Q* m_work_function *1e-6 /m_recombination_factor;
        return E;

    }


    std::vector<double> SinglePhoton::CalcdEdxFromdQdx(std::vector<double> dqdx){
        int n = dqdx.size();
        std::vector<double> dedx(n,0.0);
        for (int i = 0; i < n; i++){
            //std::cout<<"The dQ/dx is "<<dqdx[i]<<std::endl;
            dedx[i] = QtoEConversion(dqdx[i]);
            //std::cout<<"The dE/dx is "<<dedx[i]<<std::endl;
        }
        return dedx;
    }


    std::vector<double> SinglePhoton::CalcdQdxShower(
            const art::Ptr<recob::Shower>& shower,
            const std::vector<art::Ptr<recob::Cluster>> & clusters, 
            std::map<art::Ptr<recob::Cluster>,    std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,int plane){
        //if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The number of clusters in this shower is "<<clusters.size()<<std::endl;
        std::vector<double> dqdx;

        //get the 3D shower direction
        //note: in previous versions of the code there was a check to fix cases where the shower direction was inverted - this hasn't been implemented
        TVector3 shower_dir(shower->Direction().X(), shower->Direction().Y(),shower->Direction().Z());

        //calculate the pitch for this plane
        double pitch = getPitch(shower_dir, plane);	
        //if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The pitch between the shower and plane "<<plane<<" is "<<pitch<<std::endl;

        //for all the clusters in the shower
        for (const art::Ptr<recob::Cluster> &thiscluster: clusters){
            //keep only clusters on the plane
            if(thiscluster->View() != plane) continue;

            //calculate the cluster direction
            std::vector<double> cluster_axis = {cos(thiscluster->StartAngle()), sin(thiscluster->StartAngle())};		

            //get the cluster start and and in CM
            //std::cout<<"for plane/tpc/cryo:"<<plane<<"/"<<m_TPC<<"/"<<m_Cryostat<<", fXTicksOffset: "<<theDetector->GetXTicksOffset(plane, m_TPC, m_Cryostat)<<" fXTicksCoefficient: "<<theDetector->GetXTicksCoefficient(m_TPC, m_Cryostat)<<std::endl;

            //convert the cluster start and end positions to time and wire coordinates
			std::cout<<" CHECK "<<__LINE__<<" at "<<__FILE__<<" segmen fault warning"<<std::endl;
            std::vector<double> cluster_start = {0};//CHECK {thiscluster->StartWire() * m_wire_spacing,(thiscluster->StartTick() - theDetector->TriggerOffset())* _time2cm};
            std::vector<double> cluster_end = {0};//CHECK{thiscluster->EndWire() * m_wire_spacing,(thiscluster->EndTick() - theDetector->TriggerOffset())* _time2cm };

            //check that the cluster has non-zero length
            double length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) + pow(cluster_end[1] - cluster_start[1], 2));
            //if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The cluster length is "<<length<<std::endl;
            if (length <= 0){ 
                std::cout<<"skipping cluster on plane "<<plane<<", length = "<<length<<std::endl;
                continue;
            }


            //draw a rectangle around the cluster axis 
            std::vector<std::vector<double>> rectangle = buildRectangle(cluster_start, cluster_axis, m_width_dqdx_box, m_length_dqdx_box);	

            //get all the hits for this cluster
            std::vector<art::Ptr<recob::Hit>> hits =  clusterToHitMap[thiscluster];

            //for each hit in the cluster
            for (art::Ptr<recob::Hit> &thishit: hits){	
                //get the hit position in cm from the wire and time
                std::vector<double> thishit_pos ={0};//CHECK {thishit->WireID().Wire * m_wire_spacing, (thishit->PeakTime() - theDetector->TriggerOffset())* _time2cm};

                //check if inside the box
                bool v2 = isInsidev2(thishit_pos, rectangle);
                if (v2){
                    double q = GetQHit(thishit, plane); 
                    double this_dqdx = q/pitch; 
                    dqdx.push_back(this_dqdx);
                }//if hit falls inside the box

            }//for each hit inthe cluster
        }//for each cluster
        return dqdx;
    }

    double SinglePhoton::getPitch(TVector3 shower_dir, int plane){
        //get the wire direction for this plane - values are hardcoded which isn't great but the TPC geom object gave weird values
        TVector3 wire_dir = getWireVec(plane);

        //take dot product of shower and wire dir
        double cos = getCoswrtWires(shower_dir, wire_dir);

        //want only positive values so take abs, normalize by the lengths of the shower and wire
        cos = abs(cos)/(wire_dir.Mag() * shower_dir.Mag());	

        //If the cos is 0 shower is perpendicular and therefore get infinite distance 
        if (cos == 0){ return std::numeric_limits<double>::max(); }

        //output is always >= the wire spacing
        return m_wire_spacing/cos;
    }

    TVector3 SinglePhoton::getWireVec(int plane){
        TVector3 wire_dir;
        if (plane == 0){
            wire_dir = {0., -sqrt(3) / 2., 1 / 2.};
        } else if (plane == 1){
            wire_dir = {0., sqrt(3) / 2., 1 / 2.};
        } else if (plane == 2) {
            wire_dir = {0., 0., 1.};
        }
        return wire_dir;

    }

    double SinglePhoton::getCoswrtWires(TVector3 shower_dir, TVector3 wire_dir){
        //take the dot product between the wire direction and the shower direction
        double cos = wire_dir.Dot(shower_dir);
        return cos;
    }


    // shower_dir needs to be unit vector
    double SinglePhoton::getAnglewrtWires(TVector3 shower_dir,int plane){

        TVector3 wire_dir = getWireVec(plane);
        double cos_theta =  getCoswrtWires(shower_dir, wire_dir);

        double theta = acos(cos_theta);
        // return abs(theta);
        return abs(M_PI/2 - theta);

    }


    double SinglePhoton::getAmalgamateddEdx(double angle_wrt_plane0, double angle_wrt_plane1, double angle_wrt_plane2, double median_plane0, double median_plane1, double median_plane2, int plane0_nhits, int plane1_nhits, int plane2_nhits){
        //if the shower is within 10 degrees of the wires on plane 2, consider planes 1 and 0
        if(angle_wrt_plane2< degToRad(10)){
            //if it's too close to the wires on either of the planes, then stick with plane 2
            if (angle_wrt_plane1> degToRad(20)|| angle_wrt_plane0>degToRad(20) ){
                //but if it's outside of the range on plane 1, choose that
                if(angle_wrt_plane1> angle_wrt_plane0){
                    return median_plane1;
                } else{
                    return median_plane0;
                }
            }
        }
        if (plane2_nhits< 2){
            if (plane1_nhits >=2 ){           
                return median_plane1;
            } else if (plane0_nhits >=2 ){
                return median_plane0;
            }
        }

        return median_plane2;
    }

    int SinglePhoton::getAmalgamateddEdxNHits(double amalgamateddEdx, double median_plane0, double median_plane1, double median_plane2, int plane0_nhits, int plane1_nhits, int plane2_nhits){
        if (amalgamateddEdx == median_plane0){
            return plane0_nhits;
        }
        if (amalgamateddEdx == median_plane1){
            return plane1_nhits;
        }
        if (amalgamateddEdx == median_plane2){
            return plane2_nhits;
        }
        return -999;

    }   

    double SinglePhoton::degToRad(double deg){
        return deg * M_PI/180;
    }

    double SinglePhoton::radToDeg(double rad){
        return rad * 180/M_PI;
    }


    double SinglePhoton::getMeanHitWidthPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane){
        int nhits = 0;
        double widths = 0;
        for (art::Ptr<recob::Hit> thishitptr : hits){
            //check the plane
            int plane= thishitptr->View();

            //skip invalid planes       
            if (plane != this_plane) continue;

            widths += thishitptr->RMS(); // recob::Hit->RMS() returns RMS of the hit shape in tick units
            nhits++;


        }//for each hiti
        return   widths/(double)nhits;

    }



    int SinglePhoton::getNHitsPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane){
        int nhits = 0;
        for (art::Ptr<recob::Hit> thishitptr : hits){
            //check the plane
            int plane= thishitptr->View();

            //skip invalid planes       
            if (plane != this_plane) continue;

            nhits++;

        }//for each hiti
        return nhits;

    }

    std::vector<std::vector<double>> SinglePhoton::buildRectangle(std::vector<double> cluster_start, std::vector<double> cluster_axis, double width, double length){
        std::vector<std::vector<double>> corners;

        //get the axis perpedicular to the cluster axis
        double perp_axis[2] = {-cluster_axis[1], cluster_axis[0]};

        //create a vector for each corner of the rectangle on the plane
        //c1 = bottom left corner
        std::vector<double> c1 = {cluster_start[0] + perp_axis[0] * width / 2,  cluster_start[1] + perp_axis[1] * width / 2};
        //c2 = top left corner
        std::vector<double> c2 = {c1[0] + cluster_axis[0] * length, c1[1] + cluster_axis[1] * length};
        //c3 = bottom right corner
        std::vector<double> c3 = {cluster_start[0] - perp_axis[0] * width / 2, cluster_start[1] - perp_axis[1] * width / 2};
        //c4 = top right corner
        std::vector<double> c4 ={c3[0] + cluster_axis[0] * length, c3[1] + cluster_axis[1] * length}; 

        //save each of the vectors
        corners.push_back(c1);
        corners.push_back(c2);
        corners.push_back(c4);
        corners.push_back(c3);
        return corners;
    }

    bool SinglePhoton::insideBox(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle){
        //for a rectangle this is a known value but this is the most general
        int n_vertices = (int)rectangle.size();
        bool inside = false;
        int i, j = 0;
        //for each pair of vertices
        for (i = 0, j = n_vertices-1; i < n_vertices; j = i++) {
            //if the hit y coordinate is between the y and x coordinates of two vertices
            if ( ((rectangle[i][1]> thishit_pos[1]) != (rectangle[j][1]>thishit_pos[1])) 
                    &&(thishit_pos[0] < (rectangle[j][0]-rectangle[i][0]) * (thishit_pos[1]-rectangle[i][1]) / (rectangle[j][1]-rectangle[i][1]) + rectangle[i][0]) ){   
                if (inside == false){    
                    inside = true;
                } else{
                    inside = false;
                }
            }
        }
        return inside;
    }

    //determines if a point is inside the rectangle by summing the areas of the four triangles made by 
    //if the point is inside, the sum of the triangles should exactly equal the area of the rectangle
    //also returns true if the point is on the boundary
    bool SinglePhoton::isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle){
        int n_vertices = (int)rectangle.size();
        //bool inside = false;
        int i, j = 0;
        double areas = 0;
        //for each pair of vertices
        for (i = 0, j = n_vertices-1; i < n_vertices; j = i++) {
            //calculate the area of a triangle with the point and two vertices
            double this_area = areaTriangle(rectangle[i][0], rectangle[i][1], rectangle[j][0], rectangle[j][1], thishit_pos[0], thishit_pos[1]);
            areas += this_area;
        }        
        //calc area of the rectangle
        double area_rectangle = m_width_dqdx_box* m_length_dqdx_box;

        //check the sum of areas match
        if (abs(areas - area_rectangle) <= 0.001 ){
            return true;
        }
        return false;
    }

    //area of a triangle given three vertices
    double SinglePhoton::areaTriangle(double x1, double y1, double x2, double y2, double x3, double y3){
        double num = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2);
        return abs(num)/2;
    }

    double SinglePhoton::getMedian(std::vector<double> thisvector){
        //So return median if odd, average of median in even, if size==0, return the point. 
        
        //here the size corresponds to the max index
                
        int size = thisvector.size() - 1;
        //if no entries, return nonsense value
        if (size < 0) return NAN;
        if (size==0) return thisvector[size];

        //find index of median location
        double median;
        if (size%2 == 0){  // if vector has odd elements
            int ind = size/2;
            median = thisvector[ind];  
        } else{   // if vector has even number of elements
            int ind1 = size/2; 
            int ind2 = size/2-1;
            median = (thisvector[ind1]+thisvector[ind2])/2.0;
        }

        //return the value at median index
        return median;		
    }

}
