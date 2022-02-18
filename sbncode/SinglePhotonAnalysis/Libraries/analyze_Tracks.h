#include "SinglePhoton_module.h"
#include "TPrincipal.h"
#include "TVectorD.h"
#include "TruncMean.h"


namespace single_photon
{


    void SinglePhoton::ClearTracks(){
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

        m_reco_track_end_to_nearest_dead_wire_plane0.clear();
        m_reco_track_end_to_nearest_dead_wire_plane1.clear();
        m_reco_track_end_to_nearest_dead_wire_plane2.clear();

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

    void SinglePhoton::ResizeTracks(size_t size){
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

        m_reco_track_end_to_nearest_dead_wire_plane0.resize(size);
        m_reco_track_end_to_nearest_dead_wire_plane1.resize(size);
        m_reco_track_end_to_nearest_dead_wire_plane2.resize(size);

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

    void SinglePhoton::CreateTrackBranches(){
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

        vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane0",&m_reco_track_end_to_nearest_dead_wire_plane0);
        vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane1",&m_reco_track_end_to_nearest_dead_wire_plane1);
        vertex_tree->Branch("reco_track_end_to_nearest_dead_wire_plane2",&m_reco_track_end_to_nearest_dead_wire_plane2);

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






    void SinglePhoton::AnalyzeTracks(const std::vector<art::Ptr<recob::Track>>& tracks,
            std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToNuPFParticleMap,
            std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitsMap, 
            std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>>> & pfParticleToSpacePointsMap, 
            std::map<int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap,
            std::map<int, double> &sliceIdToNuScoreMap,
            std::map<art::Ptr<recob::PFParticle>,bool> &PFPToClearCosmicMap,
            std::map<art::Ptr<recob::PFParticle>, int> &PFPToSliceIdMap,
            std::map<art::Ptr<recob::PFParticle>,double> &PFPToTrackScoreMap,
            std::map<art::Ptr<recob::PFParticle>,bool> &PFPToNuSliceMap,
            PFParticleIdMap &pfParticleMap){


        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Starting recob::Track analysis"<<std::endl;;

        m_reco_asso_tracks = tracks.size();
        int i_trk=0;

        this->ResizeTracks(m_reco_asso_tracks);

        //const double adc2eU(5.1e-3);
        //const double adc2eV(5.2e-3);
        // const double adc2eW(5.4e-3);


        //    const double tau(theDetector->ElectronLifetime());


	//loop over each recob::Track
        for (TrackVector::const_iterator iter = tracks.begin(), iterEnd = tracks.end(); iter != iterEnd; ++iter)
        {

            const art::Ptr<recob::Track> track = *iter;
            const art::Ptr<recob::PFParticle> pfp = trackToNuPFParticleMap[track];
            const std::vector< art::Ptr<recob::SpacePoint> > trk_spacepoints = pfParticleToSpacePointsMap[pfp];
            const std::vector<art::Ptr<recob::Hit>> trk_hits  = pfParticleToHitsMap[pfp];

            int m_trkid = track->ID();
            double m_length = track->Length();
            auto m_trk_dir = track->Direction(); // type of m_trk_dir: a std::pair of two 3D vector
						//first: direction of first point, second: direction of the end of track

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t On Track: "<<i_trk<<" with TrackID: "<<m_trkid<<" and length: "<<m_length<<""<<std::endl;;

            m_reco_track_calo_energy_plane0[i_trk] = this->CalcEShowerPlane(trk_hits, 0);
            m_reco_track_calo_energy_plane1[i_trk] = this->CalcEShowerPlane(trk_hits, 1);
            m_reco_track_calo_energy_plane2[i_trk] = this->CalcEShowerPlane(trk_hits, 2);
            m_reco_track_calo_energy_max[i_trk] = std::max( m_reco_track_calo_energy_plane2[i_trk],  std::max(m_reco_track_calo_energy_plane0[i_trk],m_reco_track_calo_energy_plane1[i_trk]));

            m_reco_track_num_spacepoints[i_trk] = (int)trk_spacepoints.size();


            m_reco_track_startx[i_trk]= track->Start().X();   
            m_reco_track_starty[i_trk]= track->Start().Y();   
            m_reco_track_startz[i_trk]= track->Start().Z();   

            m_reco_track_length[i_trk] =m_length;
            m_reco_track_dirx[i_trk] = m_trk_dir.first.X();   
            m_reco_track_diry[i_trk] = m_trk_dir.first.Y();   
            m_reco_track_dirz[i_trk] = m_trk_dir.first.Z();   

            m_reco_track_endx[i_trk] = track->End().X();   
            m_reco_track_endy[i_trk]= track->End().Y();   
            m_reco_track_endz[i_trk]= track->End().Z();   
            
            std::vector<double> hend = {m_reco_track_endx[i_trk],m_reco_track_endy[i_trk],m_reco_track_endz[i_trk]};
            std::vector<double> hstart = {m_reco_track_startx[i_trk],m_reco_track_starty[i_trk],m_reco_track_startz[i_trk]};

            m_reco_track_end_dist_to_active_TPC[i_trk] = distToTPCActive(hend);
            m_reco_track_start_dist_to_active_TPC[i_trk] = distToTPCActive(hstart);

            m_reco_track_end_in_SCB[i_trk] = this->distToSCB(m_reco_track_end_dist_to_SCB[i_trk],hend);
            m_reco_track_start_in_SCB[i_trk] = this->distToSCB(m_reco_track_start_dist_to_SCB[i_trk],hstart);


            m_reco_track_theta_yz[i_trk] = atan2(m_reco_track_diry[i_trk],m_reco_track_dirz[i_trk]);
            m_reco_track_phi_yx[i_trk] = atan2(m_reco_track_diry[i_trk],m_reco_track_dirx[i_trk]);

            std::vector<double> tmp_trk_start = {m_reco_track_startx[i_trk],m_reco_track_starty[i_trk],m_reco_track_startz[i_trk]};
            std::vector<double> tmp_trk_end = {m_reco_track_endx[i_trk],m_reco_track_endy[i_trk],m_reco_track_endz[i_trk]};
            double max_dist_from_line = -9999999;

            m_reco_track_spacepoint_chi[i_trk] = 0.0;
            //Principal componant analysis of SPACEPOINTS and not trajectory points. Will add a few things here in future.

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Beginining PCA analysis of track"<<std::endl;;

            TPrincipal* principal = new TPrincipal(3,"ND");
            for(int x = 0; x < m_reco_track_num_spacepoints[i_trk]; x++){
		// get the position of spacepoint in xyz
                std::vector<double> tmp_spacepoints = {trk_spacepoints[x]->XYZ()[0],trk_spacepoints[x]->XYZ()[1] , trk_spacepoints[x]->XYZ()[2]};
                principal->AddRow(&tmp_spacepoints[0]);

		// distance between track direction and spacepoint
                double dist = dist_line_point(tmp_trk_start,tmp_trk_end,tmp_spacepoints);
                if(dist> max_dist_from_line) max_dist_from_line = dist;
                m_reco_track_spacepoint_chi[i_trk] += dist*dist;
            }
            m_reco_track_spacepoint_max_dist[i_trk]= max_dist_from_line;

            principal->MakePrincipals();
            TVectorD * eigen = (TVectorD*) principal->GetEigenValues();

            m_reco_track_spacepoint_principal0[i_trk]=(*eigen)(0);
            m_reco_track_spacepoint_principal1[i_trk]=(*eigen)(1);
            m_reco_track_spacepoint_principal2[i_trk]=(*eigen)(2);

            delete principal;

            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Completed PCA analysis of track. Primary componant: "<<m_reco_track_spacepoint_principal0.back()<<""<<std::endl;;

            //range based energy calculation assuming

            if(!m_run_pi0_filter){
		// assume this track is a proton track, get its energy
                m_reco_track_proton_kinetic_energy[i_trk] = proton_length2energy_tgraph.Eval(m_length)/1000.0; 
            }else{
                m_reco_track_proton_kinetic_energy[i_trk] = -9999;
            }

            if(m_length == 0.0) m_reco_track_proton_kinetic_energy[i_trk]=0.0;

            // Dead Wire Approximity
            m_reco_track_end_to_nearest_dead_wire_plane0[i_trk] = distanceToNearestDeadWire(0, m_reco_track_endy[i_trk], m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);
            m_reco_track_end_to_nearest_dead_wire_plane1[i_trk] = distanceToNearestDeadWire(1, m_reco_track_endy[i_trk], m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);
            m_reco_track_end_to_nearest_dead_wire_plane2[i_trk] = distanceToNearestDeadWire(2, m_reco_track_endy[i_trk], m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);

            m_reco_track_sliceId[i_trk] = PFPToSliceIdMap[pfp];
	    // Guanqun: how do you make sure the sliceId is positive, not -1, as for cosmic tracks?
	    // sliceIdToNuScoreMap seems to only have sliceID:nuScore pairs for these with actual nuScores.
            m_reco_track_nuscore[i_trk] = sliceIdToNuScoreMap[ m_reco_track_sliceId[i_trk]] ;
            m_reco_track_isclearcosmic[i_trk] = PFPToClearCosmicMap[pfp];

            //std::cout<<"checking track nuslice"<<std::endl;
           // std::cout<<"is nuslice for track with pfp "<<pfp->Self()<<" is: "<<PFPToNuSliceMap[pfp]<<std::endl;
            m_reco_track_is_nuslice[i_trk] = PFPToNuSliceMap[pfp];

            m_reco_track_num_daughters[i_trk] = pfp->NumDaughters();
            if(m_reco_track_num_daughters[i_trk]>0){
                //currently just look at 1 daughter
                m_reco_track_daughter_trackscore[i_trk] = PFPToTrackScoreMap[pfParticleMap[pfp->Daughters().front()]];
            }


          //  m_reco_track_trackscore[i_trk] = PFPToTrackScoreMap[pfp];
            if ( PFPToTrackScoreMap.find(pfp) != PFPToTrackScoreMap.end() ) {
                m_reco_track_trackscore[i_trk] = PFPToTrackScoreMap[pfp];
                m_reco_track_pfparticle_pdg[i_trk] = pfp->PdgCode();
            } else{
                m_reco_track_trackscore[i_trk] = -999; 
                m_reco_track_pfparticle_pdg[i_trk] = -999; 
            }

            //A loop over the trajectory points
            size_t const traj_size = track->CountValidPoints();
            m_reco_track_num_trajpoints[i_trk] = (int)traj_size;

            for(unsigned int  p = 0; p < traj_size; ++p) {
                //recob::Track::TrajectoryPoint_t const & trajp = track->TrajectoryPoint(j);
                //recob::Track::Point_t const & pos = trajp.position;
                //recob::Track::Vector_t const & mom = trajp.momentum;

            }


            i_trk++;

        } // end of recob::Track loop

        //Lets sort and order the showers
        m_reco_track_ordered_energy_index = sort_indexes(m_reco_track_proton_kinetic_energy);
        m_reco_track_ordered_displacement_index = sort_indexes(m_reco_track_length);


        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeTracks()\t||\t Finished."<<std::endl;;
    }

    void SinglePhoton::RecoMCTracks(const std::vector<art::Ptr<recob::Track>>& tracks,  
            std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap, 
            std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > & trackToMCParticleMap,
            std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,
            std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector,  
            std::map< int, art::Ptr<simb::MCParticle> > &      MCParticleToTrackIdMap, 
            std::map<int, double>& sliceIdToNuScoreMap,
            std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
            std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
            std::vector<double> & vfrac
            ){


        //if(m_is_verbose)           
            std::cout<<"SinglePhoton::RecoMCTracks()\t||\t Begininning recob::Track Reco-MC suite on: "<<tracks.size()<<" tracks."<<std::endl;

        int i_trk = 0;

        for(size_t k =0; k< tracks.size();++k){
            const art::Ptr<recob::Track> track = tracks[k];
            m_sim_track_matched[i_trk] = 0;

            if(trackToMCParticleMap.count(track)>0){

                const art::Ptr<simb::MCParticle> mcparticle = trackToMCParticleMap[track];
                std::cout<<"count2: "<<MCParticleToMCTruthMap.count(mcparticle)<<std::endl;
                const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruthMap[mcparticle];
                const art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap[track];

                std::vector<double> correctedstart(3);
                std::vector<double> correctedend(3);
                std::vector<double> raw_End  ={mcparticle->EndX(), mcparticle->EndY(), mcparticle->EndZ()};
               // std::cout<<"the raw end of this mcparticle is "<<raw_End[0]<<", "<<raw_End[1]<<", "<<raw_End[2]<<std::endl;
                this->spacecharge_correction(mcparticle, correctedstart);
                this->spacecharge_correction(mcparticle, correctedend, raw_End);
                
                //std::cout<<"the corrected end of this mcparticle is "<<correctedend[0]<<", "<<correctedend[1]<<", "<<correctedend[2]<<std::endl;

           
                m_sim_track_matched[i_trk] = 1;
                m_sim_track_energy[i_trk] = mcparticle->E();
                m_sim_track_mass[i_trk] = mcparticle->Mass();
                m_sim_track_kinetic_energy[i_trk] = m_sim_track_energy[i_trk]-m_sim_track_mass[i_trk];
                m_sim_track_pdg[i_trk] = mcparticle->PdgCode();
                m_sim_track_process[i_trk] = mcparticle->Process();
                m_sim_track_startx[i_trk] = correctedstart[0];
                m_sim_track_starty[i_trk] = correctedstart[1];
                m_sim_track_startz[i_trk] = correctedstart[2];
              
                m_sim_track_endx[i_trk]= correctedend[0];
                m_sim_track_endy[i_trk]= correctedend[1];
                m_sim_track_endz[i_trk]= correctedend[2];
            
                m_sim_track_length[i_trk]= sqrt(pow( m_sim_track_endx[i_trk] -  m_sim_track_startx[i_trk], 2)+ pow( m_sim_track_endy[i_trk] -  m_sim_track_starty[i_trk], 2) + pow( m_sim_track_endz[i_trk] -  m_sim_track_startz[i_trk], 2));
             
                m_sim_track_px[i_trk]=  mcparticle->Px();
                m_sim_track_py[i_trk]=  mcparticle->Py();
                m_sim_track_pz[i_trk]=  mcparticle->Pz();
             

                m_sim_track_origin[i_trk] = mctruth->Origin();
                m_sim_track_trackID[i_trk] = mcparticle->TrackId();
                m_sim_track_overlay_fraction[i_trk] = vfrac[i_trk];

                m_sim_track_sliceId[i_trk] = PFPToSliceIdMap[pfp];
                m_sim_track_nuscore[i_trk] = sliceIdToNuScoreMap[ m_sim_track_sliceId[i_trk]] ;
                m_sim_track_isclearcosmic[i_trk] = PFPToClearCosmicMap[pfp]; 


                if(mcparticle->Mother()>=(int)mcParticleVector.size()){
                    m_sim_track_parent_pdg[i_trk] = -1;
                }else{
                    m_sim_track_parent_pdg[i_trk] = mcParticleVector[mcparticle->Mother()]->PdgCode();
                }



                }
                i_trk++;
            }

         return;
        }





        void SinglePhoton::AnalyzeTrackCalo(const std::vector<art::Ptr<recob::Track>> &tracks, std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>> > &trackToCaloMap)
        {

            if(m_is_verbose) std::cout<<"SinglePhoton::CollectCalo(recob::Track)\t||\t Starting calo module for recob::Track"<<std::endl;;

            for(size_t i_trk = 0; i_trk<tracks.size(); ++i_trk){
                const art::Ptr<recob::Track>      track = tracks[i_trk];

                if(trackToCaloMap[track].size()!=3){
                    std::cout<<"Singlephoton::Tracks\t||\tERROR!! ERROR!!! anab::Calorimetery vector is not of length 3!!! "<<trackToCaloMap[track].size()<<". Skipping this track!!"<<std::endl;
                    continue;
                }

                const art::Ptr<anab::Calorimetry> calo_p0  = trackToCaloMap[track][0];
                const art::Ptr<anab::Calorimetry> calo_p1  = trackToCaloMap[track][1];
                const art::Ptr<anab::Calorimetry> calo_p2  = trackToCaloMap[track][2];


                size_t calo_length_p0 = calo_p0->dEdx().size();
                size_t calo_length_p1 = calo_p1->dEdx().size();
                size_t calo_length_p2 = calo_p2->dEdx().size();

                TruncMean tm_p0;
                TruncMean tm_p1;
                TruncMean tm_p2;

                std::vector<double> trunc_dEdx_p0;
                std::vector<double> res_range_good_p0;
                std::vector<double> dEdx_good_p0;

                std::vector<double> trunc_dEdx_p1;
                std::vector<double> res_range_good_p1;
                std::vector<double> dEdx_good_p1;

                std::vector<double> trunc_dEdx_p2;
                std::vector<double> res_range_good_p2;
                std::vector<double> dEdx_good_p2;

                m_reco_track_num_calo_hits_p0[i_trk] = (int)calo_length_p0;
                m_reco_track_num_calo_hits_p1[i_trk] = (int)calo_length_p1;
                m_reco_track_num_calo_hits_p2[i_trk] = (int)calo_length_p2;

                m_reco_track_best_calo_plane[i_trk]=-1;

		// guanqun: vectors have been clear and resized, so probably not need to reset their values?
                m_reco_track_good_calo_p0[i_trk] =  0;
                m_reco_track_mean_dEdx_p0[i_trk] =  0.0;
                m_reco_track_mean_dEdx_start_half_p0[i_trk] =  0.0;
                m_reco_track_mean_dEdx_end_half_p0[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_p0[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_start_half_p0[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_end_half_p0[i_trk] =  0.0;
                m_reco_track_trunc_PIDA_p0[i_trk] =  0.0;

                m_reco_track_good_calo_p1[i_trk] =  0;
                m_reco_track_mean_dEdx_p1[i_trk] =  0.0;
                m_reco_track_mean_dEdx_start_half_p1[i_trk] =  0.0;
                m_reco_track_mean_dEdx_end_half_p1[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_p1[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_start_half_p1[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_end_half_p1[i_trk] =  0.0;
                m_reco_track_trunc_PIDA_p1[i_trk] =  0.0;

                m_reco_track_good_calo_p2[i_trk] =  0;
                m_reco_track_mean_dEdx_p2[i_trk] =  0.0;
                m_reco_track_mean_dEdx_start_half_p2[i_trk] =  0.0;
                m_reco_track_mean_dEdx_end_half_p2[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_p2[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_start_half_p2[i_trk] =  0.0;
                m_reco_track_mean_trunc_dEdx_end_half_p2[i_trk] =  0.0;
                m_reco_track_trunc_PIDA_p2[i_trk] =  0.0;



                //First off look over ALL points
                //--------------------------------- plane 0 ----------- Induction
                for (size_t k = 0; k < calo_length_p0; ++k) {
                    double res_range =    calo_p0->ResidualRange()[k];  //ResidualRange() returns range from end of track
                    double dEdx =         calo_p0->dEdx()[k];

                    m_reco_track_mean_dEdx_p0[i_trk] += dEdx;
                    if(k <= calo_length_p0/2){ 
                        m_reco_track_mean_dEdx_start_half_p0[i_trk]+=dEdx;
                    }else{
                        m_reco_track_mean_dEdx_end_half_p0[i_trk]+=dEdx;
                    }

                    bool is_sensible = dEdx < m_track_calo_max_dEdx; 
                    bool is_nan =dEdx != dEdx; // != has higher precedence than = 
                    bool is_inf = std::isinf(dEdx);
                    bool is_nonzero = dEdx> m_track_calo_min_dEdx;

                    if(is_sensible && !is_nan && !is_inf && is_nonzero && k != 0 && k != calo_length_p0-1){
                        res_range_good_p0.push_back(res_range);
                        dEdx_good_p0.push_back(dEdx);
                    }

                    //    std::cout<<"\t"<<k<<" "<<calo->dEdx()[k]<<" "<<calo->ResidualRange()[k]<<" "<< ""<<std::endl;;
                }// End of first loop.

                m_reco_track_good_calo_p0[i_trk] = 0;
                if(res_range_good_p0.size() >= m_track_calo_min_dEdx_hits){
                    m_reco_track_good_calo_p0[i_trk] = res_range_good_p0.size();

                    //The radius we truncate over is going to be the max of either 1/frac of a track or 2x the minimum_dx in the res_range
                    double tenth_track = std::max(res_range_good_p0.front(), res_range_good_p0.back())/m_track_calo_trunc_fraction;
                    double min_dx = 999;
                    for(int j = res_range_good_p0.size()-1; j>1; j--){
                        double dx = fabs(res_range_good_p0[j]-res_range_good_p0[j-1]);
                        if(dx < min_dx) min_dx = dx;
                    }
                    double rad = std::max( min_dx*2, tenth_track); 

                    //Calculate the residual range
                    tm_p0.setRadius(rad);
                    tm_p0.CalcTruncMeanProfile(res_range_good_p0,dEdx_good_p0, trunc_dEdx_p0);			

                    double pida_sum_trunc=0.0;
                    //Calculate the mean truncated mean dEdx
                    for(size_t k=0; k< trunc_dEdx_p0.size(); k++){
                        double dEdx = trunc_dEdx_p0[k];
                        m_reco_track_mean_trunc_dEdx_p0[i_trk] += dEdx;
                        if(k <= trunc_dEdx_p0.size()/2){ 
                            m_reco_track_mean_trunc_dEdx_start_half_p0[i_trk]+=dEdx;
                        }else{
                            m_reco_track_mean_trunc_dEdx_end_half_p0[i_trk]+=dEdx;
                        }


                        if(trunc_dEdx_p0[k] != trunc_dEdx_p0[k] || std::isinf(trunc_dEdx_p0[k]) || trunc_dEdx_p0[k]<0){
                            std::cout<<"Truncated dedx is either inf or nan (or negative) @ "<<k<<" "<<trunc_dEdx_p0[k]<<std::endl;
                            std::cout<<"Vector Length : "<<trunc_dEdx_p0.size()<<std::endl;
                            std::cout<<"i\t range \t dedx \t trunc dedx"<<std::endl;
                            //for(int m=0; m<trunc_dEdx.size(); m++){
                            //    std::cout<<m<<"\t"<<c_resrange.at(m)<<"  "<<c_dEdx.at(m)<<"  "<<trunc_dEdx.at(m)<<std::endl;
                            //}
                            std::cout<<"Using Radius: "<<rad<<std::endl;
                            //exit(EXIT_FAILURE);
                            m_reco_track_good_calo_p0[i_trk] = 0; 
                        }
		
			// dEdx/pow(residual_range, -0.42) is the constant A in residual range formula
                        pida_sum_trunc += trunc_dEdx_p0[k]/(pow(res_range_good_p0[k],-0.42));
                    }
                    m_reco_track_trunc_PIDA_p0[i_trk] = pida_sum_trunc;           
                    m_reco_track_resrange_p0[i_trk] = res_range_good_p0;
                    m_reco_track_trunc_dEdx_p0[i_trk] = trunc_dEdx_p0;
                    m_reco_track_dEdx_p0[i_trk] = dEdx_good_p0;

                    //std::cout<<"the residual range at the start is "<<res_range_good[0]<<std::endl;
                }

                m_reco_track_mean_dEdx_p0[i_trk]            *=1.0/((double)calo_length_p0);
                m_reco_track_mean_dEdx_start_half_p0[i_trk] *=2.0/((double)calo_length_p0);
                m_reco_track_mean_dEdx_end_half_p0[i_trk]   *=2.0/((double)calo_length_p0);
                m_reco_track_mean_trunc_dEdx_p0[i_trk]            *=1.0/((double)trunc_dEdx_p0.size());
                m_reco_track_mean_trunc_dEdx_start_half_p0[i_trk] *=2.0/((double)trunc_dEdx_p0.size());
                m_reco_track_mean_trunc_dEdx_end_half_p0[i_trk]   *=2.0/((double)trunc_dEdx_p0.size());
                m_reco_track_trunc_PIDA_p0[i_trk]  *=1.0/((double)trunc_dEdx_p0.size());

                //First off look over ALL points
                //--------------------------------- plane 1 ----------- Induction
                for (size_t k = 0; k < calo_length_p1; ++k) {
                    double res_range =    calo_p1->ResidualRange()[k];
                    double dEdx =         calo_p1->dEdx()[k];

                    m_reco_track_mean_dEdx_p1[i_trk] += dEdx;
                    if(k <= calo_length_p1/2){ 
                        m_reco_track_mean_dEdx_start_half_p1[i_trk]+=dEdx;
                    }else{
                        m_reco_track_mean_dEdx_end_half_p1[i_trk]+=dEdx;
                    }

                    bool is_sensible = dEdx < m_track_calo_max_dEdx; 
                    bool is_nan =dEdx != dEdx; 
                    bool is_inf = std::isinf(dEdx);
                    bool is_nonzero = dEdx> m_track_calo_min_dEdx;

                    if(is_sensible && !is_nan && !is_inf && is_nonzero && k != 0 && k != calo_length_p1-1){
                        res_range_good_p1.push_back(res_range);
                        dEdx_good_p1.push_back(dEdx);
                    }

                    //    std::cout<<"\t"<<k<<" "<<calo->dEdx()[k]<<" "<<calo->ResidualRange()[k]<<" "<< ""<<std::endl;;
                }// End of first loop.

                m_reco_track_good_calo_p1[i_trk] = 0;
                if(res_range_good_p1.size() >= m_track_calo_min_dEdx_hits){
                    m_reco_track_good_calo_p1[i_trk] = res_range_good_p1.size();

                    //The radius we truncate over is going to be the max of either 1/frac of a track or 2x the minimum_dx in the res_range
                    double tenth_track = std::max(res_range_good_p1.front(), res_range_good_p1.back())/m_track_calo_trunc_fraction;
                    double min_dx = 999;
                    for(int j = res_range_good_p1.size()-1; j>1; j--){
                        double dx = fabs(res_range_good_p1[j]-res_range_good_p1[j-1]);
                        if(dx < min_dx) min_dx = dx;
                    }
                    double rad = std::max( min_dx*2, tenth_track); 

                    //Calculate the residual range
                    tm_p1.setRadius(rad);
                    tm_p1.CalcTruncMeanProfile(res_range_good_p1,dEdx_good_p1, trunc_dEdx_p1);			

                    double pida_sum_trunc=0.0;
                    //Calculate the mean truncated mean dEdx
                    for(size_t k=0; k< trunc_dEdx_p1.size(); k++){
                        double dEdx = trunc_dEdx_p1[k];
                        m_reco_track_mean_trunc_dEdx_p1[i_trk] += dEdx;
                        if(k <= trunc_dEdx_p1.size()/2){ 
                            m_reco_track_mean_trunc_dEdx_start_half_p1[i_trk]+=dEdx;
                        }else{
                            m_reco_track_mean_trunc_dEdx_end_half_p1[i_trk]+=dEdx;
                        }


                        if(trunc_dEdx_p1[k] != trunc_dEdx_p1[k] || std::isinf(trunc_dEdx_p1[k]) || trunc_dEdx_p1[k]<0){
                            std::cout<<"Truncated dedx is either inf or nan (or negative) @ "<<k<<" "<<trunc_dEdx_p1[k]<<std::endl;
                            std::cout<<"Vector Length : "<<trunc_dEdx_p1.size()<<std::endl;
                            std::cout<<"i\t range \t dedx \t trunc dedx"<<std::endl;
                            //for(int m=0; m<trunc_dEdx.size(); m++){
                            //    std::cout<<m<<"\t"<<c_resrange.at(m)<<"  "<<c_dEdx.at(m)<<"  "<<trunc_dEdx.at(m)<<std::endl;
                            //}
                            std::cout<<"Using Radius: "<<rad<<std::endl;
                            //exit(EXIT_FAILURE);
                            m_reco_track_good_calo_p1[i_trk] = 0; 
                        }

                        pida_sum_trunc += trunc_dEdx_p1[k]/(pow(res_range_good_p1[k],-0.42));
                    }
                    m_reco_track_trunc_PIDA_p1[i_trk] = pida_sum_trunc;           
                    m_reco_track_resrange_p1[i_trk] = res_range_good_p1;
                    m_reco_track_trunc_dEdx_p1[i_trk] = trunc_dEdx_p1;
                    m_reco_track_dEdx_p1[i_trk] = dEdx_good_p1;

                    //std::cout<<"the residual range at the start is "<<res_range_good[0]<<std::endl;
                }

                m_reco_track_mean_dEdx_p1[i_trk]            *=1.0/((double)calo_length_p1);
                m_reco_track_mean_dEdx_start_half_p1[i_trk] *=2.0/((double)calo_length_p1);
                m_reco_track_mean_dEdx_end_half_p1[i_trk]   *=2.0/((double)calo_length_p1);
                m_reco_track_mean_trunc_dEdx_p1[i_trk]            *=1.0/((double)trunc_dEdx_p1.size());
                m_reco_track_mean_trunc_dEdx_start_half_p1[i_trk] *=2.0/((double)trunc_dEdx_p1.size());
                m_reco_track_mean_trunc_dEdx_end_half_p1[i_trk]   *=2.0/((double)trunc_dEdx_p1.size());
                m_reco_track_trunc_PIDA_p1[i_trk]  *=1.0/((double)trunc_dEdx_p1.size());

                //First off look over ALL points
                //--------------------------------- plane 2 ----------- Collection
                for (size_t k = 0; k < calo_length_p2; ++k) {
                    double res_range =    calo_p2->ResidualRange()[k];
                    double dEdx =         calo_p2->dEdx()[k];

                    m_reco_track_mean_dEdx_p2[i_trk] += dEdx;
                    if(k <= calo_length_p2/2){ 
                        m_reco_track_mean_dEdx_start_half_p2[i_trk]+=dEdx;
                    }else{
                        m_reco_track_mean_dEdx_end_half_p2[i_trk]+=dEdx;
                    }

                    bool is_sensible = dEdx < m_track_calo_max_dEdx; 
                    bool is_nan =dEdx != dEdx; 
                    bool is_inf = std::isinf(dEdx);
                    bool is_nonzero = dEdx> m_track_calo_min_dEdx;

                    if(is_sensible && !is_nan && !is_inf && is_nonzero && k != 0 && k != calo_length_p2-1){
                        res_range_good_p2.push_back(res_range);
                        dEdx_good_p2.push_back(dEdx);
                    }

                    //    std::cout<<"\t"<<k<<" "<<calo->dEdx()[k]<<" "<<calo->ResidualRange()[k]<<" "<< ""<<std::endl;;
                }// End of first loop.

                m_reco_track_good_calo_p2[i_trk] = 0;
                if(res_range_good_p2.size() >= m_track_calo_min_dEdx_hits){
                    m_reco_track_good_calo_p2[i_trk] = res_range_good_p2.size();

                    //The radius we truncate over is going to be the max of either 1/frac of a track or 2x the minimum_dx in the res_range
                    double tenth_track = std::max(res_range_good_p2.front(), res_range_good_p2.back())/m_track_calo_trunc_fraction;
                    double min_dx = 999;
                    for(int j = res_range_good_p2.size()-1; j>1; j--){
                        double dx = fabs(res_range_good_p2[j]-res_range_good_p2[j-1]);
                        if(dx < min_dx) min_dx = dx;
                    }
                    double rad = std::max( min_dx*2, tenth_track); 

                    //Calculate the residual range
                    tm_p2.setRadius(rad);
                    tm_p2.CalcTruncMeanProfile(res_range_good_p2,dEdx_good_p2, trunc_dEdx_p2);			

                    double pida_sum_trunc=0.0;
                    //Calculate the mean truncated mean dEdx
                    for(size_t k=0; k< trunc_dEdx_p2.size(); k++){
                        double dEdx = trunc_dEdx_p2[k];
                        m_reco_track_mean_trunc_dEdx_p2[i_trk] += dEdx;
                        if(k <= trunc_dEdx_p2.size()/2){ 
                            m_reco_track_mean_trunc_dEdx_start_half_p2[i_trk]+=dEdx;
                        }else{
                            m_reco_track_mean_trunc_dEdx_end_half_p2[i_trk]+=dEdx;
                        }


                        if(trunc_dEdx_p2[k] != trunc_dEdx_p2[k] || std::isinf(trunc_dEdx_p2[k]) || trunc_dEdx_p2[k]<0){
                            std::cout<<"Truncated dedx is either inf or nan (or negative) @ "<<k<<" "<<trunc_dEdx_p2[k]<<std::endl;
                            std::cout<<"Vector Length : "<<trunc_dEdx_p2.size()<<std::endl;
                            std::cout<<"i\t range \t dedx \t trunc dedx"<<std::endl;
                            //for(int m=0; m<trunc_dEdx.size(); m++){
                            //    std::cout<<m<<"\t"<<c_resrange.at(m)<<"  "<<c_dEdx.at(m)<<"  "<<trunc_dEdx.at(m)<<std::endl;
                            //}
                            std::cout<<"Using Radius: "<<rad<<std::endl;
                            //exit(EXIT_FAILURE);
                            m_reco_track_good_calo_p2[i_trk] = 0; 
                        }

                        pida_sum_trunc += trunc_dEdx_p2[k]/(pow(res_range_good_p2[k],-0.42));
                    }
                    m_reco_track_trunc_PIDA_p2[i_trk] = pida_sum_trunc;           
                    m_reco_track_resrange_p2[i_trk] = res_range_good_p2;
                    m_reco_track_trunc_dEdx_p2[i_trk] = trunc_dEdx_p2;
                    m_reco_track_dEdx_p2[i_trk] = dEdx_good_p2;

                    //std::cout<<"the residual range at the start is "<<res_range_good[0]<<std::endl;
                }

                m_reco_track_mean_dEdx_p2[i_trk]            *=1.0/((double)calo_length_p2);
                m_reco_track_mean_dEdx_start_half_p2[i_trk] *=2.0/((double)calo_length_p2);
                m_reco_track_mean_dEdx_end_half_p2[i_trk]   *=2.0/((double)calo_length_p2);
                m_reco_track_mean_trunc_dEdx_p2[i_trk]            *=1.0/((double)trunc_dEdx_p2.size());
                m_reco_track_mean_trunc_dEdx_start_half_p2[i_trk] *=2.0/((double)trunc_dEdx_p2.size());
                m_reco_track_mean_trunc_dEdx_end_half_p2[i_trk]   *=2.0/((double)trunc_dEdx_p2.size());
                m_reco_track_trunc_PIDA_p2[i_trk]  *=1.0/((double)trunc_dEdx_p2.size());


                //************************ Now
                //lets pick one as a default?

                if(m_reco_track_good_calo_p2[i_trk]!=0){
                    m_reco_track_best_calo_plane[i_trk] = 2;

                    m_reco_track_mean_dEdx_best_plane[i_trk] = m_reco_track_mean_dEdx_p2[i_trk];
                    m_reco_track_mean_dEdx_start_half_best_plane[i_trk] = m_reco_track_mean_dEdx_start_half_p2[i_trk];
                    m_reco_track_mean_dEdx_end_half_best_plane[i_trk] = m_reco_track_mean_dEdx_end_half_p2[i_trk];
                    m_reco_track_good_calo_best_plane[i_trk] = m_reco_track_good_calo_p2[i_trk];
                    m_reco_track_trunc_dEdx_best_plane[i_trk] = m_reco_track_trunc_dEdx_p2[i_trk];
                    m_reco_track_mean_trunc_dEdx_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_p2[i_trk];
                    m_reco_track_mean_trunc_dEdx_start_half_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_start_half_p2[i_trk];
                    m_reco_track_mean_trunc_dEdx_end_half_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_end_half_p2[i_trk];
                    m_reco_track_trunc_PIDA_best_plane[i_trk] = m_reco_track_trunc_PIDA_p2[i_trk];
                    m_reco_track_resrange_best_plane[i_trk] = m_reco_track_resrange_p2[i_trk];
                    m_reco_track_dEdx_best_plane [i_trk] = m_reco_track_dEdx_p2[i_trk];



                }else if(m_reco_track_good_calo_p0[i_trk] > m_reco_track_good_calo_p1[i_trk] ){
                    m_reco_track_best_calo_plane[i_trk] = 0;

                    m_reco_track_mean_dEdx_best_plane[i_trk] = m_reco_track_mean_dEdx_p0[i_trk];
                    m_reco_track_mean_dEdx_start_half_best_plane[i_trk] = m_reco_track_mean_dEdx_start_half_p0[i_trk];
                    m_reco_track_mean_dEdx_end_half_best_plane[i_trk] = m_reco_track_mean_dEdx_end_half_p0[i_trk];
                    m_reco_track_good_calo_best_plane[i_trk] = m_reco_track_good_calo_p0[i_trk];
                    m_reco_track_trunc_dEdx_best_plane[i_trk] = m_reco_track_trunc_dEdx_p0[i_trk];
                    m_reco_track_mean_trunc_dEdx_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_p0[i_trk];
                    m_reco_track_mean_trunc_dEdx_start_half_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_start_half_p0[i_trk];
                    m_reco_track_mean_trunc_dEdx_end_half_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_end_half_p0[i_trk];
                    m_reco_track_trunc_PIDA_best_plane[i_trk] = m_reco_track_trunc_PIDA_p0[i_trk];
                    m_reco_track_resrange_best_plane[i_trk] = m_reco_track_resrange_p0[i_trk];
                    m_reco_track_dEdx_best_plane [i_trk] = m_reco_track_dEdx_p0[i_trk];


                }else if(m_reco_track_good_calo_p1[i_trk]!=0){
                    m_reco_track_best_calo_plane[i_trk] = 1;

                    m_reco_track_mean_dEdx_best_plane[i_trk] = m_reco_track_mean_dEdx_p1[i_trk];
                    m_reco_track_mean_dEdx_start_half_best_plane[i_trk] = m_reco_track_mean_dEdx_start_half_p1[i_trk];
                    m_reco_track_mean_dEdx_end_half_best_plane[i_trk] = m_reco_track_mean_dEdx_end_half_p1[i_trk];
                    m_reco_track_good_calo_best_plane[i_trk] = m_reco_track_good_calo_p1[i_trk];
                    m_reco_track_trunc_dEdx_best_plane[i_trk] = m_reco_track_trunc_dEdx_p1[i_trk];
                    m_reco_track_mean_trunc_dEdx_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_p1[i_trk];
                    m_reco_track_mean_trunc_dEdx_start_half_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_start_half_p1[i_trk];
                    m_reco_track_mean_trunc_dEdx_end_half_best_plane[i_trk] = m_reco_track_mean_trunc_dEdx_end_half_p1[i_trk];
                    m_reco_track_trunc_PIDA_best_plane[i_trk] = m_reco_track_trunc_PIDA_p1[i_trk];
                    m_reco_track_resrange_best_plane[i_trk] = m_reco_track_resrange_p1[i_trk];
                    m_reco_track_dEdx_best_plane [i_trk] = m_reco_track_dEdx_p1[i_trk];



                }else{
                    m_reco_track_best_calo_plane[i_trk] = -1; 
                }


            }
        }



        void SinglePhoton::CollectPID( std::vector<art::Ptr<recob::Track>> & tracks,
                std::map< art::Ptr<recob::Track>, art::Ptr<anab::ParticleID>> & trackToPIDMap){

            for(size_t i_trk=0; i_trk<tracks.size(); ++i_trk){
                art::Ptr<recob::Track> track = tracks[i_trk];
                art::Ptr<anab::ParticleID> pid = trackToPIDMap[track];
                if (!pid) {
                    std::cout << "[analyze_Tracks] bad PID object" << std::endl;
                    continue;
                }

                // For each PID object, create vector of PID scores for each algorithm
                // Loop over this and get scores for algorithm of choice
                // But first, prepare garbage values, just in case
                std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pid->ParticleIDAlgScores();
                double pidScore_BL_mu_plane0 = -999;
                double pidScore_BL_mu_plane1 = -999;
                double pidScore_BL_mu_plane2 = -999;
                double pidScore_BL_p_plane0 = -999;
                double pidScore_BL_p_plane1 = -999;
                double pidScore_BL_p_plane2 = -999;
                double pidScore_BL_mip_plane0 = -999;
                double pidScore_BL_mip_plane1 = -999;
                double pidScore_BL_mip_plane2 = -999;
                double pidScore_PIDA_plane0 = -999;
                double pidScore_PIDA_plane1 = -999;
                double pidScore_PIDA_plane2 = -999;
                double pidScore_chi2_mu_plane0 = -999;
                double pidScore_chi2_mu_plane1 = -999;
                double pidScore_chi2_mu_plane2 = -999;
                double pidScore_chi2_p_plane0 = -999;
                double pidScore_chi2_p_plane1 = -999;
                double pidScore_chi2_p_plane2 = -999;
                double pidScore_three_plane_proton = -999;

                //int planeid = 2;
                for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++) {

                    anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
                    //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);
                    int planeid = 0;//CHECK UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
                    /*
                    if (planeid != 2){
                        std::cout << "[ParticleIDValidation] Not using information for plane " 
                            << planeid << " (using plane 2 calorimetry only)" << std::endl;
                        continue;
                    }
                    */
                    if (AlgScore.fAlgName == "BraggPeakLLH" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood){
                        if (TMath::Abs(AlgScore.fAssumedPdg == 13)) {
                            if (planeid == 0) pidScore_BL_mu_plane0 = AlgScore.fValue; //AlgScore.fValue: value produced by algorithm
                            if (planeid == 1) pidScore_BL_mu_plane1 = AlgScore.fValue;
                            if (planeid == 2) pidScore_BL_mu_plane2 = AlgScore.fValue;
                        }
                        if (TMath::Abs(AlgScore.fAssumedPdg == 2212)) {
                            if (planeid == 0) pidScore_BL_p_plane0 = AlgScore.fValue;
                            if (planeid == 1) pidScore_BL_p_plane1 = AlgScore.fValue;
                            if (planeid == 2) pidScore_BL_p_plane2 = AlgScore.fValue;
                        }
                        if (TMath::Abs(AlgScore.fAssumedPdg == 0)) {
                            if (planeid == 0) pidScore_BL_mip_plane0 = AlgScore.fValue;
                            if (planeid == 1) pidScore_BL_mip_plane1 = AlgScore.fValue;
                            if (planeid == 2) pidScore_BL_mip_plane2 = AlgScore.fValue;
                        }
                    }
                    if (AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF){
                        if (TMath::Abs(AlgScore.fAssumedPdg == 13)) {
                            if (planeid == 0) pidScore_chi2_mu_plane0 = AlgScore.fValue;
                            if (planeid == 1) pidScore_chi2_mu_plane1 = AlgScore.fValue;
                            if (planeid == 2) pidScore_chi2_mu_plane2 = AlgScore.fValue;
                        }
                        if (TMath::Abs(AlgScore.fAssumedPdg == 2212)) {
                            if (planeid == 0) pidScore_chi2_p_plane0 = AlgScore.fValue;
                            if (planeid == 1) pidScore_chi2_p_plane1 = AlgScore.fValue;
                            if (planeid == 2) pidScore_chi2_p_plane2 = AlgScore.fValue;
                        }
                    }
                    if (AlgScore.fAlgName == "PIDA_mean" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
                        if (planeid == 0) pidScore_PIDA_plane0 = AlgScore.fValue;
                        if (planeid == 1) pidScore_PIDA_plane1 = AlgScore.fValue;
                        if (planeid == 2) pidScore_PIDA_plane2 = AlgScore.fValue;
                    }
//CHECK                    if (AlgScore.fAlgName == "ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && TMath::Abs(AlgScore.fAssumedPdg) == 2212 && AlgScore.fPlaneMask == UBPID::uB_SinglePlaneGetBitset(2) ){
                    if (AlgScore.fAlgName == "ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && TMath::Abs(AlgScore.fAssumedPdg) == 2212 && AlgScore.fPlaneMask == 0 ){
                    //if (AlgScore.fAlgName == "ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && TMath::Abs(AlgScore.fAssumedPdg) == 2212 && AlgScore.fPlaneMask == UBPID::uB_SinglePlaneGetBitset(2) && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
                        pidScore_three_plane_proton = AlgScore.fValue;
                    }
                }  // end looping over AlgScoresVec

                m_reco_track_pid_bragg_likelihood_mu_plane0[i_trk] = pidScore_BL_mu_plane0;
                m_reco_track_pid_bragg_likelihood_mu_plane1[i_trk] = pidScore_BL_mu_plane1;
                m_reco_track_pid_bragg_likelihood_mu_plane2[i_trk] = pidScore_BL_mu_plane2;
                m_reco_track_pid_bragg_likelihood_p_plane0[i_trk] = pidScore_BL_p_plane0;
                m_reco_track_pid_bragg_likelihood_p_plane1[i_trk] = pidScore_BL_p_plane1;
                m_reco_track_pid_bragg_likelihood_p_plane2[i_trk] = pidScore_BL_p_plane2;
                m_reco_track_pid_bragg_likelihood_mip_plane0[i_trk] = pidScore_BL_mip_plane0;
                m_reco_track_pid_bragg_likelihood_mip_plane1[i_trk] = pidScore_BL_mip_plane1;
                m_reco_track_pid_bragg_likelihood_mip_plane2[i_trk] = pidScore_BL_mip_plane2;
                m_reco_track_pid_chi2_mu_plane0[i_trk] = pidScore_chi2_mu_plane0;
                m_reco_track_pid_chi2_mu_plane1[i_trk] = pidScore_chi2_mu_plane1;
                m_reco_track_pid_chi2_mu_plane2[i_trk] = pidScore_chi2_mu_plane2;
                m_reco_track_pid_chi2_p_plane0[i_trk] = pidScore_chi2_p_plane0;
                m_reco_track_pid_chi2_p_plane1[i_trk] = pidScore_chi2_p_plane1;
                m_reco_track_pid_chi2_p_plane2[i_trk] = pidScore_chi2_p_plane2;
                m_reco_track_pid_pida_plane0[i_trk] = pidScore_PIDA_plane0;
                m_reco_track_pid_pida_plane1[i_trk] = pidScore_PIDA_plane1;
                m_reco_track_pid_pida_plane2[i_trk] = pidScore_PIDA_plane2;
                m_reco_track_pid_three_plane_proton_pid[i_trk] = pidScore_three_plane_proton;

            }
            return;
        }

    }
