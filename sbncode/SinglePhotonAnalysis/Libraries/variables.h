#ifndef SBNCODE_SINGLEPHOTONANALYSIS_VARIABLES_H
#define SBNCODE_SINGLEPHOTONANALYSIS_VARIABLES_H

#include "art/Framework/Principal/Handle.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "sbnobj/Common/Reco/OpT0FinderResult.h"

#include "TTree.h"
#include "TRandom3.h"

#include <vector>
#include <iostream>
#include <string>
#include <map>


//DEFINE GLOBAL VARIABLES, do it only once
namespace single_photon
{
      // name alias from pandora
      typedef art::ValidHandle< std::vector<recob::PFParticle> > PFParticleHandle;
      typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
      typedef std::vector< art::Ptr<recob::Track> > TrackVector;
      typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
      typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;

    extern bool g_is_verbose;

  struct para_all{
    std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
    std::vector<geo::BoxBoundedGeo> fActiveVolumes;
    double _time2cm;//value modeled from David's shower code

    bool s_use_PID_algorithms;
    bool s_use_delaunay;
    int  s_delaunay_max_hits;
    bool s_print_out_event;
    bool s_is_data; // value provided by pset
    bool s_is_overlayed;
    bool s_is_textgen;
    bool s_run_all_pfps;
    bool s_has_CRT;
    bool s_fill_trees;
    bool s_run_pi0_filter; //value provided by pset
    bool s_run_pi0_filter_2g1p;
    bool s_run_pi0_filter_2g0p;

    bool s_runPhotoNuTruth;
    bool s_runTrueEventweight;

    bool s_runSelectedEvent;  //if it should run only selected events
    bool s_runSEAview;
    bool s_runSEAviewStub;
    bool s_make_sss_plots;

    //SEAviwer bits
    double s_SEAviewPlotDistance;   //parameters related to shower-like object finding
    double s_SEAviewHitThreshold;
    double  s_SEAviewDbscanMinPts;
    double s_SEAviewDbscanEps;
    double s_SEAviewMaxPtsLinFit;
    bool   s_SEAviewMakePDF;
    int s_SEAviewNumRecoShower;
    int s_SEAviewNumRecoTrack;

    double s_SEAviewStubHitThreshold; //parameters related to track-like object finding
    double s_SEAviewStubPlotDistance;
    double s_SEAviewStubDbscanMinPts;
    double s_SEAviewStubDbscanEps;
    bool s_SEAviewStubMakePDF;
    int s_SEAviewStubNumRecoShower;
    int s_SEAviewStubNumRecoTrack;


    std::string s_selected_event_list; //full path for the file containing run/subrun/event number of selected events
    std::string s_shower3dLabel;
    std::string s_showerKalmanLabel;
    std::string s_showerKalmanCaloLabel;
    std::string s_pandoraLabel;         ///< The label for the pandora producer
    std::string s_trackLabel;           ///< The label for the track producer from PFParticles
    std::string s_showerLabel;          ///< The label for the shower producer from PFParticles
    std::string s_caloLabel;            ///< The label for calorimetry associations producer
    std::string s_flashLabel;
    std::string s_flashmatchLabel;
    std::string s_geantModuleLabel;
    //            std::string m_backtrackerLabel;
    std::string s_hitfinderLabel;
    std::string s_hitMCParticleAssnsLabel;
    std::string s_potLabel;
    std::string s_generatorLabel;

    std::string s_pidLabel;            ///< For PID stuff
    std::string s_CRTVetoLabel;
    std::string s_CRTTzeroLabel;
    std::string s_CRTHitProducer;
    std::string s_true_eventweight_label;

    std::string s_Spline_CV_label; //"eventweight4to4aFix"
    std::string s_truthmatching_signaldef;

    double s_max_conv_dist;
    double s_mass_pi0_mev;
    double s_beamgate_flash_start;
    double s_beamgate_flash_end;

    TRandom3 *rangen;

    double s_work_function;  //value provided by pset
    double s_recombination_factor; // value provided by pset

    std::vector<double> s_gain_mc; // value provided by pset 
    std::vector<double> s_gain_data; 
    double s_wire_spacing;

    //      int m_Cryostat;
    //      int m_TPC;

    bool   s_runCRT;
    double s_DTOffset;
    double s_Resolution;
    std::string  s_DAQHeaderProducer;//"daq"
    //      std::ofstream out_stream;

    spacecharge::SpaceCharge const * s_SCE;
    geo::GeometryCore const * s_geom;
    double s_tpc_active_XMin;
    double s_tpc_active_YMin;
    double s_tpc_active_ZMin;
    double s_tpc_active_XMax;
    double s_tpc_active_YMax;
    double s_tpc_active_ZMax;

    double s_width_dqdx_box; // value provided by pset
    double s_length_dqdx_box;

    double s_track_calo_min_dEdx;
    double s_track_calo_max_dEdx;
    double s_track_calo_min_dEdx_hits;
    double s_track_calo_trunc_fraction;


    double s_exiting_photon_energy_threshold ;
    double s_exiting_proton_energy_threshold ;


  };


  struct var_all{
    int pfp_w_bestnuID;
//    std::map<std::string, std::vector<double>> fmcweight;

    //Geant4
    std::vector<int> m_geant4_pdg;
    std::vector<int>          m_geant4_trackid;
    std::vector<int>          m_geant4_mother;
    std::vector<int>         m_geant4_statuscode;
    std::vector<double>          m_geant4_E;
    std::vector<double>          m_geant4_mass;
    std::vector<double>          m_geant4_px;
    std::vector<double>          m_geant4_py;
    std::vector<double>          m_geant4_pz;
    std::vector<double>          m_geant4_vx;
    std::vector<double>          m_geant4_vy;
    std::vector<double>          m_geant4_vz;
    std::vector<double>          m_geant4_dx;
    std::vector<double>          m_geant4_dy;
    std::vector<double>          m_geant4_dz;
    std::vector<std::string>          m_geant4_process;
    std::vector<std::string>          m_geant4_end_process;

    std::vector<double>          m_geant4_costheta;

    //------- Second shower related variables ----
    int m_sss_num_unassociated_hits; /* number of hits in the slice that're associated with neither shower nor tracks */
    int m_sss_num_unassociated_hits_below_threshold; /*number of unassociated hits that also didn't pass hit threshold,in the slice*/
    int m_sss_num_associated_hits; /* total number of hits from showers and tracks in the slice */


    //currently commenting this out for speed as its not used
    //ReadBDT * sssVetov1;

    int m_sss_num_candidates; /* number of unasso hit clusters which are not close enough to reco showers */
    std::vector<int> m_sss_candidate_in_nu_slice;
    std::vector<int> m_sss_candidate_num_hits;
    std::vector<int> m_sss_candidate_num_wires; //number of wires spanned by the candidate cluster
    std::vector<int>  m_sss_candidate_num_ticks;
    std::vector<int>  m_sss_candidate_plane; /* on which plan the unasso cluster is */
    std::vector<double> m_sss_candidate_PCA;
    std::vector<double> m_sss_candidate_mean_ADC;
    std::vector<double> m_sss_candidate_ADC_RMS;
    std::vector<double> m_sss_candidate_impact_parameter;
    std::vector<double> m_sss_candidate_fit_slope; //slope of the cluster direction
    std::vector<double> m_sss_candidate_veto_score;
    std::vector<double> m_sss_candidate_fit_constant; //intercept of the cluster direction
    std::vector<double> m_sss_candidate_mean_tick;
    std::vector<double> m_sss_candidate_max_tick;
    std::vector<double> m_sss_candidate_min_tick;
    std::vector<double> m_sss_candidate_min_wire;
    std::vector<double> m_sss_candidate_max_wire;
    std::vector<double> m_sss_candidate_mean_wire;
    std::vector<double> m_sss_candidate_min_dist;  // min distance from unasso cluter to the vertex */
    std::vector<double> m_sss_candidate_wire_tick_based_length;    //length of the cluster, calculated based on the wire & tick span of the cluster
    std::vector<double> m_sss_candidate_energy;
    std::vector<double> m_sss_candidate_angle_to_shower;
    std::vector<double> m_sss_candidate_closest_neighbour;
    std::vector<int>    m_sss_candidate_remerge; // index of the recob::shower candidate cluster is close to (expect it to be -1)
    std::vector<int>    m_sss_candidate_matched; /* has matched this unasso cluter to a primary MCParticle: 0-No, 1-Yes */
    std::vector<double> m_sss_candidate_matched_energy_fraction_best_plane; /* matched energy fraction of the best-matched MCParticle on best-plane */ 
    std::vector<int>    m_sss_candidate_pdg;   /* pdg of the matched MCParticle */
    std::vector<int>    m_sss_candidate_parent_pdg;
    std::vector<int>    m_sss_candidate_trackid; /* track ID of the matched MCParticle */
    std::vector<double> m_sss_candidate_true_energy;
    std::vector<double> m_sss_candidate_overlay_fraction; /* fraction of overlay in the unasso cluster hits */


    //------------ sss3d_showers variables are for reco::showers which are in the events, but not in the slice ----

    int m_sss3d_num_showers;  /* number of showers in the event but not in the slice */
    std::vector<double> m_sss3d_shower_start_x; /* shower start in X axis, for all showers in the event but not in the slice*/
    std::vector<double> m_sss3d_shower_start_y;
    std::vector<double> m_sss3d_shower_start_z;
    std::vector<double> m_sss3d_shower_dir_x; /* shower direction projection on X axis */
    std::vector<double> m_sss3d_shower_dir_y;
    std::vector<double> m_sss3d_shower_dir_z;
    std::vector<double> m_sss3d_shower_length;
    std::vector<double> m_sss3d_shower_conversion_dist; /* dist between shower start and vertex*/

    std::vector<double> m_sss3d_shower_invariant_mass; /* invariant mass of primary recob::shower, and each shower in the event, 
                              * calculated assuming vertex is where their mother particle decays */

    std::vector<double> m_sss3d_shower_implied_invariant_mass; /* similar to invariance mass, except this invariant mass  
                                  * is calced direclty using shower direction of two showers */

    std::vector<double> m_sss3d_shower_impact_parameter; /* dist between vertex and shower direction line */
    std::vector<double> m_sss3d_shower_ioc_ratio; /* ratio of impact parameter over conversion dist 
                           * 0 if the conversion distance is 0*/
    std::vector<double> m_sss3d_shower_energy_max; /* max energy of all planes (E summed from hits) */
    std::vector<double> m_sss3d_shower_score;
    std::vector<int> m_sss3d_slice_nu;
    std::vector<int> m_sss3d_slice_clear_cosmic;

    double m_sss3d_ioc_ranked_en;
    double m_sss3d_ioc_ranked_conv;
    double m_sss3d_ioc_ranked_invar;
    double m_sss3d_ioc_ranked_implied_invar;
    double m_sss3d_ioc_ranked_ioc;
    double m_sss3d_ioc_ranked_opang;
    double m_sss3d_ioc_ranked_implied_opang;
    int m_sss3d_ioc_ranked_id; //index of the sss3d_shower that has the smallest ioc.

    //same parameters, of the sss3d shower whose implied invariant mass together with primary recob::shower is closest to pi0 mass --
    double m_sss3d_invar_ranked_en;
    double m_sss3d_invar_ranked_conv;
    double m_sss3d_invar_ranked_invar;
    double m_sss3d_invar_ranked_implied_invar;
    double m_sss3d_invar_ranked_ioc;
    double m_sss3d_invar_ranked_opang;
    double m_sss3d_invar_ranked_implied_opang;
    int m_sss3d_invar_ranked_id;



    //--------------- sss2d showers are essentially group of cluters on 3 planes, that have the potential to be a shower -------
    //--------------- they are not recob::showers --------------------------

    // sss2d_ioc_ranked variables are the varaibles (mean shower energy, conv. dist, ioc, etc) of the sss2d shower that has the smallest ioc 
    // sss2d_conv_ranked variables are the varaibles (energy, conv. dist, ioc, etc) of the sss2d shower that has the smallest conv. distance 
    // sss2d_invar_ranked variables are the varaibles of the sss2d shower whose invariant mass together with primary shower is closest to pi0. 
    double m_sss2d_ioc_ranked_en;
    double m_sss2d_ioc_ranked_conv;
    double m_sss2d_ioc_ranked_ioc;
    double m_sss2d_ioc_ranked_pca;
    double m_sss2d_ioc_ranked_invar;
    double m_sss2d_ioc_ranked_angle_to_shower;
    int m_sss2d_ioc_ranked_num_planes;

    double m_sss2d_conv_ranked_en;
    double m_sss2d_conv_ranked_conv;
    double m_sss2d_conv_ranked_ioc;
    double m_sss2d_conv_ranked_pca;
    double m_sss2d_conv_ranked_invar;
    double m_sss2d_conv_ranked_angle_to_shower;
    int m_sss2d_conv_ranked_num_planes;

    double m_sss2d_invar_ranked_en;
    double m_sss2d_invar_ranked_conv;
    double m_sss2d_invar_ranked_ioc;
    double m_sss2d_invar_ranked_pca;
    double m_sss2d_invar_ranked_invar;
    double m_sss2d_invar_ranked_angle_to_shower;
    int m_sss2d_invar_ranked_num_planes;



    std::set<std::vector<int>> m_selected_set;  //set of selected events     


    //SSS parameters
    TTree* run_subrun_tree;
    TTree* pot_tree;
    TTree* vertex_tree;
    TTree* eventweight_tree;
    TTree* ncdelta_slice_tree;

    TTree* geant4_tree;

    TTree* true_eventweight_tree;

    //------------ POT related variables --------------
    int m_number_of_events;
    int m_number_of_events_in_subrun;
    double m_pot_count;
    int m_number_of_vertices;

    int m_run;
    int m_subrun;
    double m_subrun_pot;
    int m_subrun_counts;

    //------------ Event Related Variables -------------
    int m_run_number;
    int m_subrun_number;
    int m_event_number;
    double m_pot_per_event;
    double m_pot_per_subrun;

    int m_test_matched_hits;
    int m_reco_slice_objects;


    //------- Potential Unreconstructed Track Stub related variables ----
    int m_trackstub_num_unassociated_hits; /* number of hits in the slice that're associated with neither shower nor tracks */
    int m_trackstub_unassociated_hits_below_threshold; /*number of unassociated hits that also didn't pass hit threshold,in the slice*/
    int m_trackstub_associated_hits; /* total number of hits from showers and tracks in the slice */


    int m_trackstub_num_candidates; /* number of unasso hit clusters which are not close enough to reco showers */
    std::vector<int> m_trackstub_candidate_in_nu_slice; /* check if candidate is in neutrino slice: 1->YES, 0->Parts in neutrino slice, -1->Not at all */
    std::vector<int> m_trackstub_candidate_num_hits;
    std::vector<int> m_trackstub_candidate_num_wires; //number of wires spanned by the candidate cluster
    std::vector<int>  m_trackstub_candidate_num_ticks;
    std::vector<int>  m_trackstub_candidate_plane; /* on which plan the unasso cluster is */
    std::vector<double> m_trackstub_candidate_PCA;
    std::vector<double> m_trackstub_candidate_mean_ADC;
    std::vector<double> m_trackstub_candidate_ADC_RMS;
    std::vector<double> m_trackstub_candidate_veto_score;
    std::vector<double> m_trackstub_candidate_mean_tick;
    std::vector<double> m_trackstub_candidate_max_tick;
    std::vector<double> m_trackstub_candidate_min_tick;
    std::vector<double> m_trackstub_candidate_min_wire;
    std::vector<double> m_trackstub_candidate_max_wire;
    std::vector<double> m_trackstub_candidate_mean_wire;
    std::vector<double> m_trackstub_candidate_min_dist;  // min distance from unasso cluter to the vertex */
    std::vector<double> m_trackstub_candidate_min_impact_parameter_to_shower; //min impact parameter of all hits in cluster to the recob::shower direction line (on 2D plane)
    std::vector<double> m_trackstub_candidate_min_conversion_dist_to_shower_start;  //min distance between hits and recob::shower start (on 2D plane)
    std::vector<double> m_trackstub_candidate_min_ioc_to_shower_start;        //min ratio of impact_parameter_to_shower/conversion_dist_to_shower_start of all hits in the cluster
    std::vector<double> m_trackstub_candidate_ioc_based_length;    //length of the cluster, calculated based on the IOC of hit
    std::vector<double> m_trackstub_candidate_wire_tick_based_length;    //length of the cluster, calculated based on the wire & tick span of the cluster
    std::vector<double> m_trackstub_candidate_mean_ADC_first_half;    // mean ADC per hit for the first half of cluster (cluster divided into halves based on hit IOC)
    std::vector<double> m_trackstub_candidate_mean_ADC_second_half;
    std::vector<double> m_trackstub_candidate_mean_ADC_first_to_second_ratio; // ratio of the mean ADC per hit, first half of cluster over second half.
    std::vector<double> m_trackstub_candidate_track_angle_wrt_shower_direction;   //treat cluster as a track, angle between track direction and the shower direction
    std::vector<double> m_trackstub_candidate_linear_fit_chi2;    // chi2 from linear fit of the  {wire, tick} distribution of the cluster
    std::vector<double> m_trackstub_candidate_energy;
    std::vector<int>    m_trackstub_candidate_remerge; // index of the recob::shower candidate cluster is close to (expect it to be -1)
    std::vector<int>    m_trackstub_candidate_matched; /* has matched this unasso cluter to a primary MCParticle: 0-No, 1-Yes */
    std::vector<double> m_trackstub_candidate_matched_energy_fraction_best_plane; /* matched energy fraction of the best-matched MCParticle on best-plane */ 
    std::vector<int>    m_trackstub_candidate_pdg;   /* pdg of the matched MCParticle */
    std::vector<int>    m_trackstub_candidate_parent_pdg;
    std::vector<int>    m_trackstub_candidate_trackid; /* track ID of the matched MCParticle */
    std::vector<double> m_trackstub_candidate_true_energy;  /* true energy of the matched MCParticle */
    std::vector<double> m_trackstub_candidate_overlay_fraction; /* fraction of overlay in the unasso cluster hits */

    //------- grouped stub clusters --------------
    int m_trackstub_num_candidate_groups;           /* number of groups */ 
    std::vector<std::vector<double>> m_grouped_trackstub_candidate_indices; /* indices of stub clusters that are matched as a group */
    std::vector<double> m_trackstub_candidate_group_timeoverlap_fraction;   /* minimum fraction of the time overlap of grouped stub clusters */


    //-------------- Flash related variables -------------
    std::vector<double> m_reco_flashmatch_time;
    std::vector<double> m_reco_flashmatch_meas_pe;
    std::vector<double> m_reco_flashmatch_hypo_pe;
    std::vector<double> m_reco_flashmatch_score;

    std::vector<double> m_reco_flash_total_pe;
    std::vector<double> m_reco_flash_time;
    std::vector<double> m_reco_flash_time_width;
    std::vector<double> m_reco_flash_abs_time;
    std::vector<int>    m_reco_flash_frame;
    std::vector<double> m_reco_flash_ycenter;
    std::vector<double> m_reco_flash_ywidth;
    std::vector<double> m_reco_flash_zcenter;
    std::vector<double> m_reco_flash_zwidth;
    std::vector<double> m_reco_flash_total_pe_in_beamgate;
    std::vector<double> m_reco_flash_time_in_beamgate;
    std::vector<double> m_reco_flash_ycenter_in_beamgate;
    std::vector<double> m_reco_flash_zcenter_in_beamgate;

    double m_reco_flashmatch_time_bestscore;
    double m_reco_flashmatch_time_energetic;
    double m_reco_flash_time_energetic;
    double m_reco_flash_pe_peak;
    int m_reco_num_flashes;
    int m_reco_num_flashes_in_beamgate;
    //------------ Vertex Related variables -------------
    int m_reco_vertex_size;
    double m_vertex_pos_x;
    double m_vertex_pos_y;
    double m_vertex_pos_z;
    double m_vertex_pos_tick; /* time tick of vertex pos */
    double m_vertex_pos_wire_p0;
    double m_vertex_pos_wire_p2;
    double m_vertex_pos_wire_p1;
    int m_reco_vertex_in_SCB; /* is vertex in SCB: 0- No, 1- Yes */
    double m_reco_vertex_dist_to_SCB; /* dist between vertex to SCB */
    double m_reco_vertex_dist_to_active_TPC; /* dist from vertex to closest active TPC wall, -999 if not in active TPC */
    double m_reco_vertex_dist_to_CPA;


    int m_reco_asso_showers;

    //      double m_reco_vertex_to_nearest_dead_wire_plane0;
    //      double m_reco_vertex_to_nearest_dead_wire_plane1;
    //      double m_reco_vertex_to_nearest_dead_wire_plane2;

    //added eventweight
    //-------------- EventWeight related variables -------------

    int m_run_number_eventweight;
    int m_subrun_number_eventweight;
    int m_event_number_eventweight;

    double m_mcflux_nu_pos_x;
    double m_mcflux_nu_pos_y;
    double m_mcflux_nu_pos_z;
    double m_mcflux_nu_mom_x;
    double m_mcflux_nu_mom_y;
    double m_mcflux_nu_mom_z;
    double m_mcflux_nu_mom_E;
    int m_mcflux_ntype;
    int m_mcflux_ptype;
    double m_mcflux_nimpwt;
    double m_mcflux_dk2gen;
    double m_mcflux_nenergyn;
    double m_mcflux_tpx;
    double m_mcflux_tpy;
    double m_mcflux_tpz;
    double m_mcflux_vx;
    double m_mcflux_vy;
    double m_mcflux_vz;
    int m_mcflux_tptype;
    int m_mctruth_nparticles;
    int m_mctruth_particles_track_Id[100];
    int m_mctruth_particles_pdg_code[100];
    int m_mctruth_particles_mother[100];
    int m_mctruth_particles_status_code[100];
    int m_mctruth_particles_num_daughters[100]; //other similar variables
    int m_mctruth_particles_daughters[100][100];
    double m_mctruth_particles_Gvx[100];
    double m_mctruth_particles_Gvy[100];
    double m_mctruth_particles_Gvz[100];
    double m_mctruth_particles_Gvt[100];
    double m_mctruth_particles_px0[100];
    double m_mctruth_particles_py0[100];
    double m_mctruth_particles_pz0[100];
    double m_mctruth_particles_e0[100];
    int m_mctruth_particles_rescatter[100];
    double m_mctruth_particles_polx[100];
    double m_mctruth_particles_poly[100];
    double m_mctruth_particles_polz[100];
    int m_mctruth_neutrino_ccnc;
    int m_mctruth_neutrino_mode;
    int m_mctruth_neutrino_interaction_type;
    int m_mctruth_neutrino_target;
    int m_mctruth_neutrino_nucleon;
    int m_mctruth_neutrino_quark;
    double m_mctruth_neutrino_w;
    double m_mctruth_neutrino_x;
    double m_mctruth_neutrino_y;
    double m_mctruth_neutrino_qsqr;
    bool m_gtruth_is_sea_quark;
    int m_gtruth_tgt_pdg;
    int m_gtruth_tgt_Z;
    int m_gtruth_tgt_A;
    double m_gtruth_tgt_p4_x;
    double m_gtruth_tgt_p4_y;
    double m_gtruth_tgt_p4_z;
    double m_gtruth_tgt_p4_E;
    double m_gtruth_weight;
    double m_gtruth_probability;
    double m_gtruth_xsec;
    double m_gtruth_diff_xsec;
    int m_gtruth_gphase_space;
    double m_gtruth_vertex_x;
    double m_gtruth_vertex_y;
    double m_gtruth_vertex_z;
    double m_gtruth_vertex_T;
    int m_gtruth_gscatter;
    int m_gtruth_gint;
    int m_gtruth_res_num;
    int m_gtruth_num_piplus;
    int m_gtruth_num_pi0;
    int m_gtruth_num_piminus;
    int m_gtruth_num_proton;
    int m_gtruth_num_neutron;
    bool m_gtruth_is_charm;
    bool m_gtruth_is_strange;
    int m_gtruth_charm_hadron_pdg;
    int m_gtruth_strange_hadron_pdg;
    int m_gtruth_decay_mode;
    double m_gtruth_gx;
    double m_gtruth_gy;
    double m_gtruth_gt;
    double m_gtruth_gw;
    double m_gtruth_gQ2;
    double m_gtruth_gq2;
    int m_gtruth_probe_pdg;
    double m_gtruth_probe_p4_x;
    double m_gtruth_probe_p4_y;
    double m_gtruth_probe_p4_z;
    double m_gtruth_probe_p4_E;
    double m_gtruth_hit_nuc_p4_x;
    double m_gtruth_hit_nuc_p4_y;
    double m_gtruth_hit_nuc_p4_z;
    double m_gtruth_hit_nuc_p4_E;
    double m_gtruth_hit_nuc_pos;
    double m_gtruth_fs_had_syst_p4_x;
    double m_gtruth_fs_had_syst_p4_y;
    double m_gtruth_fs_had_syst_p4_z;
    double m_gtruth_fs_had_syst_p4_E;


    //----------- CRT related variables -----------------
    //for crt hits from the CRT veto product
    int m_CRT_veto_nhits;  /* number of CRT veto hits */
    std::vector<double> m_CRT_veto_hit_PE;  

    //fields storing information about the CRT hit closest to the flash
    double m_CRT_min_hit_time;
    double m_CRT_min_hit_PE;
    double m_CRT_min_hit_x;
    double m_CRT_min_hit_y;
    double m_CRT_min_hit_z;

    //Fields storing information about all CRT hits in event
    std::vector<double> m_CRT_hits_time;
    std::vector<double> m_CRT_hits_PE;
    std::vector<double> m_CRT_hits_x;
    std::vector<double> m_CRT_hits_y;
    std::vector<double> m_CRT_hits_z;
    double m_CRT_dt; //time between flash and nearest CRT hit

    double m_genie_spline_weight;
    double m_genie_CV_tune_weight;

    double m_photonu_weight_low;
    double m_photonu_weight_high;

    //------------ Track related Variables -------------

    int m_reco_asso_tracks;  /* number of track. (temp: will figure out what associate means later) */
    std::vector<int>    m_reco_track_num_daughters;
    std::vector<double> m_reco_track_daughter_trackscore;  /* track score of this reco track's first daughter */
    std::vector<double> m_reco_track_length;  /* whole length of the reco track */
    std::vector<double> m_reco_track_dirx;    /* need to understand what the pair track->Direction() returns means*/
    std::vector<double> m_reco_track_diry;
    std::vector<double> m_reco_track_dirz;
    std::vector<double> m_reco_track_startx;  /* start pos of the track in cartesian X */
    std::vector<double> m_reco_track_starty;
    std::vector<double> m_reco_track_startz;
    std::vector<double> m_reco_track_endx;    /* end of the track in cartesian X */
    std::vector<double> m_reco_track_endy;
    std::vector<double> m_reco_track_endz;
    std::vector<double> m_reco_track_end_dist_to_active_TPC;  /* min distance from track end to TPC active volume boundaries */
    std::vector<double> m_reco_track_start_dist_to_active_TPC; /* min dist from trk start to TPC active boundaries */
    std::vector<double> m_reco_track_end_dist_to_CPA;
    std::vector<double> m_reco_track_start_dist_to_CPA;
    std::vector<double> m_reco_track_end_dist_to_SCB;          /* min dist from track end to SCB */
    std::vector<double> m_reco_track_start_dist_to_SCB;
    std::vector<int>    m_reco_track_end_in_SCB;   /* if track end is in SCB boundary, 1- yes, 0- no */
    std::vector<int>    m_reco_track_start_in_SCB;
    std::vector<double> m_reco_track_calo_energy_plane0;  /* energy sum of hits on plane 0 that correspond to the reco track */
    std::vector<double> m_reco_track_calo_energy_plane1;
    std::vector<double> m_reco_track_calo_energy_plane2;
    std::vector<double> m_reco_track_calo_energy_max;    /* max energy of 3 plane for the reco track */

    std::vector<double> m_reco_track_theta_yz; /* theta, phi of the track */
    std::vector<double> m_reco_track_phi_yx;

    std::vector<int>    m_reco_track_num_trajpoints;  /* number of valid points in the track */
    std::vector<int>    m_reco_track_num_spacepoints;  /* number of recob::spacepoints coresponding to the reco track */
    std::vector<double> m_reco_track_proton_kinetic_energy; /* energy of the track, under the asssumption it's a proton track 
                                 * set to -9999 if m_run_pi0_filter is set to true */

    std::vector<size_t> m_reco_track_ordered_energy_index; /* index of m_reco_track_proton_kinetic_energy such that element values are in descending order */
    std::vector<size_t> m_reco_track_ordered_displacement_index; /* index of m_reco_track_length so that track length are in descending order */
    std::vector<double> m_reco_track_spacepoint_principal0; /* PCA of reco track (in 3D spacepoint) */
    std::vector<double> m_reco_track_spacepoint_principal1;
    std::vector<double> m_reco_track_spacepoint_principal2;

    std::vector<double> m_reco_track_spacepoint_chi;  /* essentially sum of square of distances between spacepoint and the track line*/
    std::vector<double> m_reco_track_spacepoint_max_dist; /* max distance between a track and its coresponding spacepoints */


    //corresponding variables on the best plane of reco track, which is defined as such------
    //plane 2 have good hits, then plane 2 is the best-plane
    // which plane of plane 0 and 1 has more good hits will be best plane
    //one of 3 planes has good hits, then best-plane is set to -1
    std::vector<std::vector<double>> m_reco_track_trunc_dEdx_best_plane;
    std::vector<std::vector<double>> m_reco_track_resrange_best_plane;
    std::vector<std::vector<double>> m_reco_track_dEdx_best_plane;
    std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p0;
    std::vector<std::vector<double>> m_reco_track_resrange_p0; /* vec of residual range of good hits per reco track */
    std::vector<std::vector<double>> m_reco_track_dEdx_p0;     /* vec of dEdx of good hits per reco track */
    std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p1;
    std::vector<std::vector<double>> m_reco_track_resrange_p1;
    std::vector<std::vector<double>> m_reco_track_dEdx_p1;
    std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p2;


    std::vector<int>    m_reco_track_best_calo_plane;
    std::vector<double> m_reco_track_mean_dEdx_best_plane;
    std::vector<double> m_reco_track_mean_dEdx_start_half_best_plane;
    std::vector<double> m_reco_track_mean_dEdx_end_half_best_plane;
    std::vector<int>    m_reco_track_good_calo_best_plane;
    std::vector<double> m_reco_track_mean_trunc_dEdx_best_plane;
    std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_best_plane;
    std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_best_plane;
    std::vector<double> m_reco_track_trunc_PIDA_best_plane;


    std::vector<double> m_reco_track_mean_dEdx_p0;  /* mean dEdx of hits on plane 0 of the reco track */
    std::vector<double> m_reco_track_mean_dEdx_start_half_p0; /* mean dEdx of first half of the track */
    std::vector<double> m_reco_track_mean_dEdx_end_half_p0;
    std::vector<int>    m_reco_track_good_calo_p0; /* number of good dEdx hits on plane 0 of track calorimetry */
    std::vector<double> m_reco_track_mean_trunc_dEdx_p0;  /* mean of truncated dEdx's of good hits */
    std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p0; /*mean of first half of trucated dEdx's of good hits */
    std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p0;
    std::vector<double> m_reco_track_trunc_PIDA_p0; /* mean of constant A in residual range formula, calc'd from good hits */

    std::vector<double> m_reco_track_mean_dEdx_p1;
    std::vector<double> m_reco_track_mean_dEdx_start_half_p1;
    std::vector<double> m_reco_track_mean_dEdx_end_half_p1;
    std::vector<int>    m_reco_track_good_calo_p1;
    std::vector<double> m_reco_track_mean_trunc_dEdx_p1;
    std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p1;
    std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p1;
    std::vector<double> m_reco_track_trunc_PIDA_p1;

    std::vector<double> m_reco_track_mean_dEdx_p2;
    std::vector<double> m_reco_track_mean_dEdx_start_half_p2;
    std::vector<double> m_reco_track_mean_dEdx_end_half_p2;
    std::vector<int>    m_reco_track_good_calo_p2;
    std::vector<double> m_reco_track_mean_trunc_dEdx_p2;
    std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p2;
    std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p2;
    std::vector<double> m_reco_track_trunc_PIDA_p2;
    std::vector<std::vector<double>> m_reco_track_resrange_p2;
    std::vector<std::vector<double>> m_reco_track_dEdx_p2;

    std::vector<int> m_reco_track_num_calo_hits_p0; /* number of hits in calorimetry on plane 0 of each reco track */
    std::vector<int> m_reco_track_num_calo_hits_p1;
    std::vector<int> m_reco_track_num_calo_hits_p2;

    //    vector<double> m_reco_track_end_to_nearest_dead_wire_plane0; /* distance between track end and the nearest dead wire on plane*/
    //    vector<double> m_reco_track_end_to_nearest_dead_wire_plane1;
    //    vector<double> m_reco_track_end_to_nearest_dead_wire_plane2;

    std::vector<int>    m_reco_track_sliceId; //the slice id for the slice continaing the reco track
    std::vector<double> m_reco_track_nuscore; //the neutrino score of the slice containing the reco track
    std::vector<bool>   m_reco_track_isclearcosmic;//true if reco track is in a clear cosmic slice
    std::vector<double> m_reco_track_trackscore; /* track score of reco track, -999 if track is not found in PFPToTrackScoreMap */
    std::vector<int>    m_reco_track_pfparticle_pdg; /* PDG of track's corresponding PFParticle, -999 if track is not found in PFPToTrackScoreMap*/
    std::vector<bool>   m_reco_track_is_nuslice;  /* if reco track is in a neutrino slice */




    std::vector<int> m_sim_track_matched;  /* if reco track has been matched to a MCParticle, 1-YES, 0-NO */

    //--- energy, mass, pdg ..etc.. of the matched MCParticle of reco track -----
    std::vector<double> m_sim_track_overlay_fraction;
    std::vector<double> m_sim_track_energy;
    std::vector<double> m_sim_track_mass;
    std::vector<double> m_sim_track_kinetic_energy;
    std::vector<int>    m_sim_track_pdg;
    std::vector<int>    m_sim_track_parent_pdg;

    std::vector<std::string> m_sim_track_process;
    std::vector<int>    m_sim_track_origin;   /* truth origin of the matched MCParticle */
    std::vector<double> m_sim_track_startx;  /* space-charge corrected start point of the match MCParticle */
    std::vector<double> m_sim_track_starty;
    std::vector<double> m_sim_track_startz;
    std::vector<double> m_sim_track_px;
    std::vector<double> m_sim_track_py;
    std::vector<double> m_sim_track_pz;
    std::vector<double> m_sim_track_endx;  /* space-charge corrected end-point of the matched MCParticle */
    std::vector<double> m_sim_track_endy;
    std::vector<double> m_sim_track_endz;
    std::vector<double> m_sim_track_length; /* track length calculated based on the SC-corrected start and end point of the matched MCParticle */

    std::vector<int>    m_sim_track_trackID;
    //--- energy, mass, pdg ..etc.. of the matched MCParticle of reco track -----



    std::vector<int>    m_sim_track_sliceId; //the slice id for the slice continaing the sim track, based on corresponding recob:PFP
    std::vector<double> m_sim_track_nuscore; //the neutrino score of the slice containing the sim track
    std::vector<bool>   m_sim_track_isclearcosmic;//true if sim track is in a clear cosmic slice

    std::vector<double> m_isolation_min_dist_trk_shr; /* minimum distance betwee shower hits and track hits on each plane 
                               ere is no shower hits, set to 999
                               ere is shower hits but no track hits, set to -999
                               */

    std::vector<double> m_isolation_nearest_shr_hit_to_trk_wire; /* the wire number of shower hit closest to track hits */  
    std::vector<double> m_isolation_nearest_shr_hit_to_trk_time; /* the time tick of shower hit closest to track hits in the slice */


    std::vector<double> m_isolation_num_shr_hits_win_1cm_trk; /* number of shower hits whose min distance to track hits <= 1cm 
                                   ch plane (this is a 3 element vector) */      
    std::vector<double> m_isolation_num_shr_hits_win_2cm_trk;      
    std::vector<double> m_isolation_num_shr_hits_win_5cm_trk;      
    std::vector<double> m_isolation_num_shr_hits_win_10cm_trk;      


    std::vector<double> m_isolation_min_dist_trk_unassoc; /* of all unassociated hits, min distance to closest track hits 
                                 o -999 if there is no unassociated hits or track hits on plane
                                 */

    std::vector<double> m_isolation_nearest_unassoc_hit_to_trk_wire;/* wire number of the unassociated hit that of all is nearest to track hits in the slice */  
    std::vector<double> m_isolation_nearest_unassoc_hit_to_trk_time; /* time tick of the unasso hit that is nearest to track hits in the slice */
    std::vector<double> m_isolation_num_unassoc_hits_win_1cm_trk; /* number of unasso hits whose min distance to track hits <= 1cm
                                     ch plane (this vector has 3 elements) */ 
    std::vector<double> m_isolation_num_unassoc_hits_win_2cm_trk;      
    std::vector<double> m_isolation_num_unassoc_hits_win_5cm_trk;      
    std::vector<double> m_isolation_num_unassoc_hits_win_10cm_trk;      
    std::vector<int>    m_reco_shower_num_daughters;
    std::vector<double> m_reco_shower_daughter_trackscore;

    std::vector<int>    m_reco_shower3d_exists;
    std::vector<double> m_reco_shower3d_startx;
    std::vector<double> m_reco_shower3d_starty;
    std::vector<double> m_reco_shower3d_startz;
    std::vector<double> m_reco_shower3d_dirx;
    std::vector<double> m_reco_shower3d_diry;
    std::vector<double> m_reco_shower3d_dirz;
    std::vector<double> m_reco_shower3d_theta_yz; /* theta, phi of the 3D shower (direction) */
    std::vector<double> m_reco_shower3d_phi_yx;

    std::vector<double> m_reco_shower3d_openingangle;
    std::vector<double> m_reco_shower3d_length;
    std::vector<double> m_reco_shower3d_conversion_distance;

    std::vector<double> m_reco_shower3d_impact_parameter;  /* distance between vertex and 3D shower direction */
    std::vector<double> m_reco_shower3d_implied_dirx; /* X component of the unit vector point from vertex to 3D shower start */
    std::vector<double> m_reco_shower3d_implied_diry;
    std::vector<double> m_reco_shower3d_implied_dirz;

    std::vector<double> m_reco_shower3d_energy_plane0;
    std::vector<double> m_reco_shower3d_energy_plane1;
    std::vector<double> m_reco_shower3d_energy_plane2;

    std::vector<double> m_reco_shower3d_dEdx_plane0;
    std::vector<double> m_reco_shower3d_dEdx_plane1;
    std::vector<double> m_reco_shower3d_dEdx_plane2;



    std::vector<double> m_reco_shower_startx;
    std::vector<double> m_reco_shower_starty;
    std::vector<double> m_reco_shower_startz;
    std::vector<double> m_reco_shower_start_dist_to_active_TPC; /* distance from shower start to closest TPC wall */
    std::vector<double> m_reco_shower_start_dist_to_CPA; /* distance from shower start to closest TPC wall */
    std::vector<double> m_reco_shower_start_dist_to_SCB;
    std::vector<int>    m_reco_shower_start_in_SCB;
    std::vector<double> m_reco_shower_end_dist_to_active_TPC;
    std::vector<double> m_reco_shower_end_dist_to_SCB;

    std::vector<double> m_reco_shower_dirx; /* X component of shower direction */
    std::vector<double> m_reco_shower_diry;
    std::vector<double> m_reco_shower_dirz;
    std::vector<double> m_reco_shower_theta_yz; /* theta, phi of the shower direction */
    std::vector<double> m_reco_shower_phi_yx;

    std::vector<double> m_reco_shower_openingangle;
    std::vector<double> m_reco_shower_length;
    std::vector<double> m_reco_shower_conversion_distance;  /* distance between shower start and vertex */

    std::vector<double> m_reco_shower_impact_parameter; /* distance from vertex to the shower direction */
    std::vector<double> m_reco_shower_implied_dirx; /* the X component of the unit vector pointing from vertex to shower start */
    std::vector<double> m_reco_shower_implied_diry;
    std::vector<double> m_reco_shower_implied_dirz;

    std::vector<int> m_reco_shower_delaunay_num_triangles_plane0; /* num of delaunay triangles found on plane 0 for each shower */
    std::vector<int> m_reco_shower_delaunay_num_triangles_plane1;
    std::vector<int> m_reco_shower_delaunay_num_triangles_plane2;

    //    vector<double> m_reco_shower_start_to_nearest_dead_wire_plane0;/* dist from shower start to nearest dead wire on plane 0 */
    //    vector<double> m_reco_shower_start_to_nearest_dead_wire_plane1;
    //    vector<double> m_reco_shower_start_to_nearest_dead_wire_plane2;



    std::vector<double> m_reco_shower_flash_shortest_distz;
    std::vector<double> m_reco_shower_flash_shortest_disty;
    std::vector<double> m_reco_shower_flash_shortest_distyz;

    std::vector<int> m_reco_shower_flash_shortest_index_z;
    std::vector<int> m_reco_shower_flash_shortest_index_y;
    std::vector<int> m_reco_shower_flash_shortest_index_yz;

    double  m_flash_optfltr_pe_beam;
    double  m_flash_optfltr_pe_beam_tot;
    double  m_flash_optfltr_pe_veto;
    double  m_flash_optfltr_pe_veto_tot;

    std::vector<int> m_reco_shower_num_hits_plane0; /* number of hits on plane 0 for each shower */
    std::vector<int> m_reco_shower_num_hits_plane1;
    std::vector<int> m_reco_shower_num_hits_plane2;
    std::vector<double> m_reco_shower_delaunay_area_plane0; /* total area of delaunay triangles found on plane 0 for each shower */
    std::vector<double> m_reco_shower_delaunay_area_plane1;
    std::vector<double> m_reco_shower_delaunay_area_plane2;
    std::vector<int>    m_reco_shower_sliceId; //the slice id for the slice continaing the reco shower
    std::vector<double> m_reco_shower_nuscore; //the neutrino score of the slice containing the reco shower
    std::vector<bool>   m_reco_shower_isclearcosmic;//true if reco shower is in a clear cosmic slice
    std::vector<bool>   m_reco_shower_is_nuslice;//true if reco shower is in a clear cosmic slice
    std::vector<double> m_reco_shower_trackscore;
    std::vector<double> m_reco_shower_pfparticle_pdg;
    std::vector<double> m_reco_shower_kalman_exists; /* if there is a kalman track and reco::Calo related to this shower - 0, 1 */
    std::vector<double> m_reco_shower_kalman_median_dEdx_plane0;
    std::vector<double> m_reco_shower_kalman_median_dEdx_plane1;
    std::vector<double> m_reco_shower_kalman_median_dEdx_plane2;
    std::vector<double> m_reco_shower_kalman_median_dEdx_allplane;
    std::vector<double> m_reco_shower_kalman_mean_dEdx_plane0;
    std::vector<double> m_reco_shower_kalman_mean_dEdx_plane1;
    std::vector<double> m_reco_shower_kalman_mean_dEdx_plane2;




    std::vector<int> m_sim_shower_matched;  /* whether shower has been matched to a MCParticle, 0 - False, 1 - True */

    std::vector<double> m_sim_shower_energy;
    std::vector<double> m_sim_shower_kinetic_energy;
    std::vector<double> m_sim_shower_mass;
    std::vector<int> m_sim_shower_pdg;
    std::vector<int> m_sim_shower_trackID;
    std::vector<int> m_sim_shower_parent_pdg;
    std::vector<int> m_sim_shower_parent_trackID;
    std::vector<int> m_sim_shower_origin;
    std::vector<std::string> m_sim_shower_process;
    std::vector<std::string> m_sim_shower_end_process;


    std::vector<double> m_sim_shower_start_x;  /* space charge corrected shower starting point */
    std::vector<double> m_sim_shower_start_y;
    std::vector<double> m_sim_shower_start_z;
    std::vector<double> m_sim_shower_vertex_x;  /* spacecharge corrected shower vertex */
    std::vector<double> m_sim_shower_vertex_y;
    std::vector<double> m_sim_shower_vertex_z;

    std::vector<double> m_sim_shower_px;
    std::vector<double> m_sim_shower_py;
    std::vector<double> m_sim_shower_pz;
    std::vector<int> m_sim_shower_is_true_shower;
    std::vector<int> m_sim_shower_best_matched_plane;
    std::vector<double> m_sim_shower_matched_energy_fraction_plane0; /* fraction of energy of the best-matched mother for shower on 
                                      0 over all energy deposited on plane 0 by the shower */
    std::vector<double> m_sim_shower_matched_energy_fraction_plane1;
    std::vector<double> m_sim_shower_matched_energy_fraction_plane2;
    std::vector<double> m_sim_shower_overlay_fraction; /* fraction of hits from overlay over all hits in the shower */
    std::vector<int> m_sim_shower_sliceId; //the slice id for the slice continaing the sim shower matched to reco
    std::vector<double> m_sim_shower_nuscore; //the neutrino score of the slice containing the sim shower matched to reco
    std::vector<bool> m_sim_shower_isclearcosmic;//true if sim shower matched to reco is in a clear cosmic slice
    std::vector<bool> m_sim_shower_is_nuslice;//true if sim shower matched to reco is in a clear cosmic slice



    int m_mctruth_num;
    int m_mctruth_origin;
    double m_mctruth_nu_E;
    double m_mctruth_nu_vertex_x;
    double m_mctruth_nu_vertex_y;
    double m_mctruth_nu_vertex_z;
    double m_mctruth_nu_vertex_t;
    double m_mctruth_reco_vertex_dist;
    double m_mctruth_lepton_E;
    int m_mctruth_nu_pdg;
    int m_mctruth_lepton_pdg;
    int m_mctruth_mode ;
    int m_mctruth_interaction_type ;
    int m_mctruth_ccnc;
    double m_mctruth_qsqr;
    int m_mctruth_num_daughter_particles;
    std::vector<int> m_mctruth_daughters_pdg;
    std::vector<double> m_mctruth_daughters_E;
    std::vector<int> m_mctruth_daughters_status_code;
    std::vector<int> m_mctruth_daughters_trackID;
    std::vector<int> m_mctruth_daughters_mother_trackID;
    std::vector<double> m_mctruth_daughters_px;
    std::vector<double> m_mctruth_daughters_py;
    std::vector<double> m_mctruth_daughters_pz;
    std::vector<double> m_mctruth_daughters_startx;
    std::vector<double> m_mctruth_daughters_starty;
    std::vector<double> m_mctruth_daughters_startz;
    std::vector<double> m_mctruth_daughters_time;
    std::vector<double> m_mctruth_daughters_endx;
    std::vector<double> m_mctruth_daughters_endy;
    std::vector<double> m_mctruth_daughters_endz;
    std::vector<double> m_mctruth_daughters_endtime;
    std::vector<std::string> m_mctruth_daughters_process;
    std::vector<std::string> m_mctruth_daughters_end_process;
    int   m_mctruth_num_exiting_photons ;
    int   m_mctruth_num_exiting_protons ;
    int   m_mctruth_num_exiting_pi0 ;
    int   m_mctruth_num_exiting_pipm ;
    int   m_mctruth_num_exiting_neutrons; 
    int   m_mctruth_num_exiting_delta0; 
    int   m_mctruth_num_exiting_deltapm; 
    int   m_mctruth_num_exiting_deltapp; 
    double m_mctruth_leading_exiting_proton_energy;
    int m_mctruth_is_delta_radiative;
    int m_mctruth_delta_radiative_1g1p_or_1g1n;
    double m_mctruth_delta_photon_energy;
    double m_mctruth_delta_proton_energy;
    double m_mctruth_delta_neutron_energy;
    std::vector<int> m_mctruth_exiting_delta0_num_daughters;

    std::vector<int> m_mctruth_exiting_photon_trackID;
    std::vector<int> m_mctruth_exiting_photon_mother_trackID;
    std::vector<int> m_mctruth_exiting_photon_from_delta_decay;
    std::vector<double> m_mctruth_exiting_photon_energy;
    std::vector<double> m_mctruth_exiting_photon_px;
    std::vector<double> m_mctruth_exiting_photon_py;
    std::vector<double> m_mctruth_exiting_photon_pz;

    std::vector<int> m_mctruth_exiting_proton_trackID;
    std::vector<int> m_mctruth_exiting_proton_mother_trackID;
    std::vector<int> m_mctruth_exiting_proton_from_delta_decay;
    std::vector<double> m_mctruth_exiting_proton_energy;
    std::vector<double> m_mctruth_exiting_proton_px;
    std::vector<double> m_mctruth_exiting_proton_py;
    std::vector<double> m_mctruth_exiting_proton_pz;

    std::vector<int> m_mctruth_exiting_neutron_trackID;
    std::vector<int> m_mctruth_exiting_neutron_mother_trackID;
    std::vector<int> m_mctruth_exiting_neutron_from_delta_decay;
    std::vector<double> m_mctruth_exiting_neutron_energy;
    std::vector<double> m_mctruth_exiting_neutron_px;
    std::vector<double> m_mctruth_exiting_neutron_py;
    std::vector<double> m_mctruth_exiting_neutron_pz;
    int  m_mctruth_num_reconstructable_protons;
    bool  m_mctruth_is_reconstructable_1g1p;
    bool  m_mctruth_is_reconstructable_1g0p;

    std::vector<double> m_mctruth_exiting_pi0_E;
    std::vector<double> m_mctruth_exiting_pi0_mom;
    std::vector<double> m_mctruth_exiting_pi0_px;
    std::vector<double> m_mctruth_exiting_pi0_py;
    std::vector<double> m_mctruth_exiting_pi0_pz;

    double              m_mctruth_pi0_leading_photon_energy;
    std::string         m_mctruth_pi0_leading_photon_end_process;
    double              m_mctruth_pi0_subleading_photon_energy;
    std::string         m_mctruth_pi0_subleading_photon_end_process;
    int                 m_mctruth_pi0_leading_photon_exiting_TPC;
    int                 m_mctruth_pi0_subleading_photon_exiting_TPC;
    std::vector<double> m_mctruth_pi0_subleading_photon_end;
    std::vector<double> m_mctruth_pi0_subleading_photon_start;
    std::vector<double> m_mctruth_pi0_leading_photon_end;
    std::vector<double> m_mctruth_pi0_leading_photon_start;
    std::vector<double> m_mctruth_pi0_leading_photon_mom;
    std::vector<double> m_mctruth_pi0_subleading_photon_mom;


    //the calo calculated quantities 
    std::vector<double> m_reco_shower_energy_max; //for each hit in a shower, converts Q->E, and sums. The max energy of all planes
    std::vector<double> m_reco_shower_energy_plane0; /* shower energy (summed hit energy) on plan 0 */
    std::vector<double> m_reco_shower_energy_plane1;
    std::vector<double> m_reco_shower_energy_plane2;
    std::vector<double> m_reco_shower_reclustered_energy_max;
    std::vector<double> m_reco_shower_reclustered_energy_plane0; /* total energy of the reco shower, and unassociated hit clusters 
                                    enough to it */
    std::vector<double> m_reco_shower_reclustered_energy_plane1;
    std::vector<double> m_reco_shower_reclustered_energy_plane2;
    std::vector<double> m_reco_shower_plane0;
    std::vector<double> m_reco_shower_plane1;
    std::vector<double> m_reco_shower_plane2;
    std::vector<double> m_reco_shower_plane0_nhits; /* num of shower hits on plane 0 */
    std::vector<double> m_reco_shower_plane1_nhits;
    std::vector<double> m_reco_shower_plane2_nhits;
    std::vector<double> m_reco_shower_plane0_meanRMS; /* the mean of RMS of the shower hit shape (in tick unit) on plane 0 */
    std::vector<double> m_reco_shower_plane1_meanRMS;
    std::vector<double> m_reco_shower_plane2_meanRMS;
    std::vector<int> m_reco_shower_hit_wire;
    std::vector<int> m_reco_shower_hit_plane;
    std::vector<double> m_reco_shower_hit_tick;
    std::vector<double> m_reco_shower_spacepoint_x;
    std::vector<double> m_reco_shower_spacepoint_z;
    std::vector<double> m_reco_shower_spacepoint_y;
    std::vector<size_t>  m_reco_shower_ordered_energy_index; /* indices of 'm_reco_shower_energy_max' such that energy max is in descending order */
    std::vector<std::vector<double>> m_reco_shower_dQdx_plane0; //for each shower, looks at the hits for all clusters in the plane, stores the dQ/dx for each hit 
    std::vector<std::vector<double>> m_reco_shower_dQdx_plane1;
    std::vector<std::vector<double>> m_reco_shower_dQdx_plane2;
    std::vector<std::vector<double>> m_reco_shower_dEdx_plane0; //dE/dx from the calculated dQ/dx for each hit of all clusters in shower on plane   
    std::vector<std::vector<double>> m_reco_shower_dEdx_plane1;
    std::vector<std::vector<double>> m_reco_shower_dEdx_plane2;
    std::vector<double> m_reco_shower_dEdx_plane0_mean; /* mean of dE/dx of each hit in shower */
    std::vector<double> m_reco_shower_dEdx_plane1_mean;
    std::vector<double> m_reco_shower_dEdx_plane2_mean;
    std::vector<double> m_reco_shower_dEdx_plane0_max;
    std::vector<double> m_reco_shower_dEdx_plane1_max;
    std::vector<double> m_reco_shower_dEdx_plane2_max;
    std::vector<double> m_reco_shower_dEdx_plane0_min;
    std::vector<double> m_reco_shower_dEdx_plane1_min;
    std::vector<double> m_reco_shower_dEdx_plane2_min;
    std::vector<double> m_reco_shower_dEdx_plane0_median;/* median of dE/dx of each hit in shower (median of vector element of m_reco_shower_dEdx_plane0) */
    std::vector<double> m_reco_shower_dEdx_plane1_median;
    std::vector<double> m_reco_shower_dEdx_plane2_median;
    std::vector<double>  m_reco_shower_angle_wrt_wires_plane0; /* angle between shower direction and wire dir on plane, in radian*/
    std::vector<double>  m_reco_shower_angle_wrt_wires_plane1;
    std::vector<double>  m_reco_shower_angle_wrt_wires_plane2;
    std::vector<double>  m_reco_shower_dEdx_amalgamated;
    std::vector<int>  m_reco_shower_dEdx_amalgamated_nhits;
    std::vector<double> m_reco_shower_dQdx_plane0_median;/* median of dQ/dx of each hit in shower (median of m_reco_shower_dQdx_plane0) */
    std::vector<double> m_reco_shower_dQdx_plane1_median;
    std::vector<double> m_reco_shower_dQdx_plane2_median;
    std::vector<double> m_reco_shower_dEdx_plane0_nhits; /* number of hits of all clusters of the shower on plane 0 */
    std::vector<double> m_reco_shower_dEdx_plane1_nhits;
    std::vector<double> m_reco_shower_dEdx_plane2_nhits;

    //      related variables
    std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane0;
    std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane1;
    std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane2;
    std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane0;
    std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane1;
    std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane2;
    std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane0;
    std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane1;
    std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane2;
    std::vector<double> m_reco_track_pid_pida_plane0;
    std::vector<double> m_reco_track_pid_pida_plane1;
    std::vector<double> m_reco_track_pid_pida_plane2;
    std::vector<double> m_reco_track_pid_chi2_mu_plane0;
    std::vector<double> m_reco_track_pid_chi2_mu_plane1;
    std::vector<double> m_reco_track_pid_chi2_mu_plane2;
    std::vector<double> m_reco_track_pid_chi2_p_plane0;
    std::vector<double> m_reco_track_pid_chi2_p_plane1;
    std::vector<double> m_reco_track_pid_chi2_p_plane2;
    std::vector<double> m_reco_track_pid_three_plane_proton_pid;


    //matching variables
    int  m_reco_slice_num; //total number of slices in the event
    std::vector<double> m_reco_slice_nuscore; //vector of the neutrino score for each slice in an event
    int m_reco_slice_shower_num_matched_signal; //the number of sim showers matched an MCP in the signal def
    int m_reco_slice_track_num_matched_signal; //the number of sim showers matched an MCP in the signal def
    int m_matched_signal_shower_num = 0;  /* number of match showers (that has unique best-matched primary photon ?)  */
    std::vector<int> m_reco_slice_shower_matched_sliceId; //the slice id for each matched shower
    std::vector<int> m_reco_slice_track_matched_sliceId; //the slice id for each matched track

    std::vector<int> m_reco_slice_num_pfps; //the total number of PFP's per slice
    std::vector<int> m_reco_slice_num_showers; //the subset of PFP's that are showers, ie number of showers per slice 
    std::vector<int> m_reco_slice_num_tracks; //the subset of PFP's that are tracks

    std::vector<double> m_reco_slice_shower_matched_energy; //the energy for each matched shower
    std::vector<double> m_reco_slice_track_matched_energy; //the energy for each matched track
    std::vector<double> m_reco_slice_shower_matched_conversion; //the conversion distance for each matched shower
    std::vector<double> m_reco_slice_shower_matched_overlay_frac; //fraction of overlay hits for each matched shower
    //std::map<art::Ptr<recob::PFParticle>, double > & pfParticleToNuScoreMap;//is filled during analyze slices


    //-------  matched shower: reco shower that matches to a primary photon + max energy of 3 plane > 20 + definition being ncdelta----- 
    std::vector<double> m_matched_signal_shower_overlay_fraction;
    //std::vector<double> m_matched_signal_shower_conversion_length;
    std::vector<double> m_matched_signal_shower_true_E;  /* energy of the best-matched MCparticle for the shower */
    std::vector<double> m_matched_signal_shower_nuscore; /* the neutrino score of the slice containing the reco shower */
    std::vector<int> m_matched_signal_shower_sliceId;    /* reco shower slice ID */
    std::vector<bool> m_matched_signal_shower_is_clearcosmic;
    std::vector<bool> m_matched_signal_shower_is_nuslice;
    std::vector<int> m_matched_signal_shower_tracks_in_slice; /* number of showers in the same slice as of this reco shower */
    std::vector<int> m_matched_signal_shower_showers_in_slice; /* number of tracks in the same slice as of this reco shower */


    //-------- for reco tracks that match to a primary proton ---------
    std::vector<double> m_matched_signal_track_true_E; /*  the true energy of matched MCparticle (proton) */
    std::vector<double> m_matched_signal_track_nuscore;  /* nu score of the slice containing the reco track */
    std::vector<int> m_matched_signal_track_sliceId;
    std::vector<bool> m_matched_signal_track_is_clearcosmic; /* if reco track is in clear cosmic slice */
    //  std::vector<bool> m_matched_signal_track_is_nuslice;
    std::vector<bool> m_matched_signal_track_is_nuslice;
    std::vector<int> m_matched_signal_track_tracks_in_slice; /* num of PFP that are tracks in the slice this reco track is in */
    std::vector<int> m_matched_signal_track_showers_in_slice;

    int m_matched_signal_track_num = 0;  /* num of reco tracks matched to primary proton */ 

    //int m_matched_signal_total_num_slices;

    //---------for reco tracks that match to a primary proton ---------

    bool m_reco_1g1p_is_same_slice;
    bool m_reco_1g1p_is_multiple_slices;
    bool m_reco_1g1p_is_nuslice;
    bool m_reco_1g0p_is_nuslice;
    double m_reco_1g1p_nuscore;
    double  m_reco_1g0p_nuscore;
    bool m_is_matched_1g1p;
    bool m_is_matched_1g0p;
    bool m_no_matched_showers;
    bool m_multiple_matched_showers; //if there is more than 1 eligible shower (match to primary photon, pass energy threshold)
    bool m_multiple_matched_tracks; /* if there is more than 1 eligible track (match to primary proton) */
  };


}
#endif // SBNCODE_SINGLEPHOTONANALYSIS_VARIABLES_H
