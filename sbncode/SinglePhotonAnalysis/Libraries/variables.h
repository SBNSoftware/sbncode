#ifndef SBNCODE_SINGLEPHOTONANALYSIS_VARIABLES_H
#define SBNCODE_SINGLEPHOTONANALYSIS_VARIABLES_H

#include "art/Framework/Principal/Handle.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"

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

      extern std::map<int,std::string> is_delta_map;

//DECLARATION HERE
      //Geometry dimensions; Fiducial volume and SCB (no SCB yet?)
      extern std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
      extern std::vector<geo::BoxBoundedGeo> fActiveVolumes;
      extern double m_tpc_active_XMin;
      extern double m_tpc_active_YMin;
      extern double m_tpc_active_ZMin;
      extern double m_tpc_active_XMax;
      extern double m_tpc_active_YMax;
      extern double m_tpc_active_ZMax;

      extern TRandom3 *rangen;
      extern std::string m_shower3dLabel;
      extern std::string m_showerKalmanLabel;
      extern std::string m_showerKalmanCaloLabel;
      extern std::string m_pandoraLabel;         ///< The label for the pandora producer
      extern std::string m_trackLabel;           ///< The label for the track producer from PFParticles

      extern std::string m_showerLabel;          ///< The label for the shower producer from PFParticles
      extern std::string m_caloLabel;            ///< The label for calorimetry associations producer
      extern std::string m_flashLabel;
      extern std::string m_geantModuleLabel;
      extern std::string m_hitfinderLabel;
      extern std::string m_hitMCParticleAssnsLabel;
      extern std::string m_potLabel;
      extern std::string m_generatorLabel;


      extern std::string m_pidLabel;            ///< For PID stuff
      extern std::string m_CRTVetoLabel;
      extern std::string m_CRTTzeroLabel;
      extern std::string m_CRTHitProducer;
      extern std::string m_true_eventweight_label;

      extern bool m_use_PID_algorithms;
      extern bool m_use_delaunay;
      extern int  m_delaunay_max_hits;
      extern bool m_is_verbose;
      extern bool m_print_out_event;
      extern bool m_is_data; // value provided by pset
      extern bool m_is_overlayed;
      extern bool m_is_textgen;
      extern bool m_run_all_pfps;
      extern bool m_has_CRT;
      extern bool m_fill_trees;
      extern bool m_run_pi0_filter; //value provided by pset
      extern bool m_run_pi0_filter_2g1p;
      extern bool m_run_pi0_filter_2g0p;

      extern bool m_runPhotoNuTruth;
      extern bool m_runTrueEventweight;

      extern bool m_runSelectedEvent;  //if it should run only selected events
      extern std::string m_selected_event_list; //full path for the file containing run/subrun/event number of selected events
      extern std::set<std::vector<int>> m_selected_set;  //set of selected events     

      //SEAviwer bits
      extern bool m_runSEAview;
      extern double m_SEAviewPlotDistance;   //parameters related to shower-like object finding
      extern double m_SEAviewHitThreshold;
      extern double  m_SEAviewDbscanMinPts;
      extern double m_SEAviewDbscanEps;
      extern double m_SEAviewMaxPtsLinFit;
      extern bool   m_SEAviewMakePDF;
      extern int m_SEAviewNumRecoShower;
      extern int m_SEAviewNumRecoTrack;

      extern bool m_runSEAviewStub;
      extern double m_SEAviewStubHitThreshold; //parameters related to track-like object finding
      extern double m_SEAviewStubPlotDistance;
      extern double m_SEAviewStubDbscanMinPts;
      extern double m_SEAviewStubDbscanEps;
      extern bool m_SEAviewStubMakePDF;
      extern int m_SEAviewStubNumRecoShower;
      extern int m_SEAviewStubNumRecoTrack;

      extern std::string m_Spline_CV_label; //"eventweight4to4aFix"

      extern bool m_runCRT;
      extern double m_DTOffset;
      extern double  m_Resolution;
      extern std::string  m_DAQHeaderProducer;//"daq"

      //SSS parameters
      extern double m_max_conv_dist;
      extern double m_mass_pi0_mev;

      extern double m_exiting_photon_energy_threshold ;
      extern double m_exiting_proton_energy_threshold ;

      extern geo::GeometryCore const * geom;
      extern double m_work_function;  //value provided by pset
      extern double m_recombination_factor; // value provided by pset

      extern std::vector<double> m_gain_mc; // value provided by pset 
      extern std::vector<double> m_gain_data; 
      extern double m_wire_spacing;

      extern double m_width_dqdx_box; // value provided by pset
      extern double m_length_dqdx_box;


    //------- TTree stuff
      extern TTree* run_subrun_tree;
      extern TTree* pot_tree;
      extern TTree* vertex_tree;
      extern TTree* eventweight_tree;
      extern TTree* ncdelta_slice_tree;

      extern TTree* geant4_tree;

      extern TTree* true_eventweight_tree;
      extern std::map<std::string, std::vector<double>> fmcweight;

       //------------ POT related variables --------------
      extern int m_number_of_events;
      extern int m_number_of_events_in_subrun;
      extern double m_pot_count;
      extern int m_number_of_vertices;

      extern int m_run;
      extern int m_subrun;
      extern double m_subrun_pot;
      extern int m_subrun_counts;

      //------------ Event Related Variables -------------
      extern int m_run_number;
      extern int m_subrun_number;
      extern int m_event_number;
      extern double m_pot_per_event;
      extern double m_pot_per_subrun;

      extern int m_test_matched_hits;
      extern int m_reco_slice_objects;

      //------- Potential Unreconstructed Track Stub related variables ----
      extern int m_trackstub_num_unassociated_hits; /* number of hits in the slice that're associated with neither shower nor tracks */
      extern int m_trackstub_unassociated_hits_below_threshold; /*number of unassociated hits that also didn't pass hit threshold,in the slice*/
      extern int m_trackstub_associated_hits; /* total number of hits from showers and tracks in the slice */


      extern int m_trackstub_num_candidates; /* number of unasso hit clusters which are not close enough to reco showers */
      extern std::vector<int> m_trackstub_candidate_in_nu_slice; /* check if candidate is in neutrino slice: 1->YES, 0->Parts in neutrino slice, -1->Not at all */
      extern std::vector<int> m_trackstub_candidate_num_hits;
      extern std::vector<int> m_trackstub_candidate_num_wires; //number of wires spanned by the candidate cluster
      extern std::vector<int>  m_trackstub_candidate_num_ticks;
      extern std::vector<int>  m_trackstub_candidate_plane; /* on which plan the unasso cluster is */
      extern std::vector<double> m_trackstub_candidate_PCA;
      extern std::vector<double> m_trackstub_candidate_mean_ADC;
      extern std::vector<double> m_trackstub_candidate_ADC_RMS;
      extern std::vector<double> m_trackstub_candidate_veto_score;
      extern std::vector<double> m_trackstub_candidate_mean_tick;
      extern std::vector<double> m_trackstub_candidate_max_tick;
      extern std::vector<double> m_trackstub_candidate_min_tick;
      extern std::vector<double> m_trackstub_candidate_min_wire;
      extern std::vector<double> m_trackstub_candidate_max_wire;
      extern std::vector<double> m_trackstub_candidate_mean_wire;
      extern std::vector<double> m_trackstub_candidate_min_dist;  // min distance from unasso cluter to the vertex */
      extern std::vector<double> m_trackstub_candidate_min_impact_parameter_to_shower; //min impact parameter of all hits in cluster to the recob::shower direction line (on 2D plane)
      extern std::vector<double> m_trackstub_candidate_min_conversion_dist_to_shower_start;  //min distance between hits and recob::shower start (on 2D plane)
      extern std::vector<double> m_trackstub_candidate_min_ioc_to_shower_start;        //min ratio of impact_parameter_to_shower/conversion_dist_to_shower_start of all hits in the cluster
      extern std::vector<double> m_trackstub_candidate_ioc_based_length;    //length of the cluster, calculated based on the IOC of hit
      extern std::vector<double> m_trackstub_candidate_wire_tick_based_length;    //length of the cluster, calculated based on the wire & tick span of the cluster
      extern std::vector<double> m_trackstub_candidate_mean_ADC_first_half;    // mean ADC per hit for the first half of cluster (cluster divided into halves based on hit IOC)
      extern std::vector<double> m_trackstub_candidate_mean_ADC_second_half;
      extern std::vector<double> m_trackstub_candidate_mean_ADC_first_to_second_ratio; // ratio of the mean ADC per hit, first half of cluster over second half.
      extern std::vector<double> m_trackstub_candidate_track_angle_wrt_shower_direction;   //treat cluster as a track, angle between track direction and the shower direction
      extern std::vector<double> m_trackstub_candidate_linear_fit_chi2;    // chi2 from linear fit of the  {wire, tick} distribution of the cluster
      extern std::vector<double> m_trackstub_candidate_energy;
      extern std::vector<int>    m_trackstub_candidate_remerge; // index of the recob::shower candidate cluster is close to (expect it to be -1)
      extern std::vector<int>    m_trackstub_candidate_matched; /* has matched this unasso cluter to a primary MCParticle: 0-No, 1-Yes */
      extern std::vector<double> m_trackstub_candidate_matched_energy_fraction_best_plane; /* matched energy fraction of the best-matched MCParticle on best-plane */ 
      extern std::vector<int>    m_trackstub_candidate_pdg;   /* pdg of the matched MCParticle */
      extern std::vector<int>    m_trackstub_candidate_parent_pdg;
      extern std::vector<int>    m_trackstub_candidate_trackid; /* track ID of the matched MCParticle */
      extern std::vector<double> m_trackstub_candidate_true_energy;  /* true energy of the matched MCParticle */
      extern std::vector<double> m_trackstub_candidate_overlay_fraction; /* fraction of overlay in the unasso cluster hits */

      //------- grouped stub clusters --------------
      extern int m_trackstub_num_candidate_groups;           /* number of groups */ 
      extern std::vector<std::vector<double>> m_grouped_trackstub_candidate_indices; /* indices of stub clusters that are matched as a group */
      extern std::vector<double> m_trackstub_candidate_group_timeoverlap_fraction;   /* minimum fraction of the time overlap of grouped stub clusters */



      //------- Second shower related variables ----
      extern int m_sss_num_unassociated_hits; /* number of hits in the slice that're associated with neither shower nor tracks */
      extern int m_sss_num_unassociated_hits_below_threshold; /*number of unassociated hits that also didn't pass hit threshold,in the slice*/
      extern int m_sss_num_associated_hits; /* total number of hits from showers and tracks in the slice */


      //currently commenting this out for speed as its not used
      //ReadBDT * sssVetov1;

      extern int m_sss_num_candidates; /* number of unasso hit clusters which are not close enough to reco showers */
      extern std::vector<int> m_sss_candidate_in_nu_slice;
      extern std::vector<int> m_sss_candidate_num_hits;
      extern std::vector<int> m_sss_candidate_num_wires; //number of wires spanned by the candidate cluster
      extern std::vector<int>  m_sss_candidate_num_ticks;
      extern std::vector<int>  m_sss_candidate_plane; /* on which plan the unasso cluster is */
      extern std::vector<double> m_sss_candidate_PCA;
      extern std::vector<double> m_sss_candidate_mean_ADC;
      extern std::vector<double> m_sss_candidate_ADC_RMS;
      extern std::vector<double> m_sss_candidate_impact_parameter;
      extern std::vector<double> m_sss_candidate_fit_slope; //slope of the cluster direction
      extern std::vector<double> m_sss_candidate_veto_score;
      extern std::vector<double> m_sss_candidate_fit_constant; //intercept of the cluster direction
      extern std::vector<double> m_sss_candidate_mean_tick;
      extern std::vector<double> m_sss_candidate_max_tick;
      extern std::vector<double> m_sss_candidate_min_tick;
      extern std::vector<double> m_sss_candidate_min_wire;
      extern std::vector<double> m_sss_candidate_max_wire;
      extern std::vector<double> m_sss_candidate_mean_wire;
      extern std::vector<double> m_sss_candidate_min_dist;  // min distance from unasso cluter to the vertex */
      extern std::vector<double> m_sss_candidate_wire_tick_based_length;    //length of the cluster, calculated based on the wire & tick span of the cluster
      extern std::vector<double> m_sss_candidate_energy;
      extern std::vector<double> m_sss_candidate_angle_to_shower;
      extern std::vector<double> m_sss_candidate_closest_neighbour;
      extern std::vector<int>    m_sss_candidate_remerge; // index of the recob::shower candidate cluster is close to (expect it to be -1)
      extern std::vector<int>    m_sss_candidate_matched; /* has matched this unasso cluter to a primary MCParticle: 0-No, 1-Yes */
      extern std::vector<double> m_sss_candidate_matched_energy_fraction_best_plane; /* matched energy fraction of the best-matched MCParticle on best-plane */ 
      extern std::vector<int>    m_sss_candidate_pdg;   /* pdg of the matched MCParticle */
      extern std::vector<int>    m_sss_candidate_parent_pdg;
      extern std::vector<int>    m_sss_candidate_trackid; /* track ID of the matched MCParticle */
      extern std::vector<double> m_sss_candidate_true_energy;
      extern std::vector<double> m_sss_candidate_overlay_fraction; /* fraction of overlay in the unasso cluster hits */



      //------------ sss3d_showers variables are for reco::showers which are in the events, but not in the slice ----

      extern int m_sss3d_num_showers;  /* number of showers in the event but not in the slice */
      extern std::vector<double> m_sss3d_shower_start_x; /* shower start in X axis, for all showers in the event but not in the slice*/
      extern std::vector<double> m_sss3d_shower_start_y;
      extern std::vector<double> m_sss3d_shower_start_z;
      extern std::vector<double> m_sss3d_shower_dir_x; /* shower direction projection on X axis */
      extern std::vector<double> m_sss3d_shower_dir_y;
      extern std::vector<double> m_sss3d_shower_dir_z;
      extern std::vector<double> m_sss3d_shower_length;
      extern std::vector<double> m_sss3d_shower_conversion_dist; /* dist between shower start and vertex*/

      extern std::vector<double> m_sss3d_shower_invariant_mass; /* invariant mass of primary recob::shower, and each shower in the event, 
      extern                           * calculated assuming vertex is where their mother particle decays */

      extern std::vector<double> m_sss3d_shower_implied_invariant_mass; /* similar to invariance mass, except this invariant mass  
      extern                               * is calced direclty using shower direction of two showers */

      extern std::vector<double> m_sss3d_shower_impact_parameter; /* dist between vertex and shower direction line */
      extern std::vector<double> m_sss3d_shower_ioc_ratio; /* ratio of impact parameter over conversion dist 
      extern                          * 0 if the conversion distance is 0*/
      extern std::vector<double> m_sss3d_shower_energy_max; /* max energy of all planes (E summed from hits) */
      extern std::vector<double> m_sss3d_shower_score;
      extern std::vector<int> m_sss3d_slice_nu;
      extern std::vector<int> m_sss3d_slice_clear_cosmic;

      extern bool bool_make_sss_plots;


      //------ max_energy, conversion dist, ioc of the sss3d shower that has the smallest ioc parameter ----
      extern double m_sss3d_ioc_ranked_en;
      extern double m_sss3d_ioc_ranked_conv;
      extern double m_sss3d_ioc_ranked_invar;
      extern double m_sss3d_ioc_ranked_implied_invar;
      extern double m_sss3d_ioc_ranked_ioc;
      extern double m_sss3d_ioc_ranked_opang;
      extern double m_sss3d_ioc_ranked_implied_opang;
      extern int m_sss3d_ioc_ranked_id; //index of the sss3d_shower that has the smallest ioc.

      // --- same parameters, of the sss3d shower whose implied invariant mass together with primary recob::shower is closest to pi0 mass --
      extern double m_sss3d_invar_ranked_en;
      extern double m_sss3d_invar_ranked_conv;
      extern double m_sss3d_invar_ranked_invar;
      extern double m_sss3d_invar_ranked_implied_invar;
      extern double m_sss3d_invar_ranked_ioc;
      extern double m_sss3d_invar_ranked_opang;
      extern double m_sss3d_invar_ranked_implied_opang;
      extern int m_sss3d_invar_ranked_id;


      //--------------- sss2d showers are essentially group of cluters on 3 planes, that have the potential to be a shower -------
      //--------------- they are not recob::showers --------------------------

      // sss2d_ioc_ranked variables are the varaibles (mean shower energy, conv. dist, ioc, etc) of the sss2d shower that has the smallest ioc 
      // sss2d_conv_ranked variables are the varaibles (energy, conv. dist, ioc, etc) of the sss2d shower that has the smallest conv. distance 
      // sss2d_invar_ranked variables are the varaibles of the sss2d shower whose invariant mass together with primary shower is closest to pi0. 
      extern double m_sss2d_ioc_ranked_en;
      extern double m_sss2d_ioc_ranked_conv;
      extern double m_sss2d_ioc_ranked_ioc;
      extern double m_sss2d_ioc_ranked_pca;
      extern double m_sss2d_ioc_ranked_invar;
      extern double m_sss2d_ioc_ranked_angle_to_shower;
      extern int m_sss2d_ioc_ranked_num_planes;

      extern double m_sss2d_conv_ranked_en;
      extern double m_sss2d_conv_ranked_conv;
      extern double m_sss2d_conv_ranked_ioc;
      extern double m_sss2d_conv_ranked_pca;
      extern double m_sss2d_conv_ranked_invar;
      extern double m_sss2d_conv_ranked_angle_to_shower;
      extern int m_sss2d_conv_ranked_num_planes;

      extern double m_sss2d_invar_ranked_en;
      extern double m_sss2d_invar_ranked_conv;
      extern double m_sss2d_invar_ranked_ioc;
      extern double m_sss2d_invar_ranked_pca;
      extern double m_sss2d_invar_ranked_invar;
      extern double m_sss2d_invar_ranked_angle_to_shower;
      extern int m_sss2d_invar_ranked_num_planes;


      //------------ Vertex Related variables -------------
      extern int m_reco_vertex_size;
      extern double m_vertex_pos_x;
      extern double m_vertex_pos_y;
      extern double m_vertex_pos_z;
      extern double m_vertex_pos_tick; /* time tick of vertex pos */
      extern double m_vertex_pos_wire_p0;
      extern double m_vertex_pos_wire_p2;
      extern double m_vertex_pos_wire_p1;
      extern int m_reco_vertex_in_SCB; /* is vertex in SCB: 0- No, 1- Yes */
      extern double m_reco_vertex_dist_to_SCB; /* dist between vertex to SCB */
      extern double m_reco_vertex_dist_to_active_TPC; /* dist from vertex to closest active TPC wall, -999 if not in active TPC */
      extern double m_reco_vertex_dist_to_CPA;


      extern int m_reco_asso_showers;

//    extern   double m_reco_vertex_to_nearest_dead_wire_plane0;
//    extern   double m_reco_vertex_to_nearest_dead_wire_plane1;
//    extern   double m_reco_vertex_to_nearest_dead_wire_plane2;

      //added eventweight
      //-------------- EventWeight related variables -------------
      static const int k_max_mc_particles=100;

      extern int m_run_number_eventweight;
      extern int m_subrun_number_eventweight;
      extern int m_event_number_eventweight;

      extern double m_mcflux_nu_pos_x;
      extern double m_mcflux_nu_pos_y;
      extern double m_mcflux_nu_pos_z;
      extern double m_mcflux_nu_mom_x;
      extern double m_mcflux_nu_mom_y;
      extern double m_mcflux_nu_mom_z;
      extern double m_mcflux_nu_mom_E;
      extern int m_mcflux_ntype;
      extern int m_mcflux_ptype;
      extern double m_mcflux_nimpwt;
      extern double m_mcflux_dk2gen;
      extern double m_mcflux_nenergyn;
      extern double m_mcflux_tpx;
      extern double m_mcflux_tpy;
      extern double m_mcflux_tpz;
      extern double m_mcflux_vx;
      extern double m_mcflux_vy;
      extern double m_mcflux_vz;
      extern int m_mcflux_tptype;
      extern int m_mctruth_nparticles;
      extern int m_mctruth_particles_track_Id[k_max_mc_particles];
      extern int m_mctruth_particles_pdg_code[k_max_mc_particles];
      extern int m_mctruth_particles_mother[k_max_mc_particles];
      extern int m_mctruth_particles_status_code[k_max_mc_particles];
      extern int m_mctruth_particles_num_daughters[k_max_mc_particles]; //other similar variables
      extern int m_mctruth_particles_daughters[100][100];
      extern double m_mctruth_particles_Gvx[k_max_mc_particles];
      extern double m_mctruth_particles_Gvy[k_max_mc_particles];
      extern double m_mctruth_particles_Gvz[k_max_mc_particles];
      extern double m_mctruth_particles_Gvt[k_max_mc_particles];
      extern double m_mctruth_particles_px0[k_max_mc_particles];
      extern double m_mctruth_particles_py0[k_max_mc_particles];
      extern double m_mctruth_particles_pz0[k_max_mc_particles];
      extern double m_mctruth_particles_e0[k_max_mc_particles];
      extern int m_mctruth_particles_rescatter[k_max_mc_particles];
      extern double m_mctruth_particles_polx[k_max_mc_particles];
      extern double m_mctruth_particles_poly[k_max_mc_particles];
      extern double m_mctruth_particles_polz[k_max_mc_particles];
      extern int m_mctruth_neutrino_ccnc;
      extern int m_mctruth_neutrino_mode;
      extern int m_mctruth_neutrino_interaction_type;
      extern int m_mctruth_neutrino_target;
      extern int m_mctruth_neutrino_nucleon;
      extern int m_mctruth_neutrino_quark;
      extern double m_mctruth_neutrino_w;
      extern double m_mctruth_neutrino_x;
      extern double m_mctruth_neutrino_y;
      extern double m_mctruth_neutrino_qsqr;
      extern bool m_gtruth_is_sea_quark;
      extern int m_gtruth_tgt_pdg;
      extern int m_gtruth_tgt_Z;
      extern int m_gtruth_tgt_A;
      extern double m_gtruth_tgt_p4_x;
      extern double m_gtruth_tgt_p4_y;
      extern double m_gtruth_tgt_p4_z;
      extern double m_gtruth_tgt_p4_E;
      extern double m_gtruth_weight;
      extern double m_gtruth_probability;
      extern double m_gtruth_xsec;
      extern double m_gtruth_diff_xsec;
      extern int m_gtruth_gphase_space;
      extern double m_gtruth_vertex_x;
      extern double m_gtruth_vertex_y;
      extern double m_gtruth_vertex_z;
      extern double m_gtruth_vertex_T;
      extern int m_gtruth_gscatter;
      extern int m_gtruth_gint;
      extern int m_gtruth_res_num;
      extern int m_gtruth_num_piplus;
      extern int m_gtruth_num_pi0;
      extern int m_gtruth_num_piminus;
      extern int m_gtruth_num_proton;
      extern int m_gtruth_num_neutron;
      extern bool m_gtruth_is_charm;
      extern bool m_gtruth_is_strange;
      extern int m_gtruth_charm_hadron_pdg;
      extern int m_gtruth_strange_hadron_pdg;
      extern int m_gtruth_decay_mode;
      extern double m_gtruth_gx;
      extern double m_gtruth_gy;
      extern double m_gtruth_gt;
      extern double m_gtruth_gw;
      extern double m_gtruth_gQ2;
      extern double m_gtruth_gq2;
      extern int m_gtruth_probe_pdg;
      extern double m_gtruth_probe_p4_x;
      extern double m_gtruth_probe_p4_y;
      extern double m_gtruth_probe_p4_z;
      extern double m_gtruth_probe_p4_E;
      extern double m_gtruth_hit_nuc_p4_x;
      extern double m_gtruth_hit_nuc_p4_y;
      extern double m_gtruth_hit_nuc_p4_z;
      extern double m_gtruth_hit_nuc_p4_E;
      extern double m_gtruth_hit_nuc_pos;
      extern double m_gtruth_fs_had_syst_p4_x;
      extern double m_gtruth_fs_had_syst_p4_y;
      extern double m_gtruth_fs_had_syst_p4_z;
      extern double m_gtruth_fs_had_syst_p4_E;

      //-------------- Flash related variables -------------
      extern int m_reco_num_templates;
      extern std::vector<double> m_reco_template;  /* temp comment: does not seem to be used */


      //-------------- Flash related variables -------------
      extern std::vector<double> m_reco_flash_total_pe;
      extern std::vector<double> m_reco_flash_time;
      extern std::vector<double> m_reco_flash_time_width;
      extern std::vector<double> m_reco_flash_abs_time;
      extern std::vector<int>    m_reco_flash_frame;
      extern std::vector<double> m_reco_flash_ycenter;
      extern std::vector<double> m_reco_flash_ywidth;
      extern std::vector<double> m_reco_flash_zcenter;
      extern std::vector<double> m_reco_flash_zwidth;
      extern std::vector<double> m_reco_flash_total_pe_in_beamgate;
      extern std::vector<double> m_reco_flash_time_in_beamgate;
      extern std::vector<double> m_reco_flash_ycenter_in_beamgate;
      extern std::vector<double> m_reco_flash_zcenter_in_beamgate;

      extern int m_reco_num_flashes;
      extern int m_reco_num_flashes_in_beamgate;

      extern double m_beamgate_flash_start;
      extern double m_beamgate_flash_end;


      //----------- CRT related variables -----------------

      //for crt hits from the CRT veto product
      extern int m_CRT_veto_nhits;  /* number of CRT veto hits */
      extern std::vector<double> m_CRT_veto_hit_PE;  

      //fields storing information about the CRT hit closest to the flash
      extern double m_CRT_min_hit_time;
      extern double m_CRT_min_hit_PE;
      extern double m_CRT_min_hit_x;
      extern double m_CRT_min_hit_y;
      extern double m_CRT_min_hit_z;

      //Fields storing information about all CRT hits in event
      extern std::vector<double> m_CRT_hits_time;
      extern std::vector<double> m_CRT_hits_PE;
      extern std::vector<double> m_CRT_hits_x;
      extern std::vector<double> m_CRT_hits_y;
      extern std::vector<double> m_CRT_hits_z;
      extern double m_CRT_dt; //time between flash and nearest CRT hit

      extern double m_genie_spline_weight;
      extern double m_genie_CV_tune_weight;

      extern double m_photonu_weight_low;
      extern double m_photonu_weight_high;


      extern int pfp_w_bestnuID;     
    //------------ Track related Variables -------------
      extern double m_track_calo_min_dEdx;
      extern double m_track_calo_max_dEdx;
      extern double m_track_calo_min_dEdx_hits;
      extern double m_track_calo_trunc_fraction;

      extern int m_reco_asso_tracks;  /* number of track. (temp: will figure out what associate means later) */
      extern std::vector<int> m_reco_track_num_daughters;
      extern std::vector<double> m_reco_track_daughter_trackscore;  /* track score of this reco track's first daughter */
      extern std::vector<double> m_reco_track_length;  /* whole length of the reco track */
      extern std::vector<double> m_reco_track_dirx;    /* need to understand what the pair track->Direction() returns means*/
      extern std::vector<double> m_reco_track_diry;
      extern std::vector<double> m_reco_track_dirz;
      extern std::vector<double> m_reco_track_startx;  /* start pos of the track in cartesian X */
      extern std::vector<double> m_reco_track_starty;
      extern std::vector<double> m_reco_track_startz;
      extern std::vector<double> m_reco_track_endx;    /* end of the track in cartesian X */
      extern std::vector<double> m_reco_track_endy;
      extern std::vector<double> m_reco_track_endz;
      extern std::vector<double> m_reco_track_end_dist_to_active_TPC;  /* min distance from track end to TPC active volume boundaries */
      extern std::vector<double> m_reco_track_start_dist_to_active_TPC; /* min dist from trk start to TPC active boundaries */
      extern std::vector<double> m_reco_track_end_dist_to_CPA;
      extern std::vector<double> m_reco_track_start_dist_to_CPA;
      extern std::vector<double> m_reco_track_end_dist_to_SCB;          /* min dist from track end to SCB */
      extern std::vector<double> m_reco_track_start_dist_to_SCB;
      extern std::vector<int> m_reco_track_end_in_SCB;   /* if track end is in SCB boundary, 1- yes, 0- no */
      extern std::vector<int> m_reco_track_start_in_SCB;
      extern std::vector<double> m_reco_track_calo_energy_plane0;  /* energy sum of hits on plane 0 that correspond to the reco track */
      extern std::vector<double> m_reco_track_calo_energy_plane1;
      extern std::vector<double> m_reco_track_calo_energy_plane2;
      extern std::vector<double> m_reco_track_calo_energy_max;    /* max energy of 3 plane for the reco track */

      extern std::vector<double>   m_reco_track_theta_yz; /* theta, phi of the track */
      extern std::vector<double>   m_reco_track_phi_yx;

      extern std::vector<int> m_reco_track_num_trajpoints;  /* number of valid points in the track */
      extern std::vector<int> m_reco_track_num_spacepoints;  /* number of recob::spacepoints coresponding to the reco track */
      extern std::vector<double> m_reco_track_proton_kinetic_energy; /* energy of the track, under the asssumption it's a proton track 
                                   * set to -9999 if m_run_pi0_filter is set to true */

      extern std::vector<size_t>  m_reco_track_ordered_energy_index; /* index of m_reco_track_proton_kinetic_energy such that element values are in descending order */
      extern std::vector<size_t>  m_reco_track_ordered_displacement_index; /* index of m_reco_track_length so that track length are in descending order */
      extern std::vector<double> m_reco_track_spacepoint_principal0; /* PCA of reco track (in 3D spacepoint) */
      extern std::vector<double> m_reco_track_spacepoint_principal1;
      extern std::vector<double> m_reco_track_spacepoint_principal2;

      extern std::vector<double> m_reco_track_spacepoint_chi;  /* essentially sum of square of distances between spacepoint and the track line*/
      extern std::vector<double> m_reco_track_spacepoint_max_dist; /* max distance between a track and its coresponding spacepoints */


      // ---- corresponding variables on the best plane of reco track, which is defined as such------
      // if plane 2 have good hits, then plane 2 is the best-plane
      // otherwise, which plane of plane 0 and 1 has more good hits will be best plane
      // if none of 3 planes has good hits, then best-plane is set to -1
      extern std::vector<int> m_reco_track_best_calo_plane;
      extern std::vector<double> m_reco_track_mean_dEdx_best_plane;
      extern std::vector<double> m_reco_track_mean_dEdx_start_half_best_plane;
      extern std::vector<double> m_reco_track_mean_dEdx_end_half_best_plane;
      extern std::vector<int> m_reco_track_good_calo_best_plane;
      extern std::vector<std::vector<double>> m_reco_track_trunc_dEdx_best_plane;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_best_plane;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_best_plane;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_best_plane;
      extern std::vector<double> m_reco_track_trunc_PIDA_best_plane;
      extern std::vector<std::vector<double>> m_reco_track_resrange_best_plane;
      extern std::vector<std::vector<double>> m_reco_track_dEdx_best_plane;

      extern std::vector<double> m_reco_track_mean_dEdx_p0;  /* mean dEdx of hits on plane 0 of the reco track */
      extern std::vector<double> m_reco_track_mean_dEdx_start_half_p0; /* mean dEdx of first half of the track */
      extern std::vector<double> m_reco_track_mean_dEdx_end_half_p0;
      extern std::vector<int> m_reco_track_good_calo_p0; /* number of good dEdx hits on plane 0 of track calorimetry */
      extern std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p0;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_p0;  /* mean of truncated dEdx's of good hits */
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p0; /*mean of first half of trucated dEdx's of good hits */
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p0;
      extern std::vector<double> m_reco_track_trunc_PIDA_p0; /* mean of constant A in residual range formula, calc'd from good hits */
      extern std::vector<std::vector<double>> m_reco_track_resrange_p0; /* vec of residual range of good hits per reco track */
      extern std::vector<std::vector<double>> m_reco_track_dEdx_p0;     /* vec of dEdx of good hits per reco track */

      extern std::vector<double> m_reco_track_mean_dEdx_p1;
      extern std::vector<double> m_reco_track_mean_dEdx_start_half_p1;
      extern std::vector<double> m_reco_track_mean_dEdx_end_half_p1;
      extern std::vector<int> m_reco_track_good_calo_p1;
      extern std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p1;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_p1;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p1;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p1;
      extern std::vector<double> m_reco_track_trunc_PIDA_p1;
      extern std::vector<std::vector<double>> m_reco_track_resrange_p1;
      extern std::vector<std::vector<double>> m_reco_track_dEdx_p1;

      extern std::vector<double> m_reco_track_mean_dEdx_p2;
      extern std::vector<double> m_reco_track_mean_dEdx_start_half_p2;
      extern std::vector<double> m_reco_track_mean_dEdx_end_half_p2;
      extern std::vector<int> m_reco_track_good_calo_p2;
      extern std::vector<std::vector<double>> m_reco_track_trunc_dEdx_p2;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_p2;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_start_half_p2;
      extern std::vector<double> m_reco_track_mean_trunc_dEdx_end_half_p2;
      extern std::vector<double> m_reco_track_trunc_PIDA_p2;
      extern std::vector<std::vector<double>> m_reco_track_resrange_p2;
      extern std::vector<std::vector<double>> m_reco_track_dEdx_p2;

      extern std::vector<int> m_reco_track_num_calo_hits_p0; /* number of hits in calorimetry on plane 0 of each reco track */
      extern std::vector<int> m_reco_track_num_calo_hits_p1;
      extern std::vector<int> m_reco_track_num_calo_hits_p2;

//      std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane0; /* distance between track end and the nearest dead wire on plane*/
//      std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane1;
//      std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane2;
      
      extern std::vector<int> m_reco_track_sliceId; //the slice id for the slice continaing the reco track
      extern std::vector<double> m_reco_track_nuscore; //the neutrino score of the slice containing the reco track
      extern std::vector<bool> m_reco_track_isclearcosmic;//true if reco track is in a clear cosmic slice
      extern std::vector<double> m_reco_track_trackscore; /* track score of reco track, -999 if track is not found in PFPToTrackScoreMap */
      extern std::vector<int> m_reco_track_pfparticle_pdg; /* PDG of track's corresponding PFParticle, -999 if track is not found in PFPToTrackScoreMap*/
      extern std::vector<bool> m_reco_track_is_nuslice;  /* if reco track is in a neutrino slice */
      
      
      
      
      extern std::vector<int> m_sim_track_matched;  /* if reco track has been matched to a MCParticle, 1-YES, 0-NO */
      
      //-------- energy, mass, pdg ..etc.. of the matched MCParticle of reco track -----
      extern std::vector<double> m_sim_track_overlay_fraction;
      extern std::vector<double> m_sim_track_energy;
      extern std::vector<double> m_sim_track_mass;
      extern std::vector<double> m_sim_track_kinetic_energy;
      extern std::vector<int> m_sim_track_pdg;
      extern std::vector<int> m_sim_track_parent_pdg;
      
      /* event origin types:
       * kUnknown: ???  
       * kBeamNeutrino: Beam neutrinos.
       * kCosmicRay: Cosmic rays.
       * kSuperNovaNeutrino: Supernova neutrinos.
       * kSingleParticle: single particles thrown at the detector
       */
      extern std::vector<int> m_sim_track_origin;   /* truth origin of the matched MCParticle */
      extern std::vector<std::string> m_sim_track_process;
      extern std::vector<double> m_sim_track_startx;  /* space-charge corrected start point of the match MCParticle */
      extern std::vector<double> m_sim_track_starty;
      extern std::vector<double> m_sim_track_startz;
      extern std::vector<double> m_sim_track_px;
      extern std::vector<double> m_sim_track_py;
      extern std::vector<double> m_sim_track_pz;
      extern std::vector<double> m_sim_track_endx;  /* space-charge corrected end-point of the matched MCParticle */
      extern std::vector<double> m_sim_track_endy;
      extern std::vector<double> m_sim_track_endz;
      extern std::vector<double> m_sim_track_length; /* track length calculated based on the SC-corrected start and end point of the matched MCParticle */

      extern std::vector<int> m_sim_track_trackID;
      //-------- energy, mass, pdg ..etc.. of the matched MCParticle of reco track -----
      
      
      
      extern std::vector<int> m_sim_track_sliceId; //the slice id for the slice continaing the sim track, based on corresponding recob:PFP
      extern std::vector<double> m_sim_track_nuscore; //the neutrino score of the slice containing the sim track
      extern std::vector<bool> m_sim_track_isclearcosmic;//true if sim track is in a clear cosmic slice
      
      /*-------------------------------------------------------------------------------------*/
      extern std::vector<double> m_isolation_min_dist_trk_shr; /* minimum distance betwee shower hits and track hits on each plane 
      * if there is no shower hits, set to 999
      * if there is shower hits but no track hits, set to -999
      */  
      extern std::vector<double> m_isolation_nearest_shr_hit_to_trk_wire; /* the wire number of shower hit closest to track hits */  
      extern std::vector<double> m_isolation_nearest_shr_hit_to_trk_time; /* the time tick of shower hit closest to track hits in the slice */
      
      
      extern std::vector<double> m_isolation_num_shr_hits_win_1cm_trk; /* number of shower hits whose min distance to track hits <= 1cm 
      * of each plane (this is a 3 element vector) */      
      extern std::vector<double> m_isolation_num_shr_hits_win_2cm_trk;      
      extern std::vector<double> m_isolation_num_shr_hits_win_5cm_trk;      
      extern std::vector<double> m_isolation_num_shr_hits_win_10cm_trk;      
      
      
      extern std::vector<double> m_isolation_min_dist_trk_unassoc; /* of all unassociated hits, min distance to closest track hits 
      * set to -999 if there is no unassociated hits or track hits on plane
      */
      extern std::vector<double> m_isolation_nearest_unassoc_hit_to_trk_wire;/* wire number of the unassociated hit that of all is nearest to track hits in the slice */  
      extern std::vector<double> m_isolation_nearest_unassoc_hit_to_trk_time; /* time tick of the unasso hit that is nearest to track hits in the slice */
      extern std::vector<double> m_isolation_num_unassoc_hits_win_1cm_trk; /* number of unasso hits whose min distance to track hits <= 1cm
      * on each plane (this vector has 3 elements) */ 
      extern std::vector<double> m_isolation_num_unassoc_hits_win_2cm_trk;      
      extern std::vector<double> m_isolation_num_unassoc_hits_win_5cm_trk;      
      extern std::vector<double> m_isolation_num_unassoc_hits_win_10cm_trk;      
      
      
      /*-------------------------------------------------------------------------------------*/
      //------------ Shower related Variables  -------------
      
      extern std::vector<int> m_reco_shower_num_daughters;
      extern std::vector<double> m_reco_shower_daughter_trackscore;

      extern std::vector<int>   m_reco_shower3d_exists;
      extern std::vector<double>   m_reco_shower3d_startx;
      extern std::vector<double>   m_reco_shower3d_starty;
      extern std::vector<double>   m_reco_shower3d_startz;
      extern std::vector<double>   m_reco_shower3d_dirx;
      extern std::vector<double>   m_reco_shower3d_diry;
      extern std::vector<double>   m_reco_shower3d_dirz;
      extern std::vector<double>   m_reco_shower3d_theta_yz; /* theta, phi of the 3D shower (direction) */
      extern std::vector<double>   m_reco_shower3d_phi_yx;

      extern std::vector<double> m_reco_shower3d_openingangle;
      extern std::vector<double> m_reco_shower3d_length;
      extern std::vector<double> m_reco_shower3d_conversion_distance;

      extern std::vector<double>   m_reco_shower3d_impact_parameter;  /* distance between vertex and 3D shower direction */
      extern std::vector<double>    m_reco_shower3d_implied_dirx; /* X component of the unit vector point from vertex to 3D shower start */
      extern std::vector<double>     m_reco_shower3d_implied_diry;
      extern std::vector<double>     m_reco_shower3d_implied_dirz;

      extern std::vector<double> m_reco_shower3d_energy_plane0;
      extern std::vector<double> m_reco_shower3d_energy_plane1;
      extern std::vector<double> m_reco_shower3d_energy_plane2;

      extern std::vector<double> m_reco_shower3d_dEdx_plane0;
      extern std::vector<double> m_reco_shower3d_dEdx_plane1;
      extern std::vector<double> m_reco_shower3d_dEdx_plane2;
      
      
      
      extern std::vector<double>   m_reco_shower_startx;
      extern std::vector<double>   m_reco_shower_starty;
      extern std::vector<double>   m_reco_shower_startz;
      extern std::vector<double> m_reco_shower_start_dist_to_active_TPC; /* distance from shower start to closest TPC wall */
      extern std::vector<double> m_reco_shower_start_dist_to_CPA; /* distance from shower start to closest TPC wall */
      extern std::vector<double> m_reco_shower_start_dist_to_SCB;
      extern std::vector<int> m_reco_shower_start_in_SCB;
      extern std::vector<double> m_reco_shower_end_dist_to_active_TPC;
      extern std::vector<double> m_reco_shower_end_dist_to_SCB;

      extern std::vector<double>   m_reco_shower_dirx; /* X component of shower direction */
      extern std::vector<double>   m_reco_shower_diry;
      extern std::vector<double>   m_reco_shower_dirz;
      extern std::vector<double>   m_reco_shower_theta_yz; /* theta, phi of the shower direction */
      extern std::vector<double>   m_reco_shower_phi_yx;

      extern std::vector<double> m_reco_shower_openingangle;
      extern std::vector<double> m_reco_shower_length;
      extern std::vector<double> m_reco_shower_conversion_distance;  /* distance between shower start and vertex */

      extern std::vector<double>   m_reco_shower_impact_parameter; /* distance from vertex to the shower direction */
      extern std::vector<double>    m_reco_shower_implied_dirx; /* the X component of the unit vector pointing from vertex to shower start */
      extern std::vector<double>     m_reco_shower_implied_diry;
      extern std::vector<double>     m_reco_shower_implied_dirz;

      extern std::vector<int> m_reco_shower_delaunay_num_triangles_plane0; /* num of delaunay triangles found on plane 0 for each shower */
      extern std::vector<int> m_reco_shower_delaunay_num_triangles_plane1;
      extern std::vector<int> m_reco_shower_delaunay_num_triangles_plane2;
//
//      std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane0;/* dist from shower start to nearest dead wire on plane 0 */
//      std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane1;
//      std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane2;
      
      
      //shower flash matching
      
      extern std::vector<double> m_reco_shower_flash_shortest_distz;
      extern std::vector<double> m_reco_shower_flash_shortest_disty;
      extern std::vector<double> m_reco_shower_flash_shortest_distyz;

      extern std::vector<int> m_reco_shower_flash_shortest_index_z;
      extern std::vector<int> m_reco_shower_flash_shortest_index_y;
      extern std::vector<int> m_reco_shower_flash_shortest_index_yz;

      extern double  m_flash_optfltr_pe_beam;
      extern double  m_flash_optfltr_pe_beam_tot;
      extern double  m_flash_optfltr_pe_veto;
      extern double  m_flash_optfltr_pe_veto_tot;
      
      //end flash matching
      extern std::vector<int> m_reco_shower_num_hits_plane0; /* number of hits on plane 0 for each shower */
      extern std::vector<int> m_reco_shower_num_hits_plane1;
      extern std::vector<int> m_reco_shower_num_hits_plane2;
      extern std::vector<double> m_reco_shower_delaunay_area_plane0; /* total area of delaunay triangles found on plane 0 for each shower */
      extern std::vector<double> m_reco_shower_delaunay_area_plane1;
      extern std::vector<double> m_reco_shower_delaunay_area_plane2;
      extern std::vector<int> m_reco_shower_sliceId; //the slice id for the slice continaing the reco shower
      extern std::vector<double> m_reco_shower_nuscore; //the neutrino score of the slice containing the reco shower
      extern std::vector<bool> m_reco_shower_isclearcosmic;//true if reco shower is in a clear cosmic slice
      extern std::vector<bool> m_reco_shower_is_nuslice;//true if reco shower is in a clear cosmic slice
      extern std::vector<double> m_reco_shower_trackscore;
      extern std::vector<double> m_reco_shower_pfparticle_pdg;
      extern std::vector<double> m_reco_shower_kalman_exists; /* if there is a kalman track and reco::Calo related to this shower - 0, 1 */
      extern std::vector<double>   m_reco_shower_kalman_median_dEdx_plane0;
      extern std::vector<double>     m_reco_shower_kalman_median_dEdx_plane1;
      extern std::vector<double>   m_reco_shower_kalman_median_dEdx_plane2;
      extern std::vector<double>   m_reco_shower_kalman_median_dEdx_allplane;
      extern std::vector<double>      m_reco_shower_kalman_mean_dEdx_plane0;
      extern std::vector<double>    m_reco_shower_kalman_mean_dEdx_plane1;
      extern std::vector<double>    m_reco_shower_kalman_mean_dEdx_plane2;
      
      
      
      
      extern std::vector<int> m_sim_shower_matched;  /* whether shower has been matched to a MCParticle, 0 - False, 1 - True */
     
     // ----- energy, mass, pdg ... of the best-matched MCParticle for the shower ------
      extern std::vector<double> m_sim_shower_energy;
      extern std::vector<double> m_sim_shower_kinetic_energy;
      extern std::vector<double> m_sim_shower_mass;
      extern std::vector<int> m_sim_shower_pdg;
      extern std::vector<int> m_sim_shower_trackID;
      extern std::vector<int> m_sim_shower_parent_pdg;
      extern std::vector<int> m_sim_shower_parent_trackID;
      extern std::vector<int> m_sim_shower_origin;
      extern std::vector<std::string> m_sim_shower_process;
      extern std::vector<std::string> m_sim_shower_end_process;
      // ----- energy, mass, pdg ... of the best-matched MCParticle for the shower ------
      
      
      extern std::vector<double> m_sim_shower_start_x;  /* space charge corrected shower starting point */
      extern std::vector<double> m_sim_shower_start_y;
      extern std::vector<double> m_sim_shower_start_z;
      extern std::vector<double> m_sim_shower_vertex_x;  /* spacecharge corrected shower vertex */
      extern std::vector<double> m_sim_shower_vertex_y;
      extern std::vector<double> m_sim_shower_vertex_z;
      
      extern std::vector<double> m_sim_shower_px;
      extern std::vector<double> m_sim_shower_py;
      extern std::vector<double> m_sim_shower_pz;
      extern std::vector<int> m_sim_shower_is_true_shower;
      extern std::vector<int> m_sim_shower_best_matched_plane;
      extern std::vector<double> m_sim_shower_matched_energy_fraction_plane0; /* fraction of energy of the best-matched mother for shower on 
      * plane 0 over all energy deposited on plane 0 by the shower */
      extern std::vector<double> m_sim_shower_matched_energy_fraction_plane1;
      extern std::vector<double> m_sim_shower_matched_energy_fraction_plane2;
      extern std::vector<double> m_sim_shower_overlay_fraction; /* fraction of hits from overlay over all hits in the shower */
      extern std::vector<int> m_sim_shower_sliceId; //the slice id for the slice continaing the sim shower matched to reco
      extern std::vector<double> m_sim_shower_nuscore; //the neutrino score of the slice containing the sim shower matched to reco
      extern std::vector<bool> m_sim_shower_isclearcosmic;//true if sim shower matched to reco is in a clear cosmic slice
      extern std::vector<bool> m_sim_shower_is_nuslice;//true if sim shower matched to reco is in a clear cosmic slice
      
      
      
      //------------ MCTruth related Variables  -------------
      extern int m_mctruth_num;
      extern int m_mctruth_origin;
      extern double m_mctruth_nu_E;
      extern double m_mctruth_nu_vertex_x;
      extern double m_mctruth_nu_vertex_y;
      extern double m_mctruth_nu_vertex_z;
      extern double m_mctruth_reco_vertex_dist;
      extern double m_mctruth_lepton_E;
      extern int m_mctruth_nu_pdg;
      extern int m_mctruth_lepton_pdg;
      extern int m_mctruth_mode ;
      extern int m_mctruth_interaction_type ;
      extern int m_mctruth_ccnc;
      extern double m_mctruth_qsqr;
      extern int m_mctruth_num_daughter_particles;
      extern std::vector<int> m_mctruth_daughters_pdg;
      extern std::vector<double> m_mctruth_daughters_E;
      extern std::vector<int> m_mctruth_daughters_status_code;
      extern std::vector<int> m_mctruth_daughters_trackID;
      extern std::vector<int> m_mctruth_daughters_mother_trackID;
      extern std::vector<double> m_mctruth_daughters_px;
      extern std::vector<double> m_mctruth_daughters_py;
      extern std::vector<double> m_mctruth_daughters_pz;
      extern std::vector<double> m_mctruth_daughters_startx;
      extern std::vector<double> m_mctruth_daughters_starty;
      extern std::vector<double> m_mctruth_daughters_startz;
      extern std::vector<double> m_mctruth_daughters_time;
      extern std::vector<double> m_mctruth_daughters_endx;
      extern std::vector<double> m_mctruth_daughters_endy;
      extern std::vector<double> m_mctruth_daughters_endz;
      extern std::vector<double> m_mctruth_daughters_endtime;
      extern std::vector<std::string> m_mctruth_daughters_process;
      extern std::vector<std::string> m_mctruth_daughters_end_process;
      extern int     m_mctruth_num_exiting_photons ;
      extern int      m_mctruth_num_exiting_protons ;
      extern int    m_mctruth_num_exiting_pi0 ;
      extern int   m_mctruth_num_exiting_pipm ;
      extern int   m_mctruth_num_exiting_neutrons; 
      extern int   m_mctruth_num_exiting_delta0; 
      extern int   m_mctruth_num_exiting_deltapm; 
      extern int   m_mctruth_num_exiting_deltapp; 
      extern double m_mctruth_leading_exiting_proton_energy;
      extern int m_mctruth_is_delta_radiative;
      extern int m_mctruth_delta_radiative_1g1p_or_1g1n;
      extern double m_mctruth_delta_photon_energy;
      extern double m_mctruth_delta_proton_energy;
      extern double m_mctruth_delta_neutron_energy;
      extern std::vector<int> m_mctruth_exiting_delta0_num_daughters;

      extern std::vector<int> m_mctruth_exiting_photon_trackID;
      extern std::vector<int> m_mctruth_exiting_photon_mother_trackID;
      extern std::vector<int> m_mctruth_exiting_photon_from_delta_decay;
      extern std::vector<double> m_mctruth_exiting_photon_energy;
      extern std::vector<double> m_mctruth_exiting_photon_px;
      extern std::vector<double> m_mctruth_exiting_photon_py;
      extern std::vector<double> m_mctruth_exiting_photon_pz;

      extern std::vector<int> m_mctruth_exiting_proton_trackID;
      extern std::vector<int> m_mctruth_exiting_proton_mother_trackID;
      extern std::vector<int> m_mctruth_exiting_proton_from_delta_decay;
      extern std::vector<double> m_mctruth_exiting_proton_energy;
      extern std::vector<double> m_mctruth_exiting_proton_px;
      extern std::vector<double> m_mctruth_exiting_proton_py;
      extern std::vector<double> m_mctruth_exiting_proton_pz;

      extern std::vector<int> m_mctruth_exiting_neutron_trackID;
      extern std::vector<int> m_mctruth_exiting_neutron_mother_trackID;
      extern std::vector<int> m_mctruth_exiting_neutron_from_delta_decay;
      extern std::vector<double> m_mctruth_exiting_neutron_energy;
      extern std::vector<double> m_mctruth_exiting_neutron_px;
      extern std::vector<double> m_mctruth_exiting_neutron_py;
      extern std::vector<double> m_mctruth_exiting_neutron_pz;
      extern int  m_mctruth_num_reconstructable_protons;
      extern bool  m_mctruth_is_reconstructable_1g1p;
      extern bool  m_mctruth_is_reconstructable_1g0p;

      extern std::vector<double>        m_mctruth_exiting_pi0_E;
      extern std::vector<double>        m_mctruth_exiting_pi0_mom;
      extern std::vector<double>        m_mctruth_exiting_pi0_px;
      extern std::vector<double>        m_mctruth_exiting_pi0_py;
      extern std::vector<double>        m_mctruth_exiting_pi0_pz;

      extern double m_mctruth_pi0_leading_photon_energy;
      extern std::string m_mctruth_pi0_leading_photon_end_process;
      extern double m_mctruth_pi0_subleading_photon_energy;
      extern std::string m_mctruth_pi0_subleading_photon_end_process;
      extern std::vector<double> m_mctruth_pi0_subleading_photon_end;
      extern std::vector<double> m_mctruth_pi0_subleading_photon_start;
      extern std::vector<double> m_mctruth_pi0_leading_photon_end;
      extern std::vector<double> m_mctruth_pi0_leading_photon_start;
      extern int    m_mctruth_pi0_leading_photon_exiting_TPC;
      extern int    m_mctruth_pi0_subleading_photon_exiting_TPC;
      extern std::vector<double> m_mctruth_pi0_leading_photon_mom;
      extern std::vector<double> m_mctruth_pi0_subleading_photon_mom;
      extern std::string  m_truthmatching_signaldef;

      //the calo calculated quantities 
      extern std::vector<double> m_reco_shower_energy_max; //for each hit in a shower, converts Q->E, and sums. The max energy of all planes
      extern std::vector<double> m_reco_shower_energy_plane0; /* shower energy (summed hit energy) on plan 0 */
      extern std::vector<double> m_reco_shower_energy_plane1;
      extern std::vector<double> m_reco_shower_energy_plane2;
      extern std::vector<double> m_reco_shower_reclustered_energy_max;
      extern std::vector<double> m_reco_shower_reclustered_energy_plane0; /* total energy of the reco shower, and unassociated hit clusters 
      * close enough to it */
      extern std::vector<double> m_reco_shower_reclustered_energy_plane1;
      extern std::vector<double> m_reco_shower_reclustered_energy_plane2;
      extern std::vector<double> m_reco_shower_plane0;
      extern std::vector<double> m_reco_shower_plane1;
      extern std::vector<double> m_reco_shower_plane2;
      extern std::vector<double> m_reco_shower_plane0_nhits; /* num of shower hits on plane 0 */
      extern std::vector<double> m_reco_shower_plane1_nhits;
      extern std::vector<double> m_reco_shower_plane2_nhits;
      extern std::vector<double> m_reco_shower_plane0_meanRMS; /* the mean of RMS of the shower hit shape (in tick unit) on plane 0 */
      extern std::vector<double> m_reco_shower_plane1_meanRMS;
      extern std::vector<double> m_reco_shower_plane2_meanRMS;
      extern std::vector<int> m_reco_shower_hit_wire;
      extern std::vector<int> m_reco_shower_hit_plane;
      extern std::vector<double> m_reco_shower_hit_tick;
      extern std::vector<double> m_reco_shower_spacepoint_x;
      extern std::vector<double> m_reco_shower_spacepoint_z;
      extern std::vector<double> m_reco_shower_spacepoint_y;
      extern std::vector<size_t>  m_reco_shower_ordered_energy_index; /* indices of 'm_reco_shower_energy_max' such that energy max is in descending order */
      extern std::vector<std::vector<double>> m_reco_shower_dQdx_plane0; //for each shower, looks at the hits for all clusters in the plane, stores the dQ/dx for each hit 
      extern std::vector<std::vector<double>> m_reco_shower_dQdx_plane1;
      extern std::vector<std::vector<double>> m_reco_shower_dQdx_plane2;
      extern std::vector<std::vector<double>> m_reco_shower_dEdx_plane0; //dE/dx from the calculated dQ/dx for each hit of all clusters in shower on plane   
      extern std::vector<std::vector<double>> m_reco_shower_dEdx_plane1;
      extern std::vector<std::vector<double>> m_reco_shower_dEdx_plane2;
      extern std::vector<double> m_reco_shower_dEdx_plane0_mean; /* mean of dE/dx of each hit in shower */
      extern std::vector<double> m_reco_shower_dEdx_plane1_mean;
      extern std::vector<double> m_reco_shower_dEdx_plane2_mean;
      extern std::vector<double> m_reco_shower_dEdx_plane0_max;
      extern std::vector<double> m_reco_shower_dEdx_plane1_max;
      extern std::vector<double> m_reco_shower_dEdx_plane2_max;
      extern std::vector<double> m_reco_shower_dEdx_plane0_min;
      extern std::vector<double> m_reco_shower_dEdx_plane1_min;
      extern std::vector<double> m_reco_shower_dEdx_plane2_min;
      extern std::vector<double> m_reco_shower_dEdx_plane0_median;/* median of dE/dx of each hit in shower (median of vector element of m_reco_shower_dEdx_plane0) */
      extern std::vector<double> m_reco_shower_dEdx_plane1_median;
      extern std::vector<double> m_reco_shower_dEdx_plane2_median;
      extern std::vector<double>  m_reco_shower_angle_wrt_wires_plane0; /* angle between shower direction and wire dir on plane, in radian*/
      extern std::vector<double>  m_reco_shower_angle_wrt_wires_plane1;
      extern std::vector<double>  m_reco_shower_angle_wrt_wires_plane2;
      extern std::vector<double>  m_reco_shower_dEdx_amalgamated;
      extern std::vector<int>  m_reco_shower_dEdx_amalgamated_nhits;
      extern std::vector<double> m_reco_shower_dQdx_plane0_median;/* median of dQ/dx of each hit in shower (median of m_reco_shower_dQdx_plane0) */
      extern std::vector<double> m_reco_shower_dQdx_plane1_median;
      extern std::vector<double> m_reco_shower_dQdx_plane2_median;
      extern std::vector<double> m_reco_shower_dEdx_plane0_nhits; /* number of hits of all clusters of the shower on plane 0 */
      extern std::vector<double> m_reco_shower_dEdx_plane1_nhits;
      extern std::vector<double> m_reco_shower_dEdx_plane2_nhits;
      extern double _time2cm;//value modeled from David's shower code
      // PID-related variables
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane0;
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane1;
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_mu_plane2;
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane0;
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane1;
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_p_plane2;
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane0;
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane1;
      extern std::vector<double> m_reco_track_pid_bragg_likelihood_mip_plane2;
      extern std::vector<double> m_reco_track_pid_pida_plane0;
      extern std::vector<double> m_reco_track_pid_pida_plane1;
      extern std::vector<double> m_reco_track_pid_pida_plane2;
      extern std::vector<double> m_reco_track_pid_chi2_mu_plane0;
      extern std::vector<double> m_reco_track_pid_chi2_mu_plane1;
      extern std::vector<double> m_reco_track_pid_chi2_mu_plane2;
      extern std::vector<double> m_reco_track_pid_chi2_p_plane0;
      extern std::vector<double> m_reco_track_pid_chi2_p_plane1;
      extern std::vector<double> m_reco_track_pid_chi2_p_plane2;
      extern std::vector<double> m_reco_track_pid_three_plane_proton_pid;

      //Geant4
      extern std::vector<int> m_geant4_pdg;
      extern std::vector<int>          m_geant4_trackid;
      extern std::vector<int>          m_geant4_mother;
      extern std::vector<int>         m_geant4_statuscode;
      extern std::vector<double>          m_geant4_E;
      extern std::vector<double>          m_geant4_mass;
      extern std::vector<double>          m_geant4_px;
      extern std::vector<double>          m_geant4_py;
      extern std::vector<double>          m_geant4_pz;
      extern std::vector<double>          m_geant4_vx;
      extern std::vector<double>          m_geant4_vy;
      extern std::vector<double>          m_geant4_vz;
      extern std::vector<double>          m_geant4_dx;
      extern std::vector<double>          m_geant4_dy;
      extern std::vector<double>          m_geant4_dz;
      extern std::vector<std::string>          m_geant4_process;
      extern std::vector<std::string>          m_geant4_end_process;
      extern std::vector<double>          m_geant4_costheta;

      //matching variables
      extern int  m_reco_slice_num; //total number of slices in the event
      extern std::vector<double> m_reco_slice_nuscore; //vector of the neutrino score for each slice in an event
      extern int m_reco_slice_shower_num_matched_signal; //the number of sim showers matched an MCP in the signal def
      extern int m_reco_slice_track_num_matched_signal; //the number of sim showers matched an MCP in the signal def
      extern std::vector<int> m_reco_slice_shower_matched_sliceId; //the slice id for each matched shower
      extern std::vector<int> m_reco_slice_track_matched_sliceId; //the slice id for each matched track

      extern std::vector<int> m_reco_slice_num_pfps; //the total number of PFP's per slice
      extern std::vector<int> m_reco_slice_num_showers; //the subset of PFP's that are showers, ie number of showers per slice 
      extern std::vector<int> m_reco_slice_num_tracks; //the subset of PFP's that are tracks

      extern std::vector<double> m_reco_slice_shower_matched_energy; //the energy for each matched shower
      extern std::vector<double> m_reco_slice_track_matched_energy; //the energy for each matched track
      extern std::vector<double> m_reco_slice_shower_matched_conversion; //the conversion distance for each matched shower
      extern std::vector<double> m_reco_slice_shower_matched_overlay_frac; //fraction of overlay hits for each matched shower
      //std::map<art::Ptr<recob::PFParticle>, double > & pfParticleToNuScoreMap;//is filled during analyze slices


      //-------  matched shower: reco shower that matches to a primary photon + max energy of 3 plane > 20 + definition being ncdelta----- 
      extern std::vector<double> m_matched_signal_shower_overlay_fraction;
      //std::vector<double> m_matched_signal_shower_conversion_length;
      extern std::vector<double> m_matched_signal_shower_true_E;  /* energy of the best-matched MCparticle for the shower */
      extern std::vector<double> m_matched_signal_shower_nuscore; /* the neutrino score of the slice containing the reco shower */
      extern std::vector<int> m_matched_signal_shower_sliceId;    /* reco shower slice ID */
      extern std::vector<bool> m_matched_signal_shower_is_clearcosmic;
      extern int m_matched_signal_shower_num;  /* number of match showers (that has unique best-matched primary photon ?)  */
      extern std::vector<bool> m_matched_signal_shower_is_nuslice;
      extern std::vector<int> m_matched_signal_shower_tracks_in_slice; /* number of showers in the same slice as of this reco shower */
      extern std::vector<int> m_matched_signal_shower_showers_in_slice; /* number of tracks in the same slice as of this reco shower */


      //-------- for reco tracks that match to a primary proton ---------
      extern std::vector<double> m_matched_signal_track_true_E; /*  the true energy of matched MCparticle (proton) */
      extern std::vector<double> m_matched_signal_track_nuscore;  /* nu score of the slice containing the reco track */
      extern std::vector<int> m_matched_signal_track_sliceId;
      extern std::vector<bool> m_matched_signal_track_is_clearcosmic; /* if reco track is in clear cosmic slice */
      //  std::vector<bool> m_matched_signal_track_is_nuslice;
      extern std::vector<bool> m_matched_signal_track_is_nuslice;
      extern std::vector<int> m_matched_signal_track_tracks_in_slice; /* num of PFP that are tracks in the slice this reco track is in */
      extern std::vector<int> m_matched_signal_track_showers_in_slice;

      extern int m_matched_signal_track_num;  /* num of reco tracks matched to primary proton */ 

      //int m_matched_signal_total_num_slices;

      //---------for reco tracks that match to a primary proton ---------

      extern bool m_reco_1g1p_is_same_slice;
      extern bool m_reco_1g1p_is_multiple_slices;
      extern bool m_reco_1g1p_is_nuslice;
      extern bool m_reco_1g0p_is_nuslice;
      extern double m_reco_1g1p_nuscore;
      extern double  m_reco_1g0p_nuscore;
      extern bool m_is_matched_1g1p;
      extern bool m_is_matched_1g0p;
      extern bool m_no_matched_showers;
      extern bool m_multiple_matched_showers; //if there is more than 1 eligible shower (match to primary photon, pass energy threshold)
      extern bool m_multiple_matched_tracks; /* if there is more than 1 eligible track (match to primary proton) */


}
#endif // SBNCODE_SINGLEPHOTONANALYSIS_VARIABLES_H
