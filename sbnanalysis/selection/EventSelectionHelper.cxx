/**
 *
 * Make this a PostProcessorBase 
 *
 * This will take SBNCode processed events as input
 *
 * Can define ProcessEvent using just the SBNCode Event definition
 *
 * ---------------------------------------------------------------------------------------------
 *
 *  Author: R. Jones
 *  Email : <rsjones@fnal.gov>
 *  Date  : June 2019 
 *
 * ---------------------------------------------------------------------------------------------
 */

#include "EventSelectionHelper.hh"
#include "LoadEvents.hh"
#include <iostream>
#include <numeric>
#include "TLeaf.h"
#include "TBranch.h"
#include "TVector3.h"
#include <algorithm>
#include <iterator>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <stdexcept>

namespace selection{

  double EventSelectionHelper::GetPOT(TTree *subrun){
    // Get the pot for the individual file
    subrun->GetEntry(0);
    TBranch *b_pot = subrun->GetBranch("subrun_pot");
    return b_pot->GetLeaf("subrun_pot")->GetValue();
  }

  //------------------------------------------------------------------------------------------ 
 
  void EventSelectionHelper::GetTimeLeft(const int start_time, const int total, const unsigned int i){
    
    int now = static_cast<int>(time(NULL));
    int diff = now - start_time;
    int n_left = total - (i+1);
    int time_left = std::round(diff/double(i+1) * n_left);
    int seconds_per_minute = 60; 
    int seconds_per_hour = 3600; 
    int seconds_per_day = 86400;
    int days_left = time_left / (seconds_per_day);
    int hours_left = (time_left - (days_left * seconds_per_day)) / (seconds_per_hour);
    int minutes_left = (time_left - (days_left * seconds_per_day) - (hours_left * seconds_per_hour)) / (seconds_per_minute);
    int seconds_left = time_left - (days_left * seconds_per_day) - (hours_left * seconds_per_hour) - (minutes_left * seconds_per_minute);
    if(i%2){
      std::cout << " Estimated time left: " << std::setw(5) << days_left << " days, ";
      std::cout                             << std::setw(5) << hours_left << " hours, ";
      std::cout                             << std::setw(5) << minutes_left << " minutes, ";
      std::cout                             << std::setw(5) << seconds_left << " seconds." << '\r' << flush;
    }
  }

  //------------------------------------------------------------------------------------------ 
 
  void EventSelectionHelper::LoadEventList(const std::string &file_name, EventList &event_list, const int &file, double &pot){
 
    TFile f(file_name.c_str());
    if(!f.IsOpen()){
      std::cerr << " Error opening file " << std::endl;
      exit(1);
    }

    TTree *t_event    = (TTree*) f.Get("event_tree");
    TTree *t_subrun   = (TTree*) f.Get("subrun_tree");
    TTree *t_particle = (TTree*) f.Get("particle_tree");
    TTree *t_track    = (TTree*) f.Get("recotrack_tree");
    TTree *t_shower   = (TTree*) f.Get("recoshower_tree");

    TBranch *b_event_id        = t_event->GetBranch("event_id");
    TBranch *b_time_now        = t_event->GetBranch("time_now");
    TBranch *b_r_vertex        = t_event->GetBranch("r_vertex");
    TBranch *b_t_vertex        = t_event->GetBranch("t_vertex");
    TBranch *b_t_interaction   = t_event->GetBranch("t_interaction");
    TBranch *b_t_scatter       = t_event->GetBranch("t_scatter");
    TBranch *b_t_iscc          = t_event->GetBranch("t_iscc");
    TBranch *b_t_nu_pdgcode    = t_event->GetBranch("t_nu_pdgcode");
    TBranch *b_t_charged_pions = t_event->GetBranch("t_charged_pions");
    TBranch *b_t_neutral_pions = t_event->GetBranch("t_neutral_pions");
    TBranch *b_t_vertex_energy = t_event->GetBranch("t_vertex_energy");
    TBranch *b_t_neutrino_qsqr = t_event->GetBranch("t_qsqr");
    
    pot = EventSelectionHelper::GetPOT(t_subrun);
    
    unsigned int n_events = t_event->GetEntries();

    unsigned int start_tracks      = 0;
    unsigned int start_showers     = 0;
    unsigned int start_mcparticles = 0;

    for(unsigned int j = 0; j < n_events; ++j){

      ParticleList mcparticles;
      ParticleList recoparticles;
      TrackList    tracks;
      ShowerList   showers;

      TVector3 r_vertex, t_vertex;
      unsigned int interaction, pions_ch, pions_neu, scatter;
      int neutrino_pdg;
      bool iscc(false);
      float neu_energy;
      float neu_qsqr;

      t_event->GetEntry(j);

      int event_id = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now = b_time_now->GetLeaf("time_now")->GetValue();
      r_vertex[0]  = b_r_vertex->GetLeaf("r_vertex")->GetValue(0);
      r_vertex[1]  = b_r_vertex->GetLeaf("r_vertex")->GetValue(1);
      r_vertex[2]  = b_r_vertex->GetLeaf("r_vertex")->GetValue(2);
      t_vertex[0]  = b_t_vertex->GetLeaf("t_vertex")->GetValue(0);
      t_vertex[1]  = b_t_vertex->GetLeaf("t_vertex")->GetValue(1);
      t_vertex[2]  = b_t_vertex->GetLeaf("t_vertex")->GetValue(2);
      interaction  = b_t_interaction->GetLeaf("t_interaction")->GetValue();
      scatter      = b_t_scatter->GetLeaf("t_scatter")->GetValue();
      iscc         = b_t_iscc->GetLeaf("t_iscc")->GetValue();
      neutrino_pdg = b_t_nu_pdgcode->GetLeaf("t_nu_pdgcode")->GetValue();
      pions_ch     = b_t_charged_pions->GetLeaf("t_charged_pions")->GetValue();
      pions_neu    = b_t_neutral_pions->GetLeaf("t_neutral_pions")->GetValue();
      neu_energy   = b_t_vertex_energy->GetLeaf("t_vertex_energy")->GetValue();
      neu_qsqr     = b_t_neutrino_qsqr->GetLeaf("t_qsqr")->GetValue();
   
      std::pair<int,int> event_identification(event_id,time_now);

      EventSelectionHelper::GetTrackList(start_tracks, t_track, event_identification, tracks);
      EventSelectionHelper::GetShowerList(start_showers, t_shower, event_identification, showers);
      EventSelectionHelper::GetMCParticleList(start_mcparticles, t_particle, event_identification, mcparticles);
      
      if(tracks.size() != 0) EventSelectionHelper::GetRecoParticleFromTrack1EscapingDistanceCut(tracks, recoparticles);
      if(showers.size() != 0) EventSelectionHelper::GetRecoParticleFromShower(showers, r_vertex, recoparticles);
     
      // Check if any particles should be flipped
      EventSelectionHelper::CheckAndFlip(r_vertex, recoparticles);

      event_list.push_back(Event(mcparticles, recoparticles, interaction, scatter, neutrino_pdg, pions_ch, pions_neu, iscc, t_vertex, r_vertex, neu_energy, neu_qsqr, file, event_id));

      start_tracks      += tracks.size();
      start_showers     += showers.size();
      start_mcparticles += mcparticles.size();
    }
  }

  //------------------------------------------------------------------------------------------ 
  
  float EventSelectionHelper::GetDistanceToPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end){
    // Get the value of the unit vector of the particle dotted with the normal to the plane
    // If this is zero, the particle is parallel so throw an exception and catch it in the main
    // then continue looping through the list of planes
    TVector3 track_direction   = (end - vtx).Unit();
    float direction_from_plane = track_direction.Dot(plane.GetUnitN());

    if(std::abs(direction_from_plane) <= std::numeric_limits<float>::epsilon()){
      throw std::domain_error("The track is parallel to the plane");
    }

    return (1/direction_from_plane)*((plane.GetV() - vtx).Dot(plane.GetUnitN()));
  }

  //------------------------------------------------------------------------------------------ 

  float EventSelectionHelper::GetDistanceFromParticleToPlane(const Plane &plane, const Particle &particle){
    // Get the vertex and end of the particle and pass to the distance function
    TVector3 vtx = particle.GetVertex();
    TVector3 end = particle.GetEnd();
    return EventSelectionHelper::GetDistanceToPlane(plane, vtx, end);
  }

  //------------------------------------------------------------------------------------------ 

  float EventSelectionHelper::GetDistanceFromTrackToPlane(const Plane &plane, const Track &track){
    // Get the vertex and end of the particle and pass to the distance function
    TVector3 vtx = track.m_vertex;
    TVector3 end = track.m_end;
    return EventSelectionHelper::GetDistanceToPlane(plane, vtx, end);
  }

  //------------------------------------------------------------------------------------------ 
  
  bool EventSelectionHelper::CheckIfParticleIntersectsPlane(const Plane &plane, const Particle &particle){
    // Get the vertex and end of the particle and pass to the distance function
    TVector3 vtx = particle.GetVertex();
    TVector3 end = particle.GetEnd();
    float length = particle.GetLength();
    return EventSelectionHelper::CheckIfIntersectsPlane(plane, vtx, end, length);
  }
  
  //------------------------------------------------------------------------------------------ 
  
  bool EventSelectionHelper::CheckIfTrackIntersectsPlane(const Plane &plane, const Track &track){
    // Get the vertex and end of the particle and pass to the distance function
    TVector3 vtx = track.m_vertex;
    TVector3 end = track.m_end;
    float length = track.m_length;
    return EventSelectionHelper::CheckIfIntersectsPlane(plane, vtx, end, length);
  }

  //------------------------------------------------------------------------------------------ 
  
  bool EventSelectionHelper::CheckIfIntersectsPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end, const float &length){
    float d = -std::numeric_limits<float>::max();
    try{
      // Will throw exception if the particle is parallel to the plane
      d = GetDistanceToPlane(plane, vtx, end);
    }
    // If caught, the particle is parallel to the plane and does not intersect
    catch(const std::domain_error&) {return false;}

    if(d < 0 || d > length) return false;
    TVector3 track_direction    = (end - vtx).Unit();
    TVector3 intersection_point = vtx + d * track_direction; 

    return IsProjectedPointInPlaneBounds(intersection_point, plane);
  }

  //------------------------------------------------------------------------------------------ 

  bool EventSelectionHelper::IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane){
    // Check if the point lies within the bound plane
    return (std::abs((point-plane.GetV()).Dot(plane.GetUnitA())) <= plane.GetAlpha() && std::abs((point-plane.GetV()).Dot(plane.GetUnitB())) <= plane.GetBeta());
  }

  //------------------------------------------------------------------------------------------ 

  void EventSelectionHelper::CheckAndFlip(const TVector3 &vtx, ParticleList &particles){
  
    // Loop over reconstructed particles
    // Check if the end point of the particle is closer to the neutrino vertex than the start
    // Flip if true
    for(Particle &p : particles){
      // Make sure the particle we are looking at is a reconstructed track
      if(!p.GetFromRecoTrack()) continue;
    
      // If the end is closer to the neutrino vertex than the start, flip it
      float nu_vtx_dist = (vtx - p.GetVertex()).Mag(); 
      float nu_end_dist = (vtx - p.GetEnd()).Mag(); 
      if(nu_vtx_dist > nu_end_dist) p.FlipTrack();
    }
  }

  //------------------------------------------------------------------------------------------ 

  void EventSelectionHelper::GetSBNDAVPlanes(PlaneList &planes) {
    // Define the fiducial planes of the detector
    //  The order: +/-x, +/-y, +/-z
    //
    // To take into account the fiducial border define 
    //    Half with fiducial   = Half width  - width fiducial border 
    //    Half height fiducial = Half height - height fiducial border
    //    Half length fiducial = Half length - length fiducial border
    float w = 200; // half full detector width (x)
    float h = 200; // half full detector height (y)
    float l = 250; // half full detector length (z)
    float w_border = 0; // width fiducial border
    float h_border = 0; // height fiducial border
    float l_border = 0; // length fiducial border
    float fid_w = w - w_border; // fiducial half width
    float fid_h = h - h_border; // fiducial half height
    float fid_l = l - l_border; // fiducial half length

    // Define the ficudial planes of the detector
    planes.emplace_back(TVector3( fid_w, 0, l),         TVector3(0, fid_h, 0), TVector3(0, 0,  fid_l));
    planes.emplace_back(TVector3(-fid_w, 0, l),         TVector3(0, fid_h, 0), TVector3(0, 0, -fid_l));
    planes.emplace_back(TVector3(0,  fid_h, l),         TVector3(fid_w, 0, 0), TVector3(0, 0, -fid_l));
    planes.emplace_back(TVector3(0, -fid_h, l),         TVector3(fid_w, 0, 0), TVector3(0, 0,  fid_l));
    planes.emplace_back(TVector3(0, 0, l_border),       TVector3(fid_w, 0, 0), TVector3(0,  fid_h, 0));
    planes.emplace_back(TVector3(0, 0, (2*l)-l_border), TVector3(fid_w, 0, 0), TVector3(0, -fid_h, 0));
  }

  //------------------------------------------------------------------------------------------ 
 
  void EventSelectionHelper::GetSBNDFiducialPlanes(PlaneList &planes) {
    // Define the fiducial planes of the detector
    //  The order: +/-x, +/-y, +/-z
    //
    // To take into account the fiducial border define 
    //    Half with fiducial   = Half width  - width fiducial border 
    //    Half height fiducial = Half height - height fiducial border
    //    Half length fiducial = Half length - length fiducial border
    float w = 200; // half full detector width (x)
    float h = 200; // half full detector height (y)
    float l = 250; // half full detector length (z)
    float w_border = 10; // width fiducial border
    float h_border = 20; // height fiducial border
    float l_border = 10; // length fiducial border
    float fid_w = w - w_border; // fiducial half width
    float fid_h = h - h_border; // fiducial half height
    float fid_l = l - l_border; // fiducial half length

    // Define the ficudial planes of the detector
    planes.emplace_back(TVector3( fid_w, 0, l),         TVector3(0, fid_h, 0), TVector3(0, 0,  fid_l)); // 0 = +x
    planes.emplace_back(TVector3(-fid_w, 0, l),         TVector3(0, fid_h, 0), TVector3(0, 0, -fid_l)); // 1 = -x
    planes.emplace_back(TVector3(0,  fid_h, l),         TVector3(fid_w, 0, 0), TVector3(0, 0, -fid_l)); // 2 = +y
    planes.emplace_back(TVector3(0, -fid_h, l),         TVector3(fid_w, 0, 0), TVector3(0, 0,  fid_l)); // 3 = -y
    planes.emplace_back(TVector3(0, 0, l_border),       TVector3(fid_w, 0, 0), TVector3(0,  fid_h, 0)); // 4 = +z
    planes.emplace_back(TVector3(0, 0, (2*l)-l_border), TVector3(fid_w, 0, 0), TVector3(0, -fid_h, 0)); // 5 = -z
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionHelper::GetUniqueEventList(TTree *event_tree, UniqueEventIdList &unique_event_list){
  
    TBranch *b_event_id = event_tree->GetBranch("event_id");
    TBranch *b_time_now = event_tree->GetBranch("time_now");
    
    unsigned int n_events = event_tree->GetEntries();

    for(unsigned int i = 0; i < n_events; ++i){
   
      event_tree->GetEntry(i);

      int event_id_a = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now_a = b_time_now->GetLeaf("time_now")->GetValue();
  
      bool shouldAdd(true);
    
      for(unsigned int j = 0; j < n_events; ++j){
        event_tree->GetEntry(j);

        int event_id_b = b_event_id->GetLeaf("event_id")->GetValue();
        int time_now_b = b_time_now->GetLeaf("time_now")->GetValue();
        
        if (event_id_a == event_id_b && time_now_a == time_now_b && i != j){
          shouldAdd = false;
          break;
        }
      }
       
      if (shouldAdd)
        unique_event_list.push_back(std::pair<int, int>(event_id_a, time_now_a));
    }        
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionHelper::GetTrackList(unsigned int start, TTree *track_tree, const std::pair<int, int> &unique_event, TrackList &track_list){
   
    TBranch *b_event_id         = track_tree->GetBranch("event_id");
    TBranch *b_time_now         = track_tree->GetBranch("time_now");
    TBranch *b_id_charge        = track_tree->GetBranch("tr_id_charge");
    TBranch *b_id_energy        = track_tree->GetBranch("tr_id_energy");
    TBranch *b_id_hits          = track_tree->GetBranch("tr_id_hits");
    TBranch *b_n_hits           = track_tree->GetBranch("tr_n_hits");
    TBranch *b_vertex           = track_tree->GetBranch("tr_vertex");
    TBranch *b_end              = track_tree->GetBranch("tr_end");
    TBranch *b_pida             = track_tree->GetBranch("tr_pida");
    TBranch *b_chi2_mu          = track_tree->GetBranch("tr_chi2_mu");
    TBranch *b_chi2_pi          = track_tree->GetBranch("tr_chi2_pi");
    TBranch *b_chi2_pr          = track_tree->GetBranch("tr_chi2_pr");
    TBranch *b_chi2_ka          = track_tree->GetBranch("tr_chi2_ka");
    TBranch *b_length           = track_tree->GetBranch("tr_length");
    TBranch *b_kinetic_energy   = track_tree->GetBranch("tr_kinetic_energy");
    TBranch *b_mcs_mom_muon     = track_tree->GetBranch("tr_mcs_mom_muon");
    TBranch *b_range_mom_muon   = track_tree->GetBranch("tr_range_mom_muon");
    TBranch *b_range_mom_proton = track_tree->GetBranch("tr_range_mom_proton");
    TBranch *b_size             = track_tree->GetBranch("tr_dedx_size");      
    TBranch *b_residual_range   = track_tree->GetBranch("tr_residual_range"); 
    TBranch *b_dedx             = track_tree->GetBranch("tr_dedx");           
    
    unsigned int n_entries = track_tree->GetEntries();
    unsigned int n_dedx = 0;

    for(unsigned int i = 0; i < n_entries; ++i){
    
      track_tree->GetEntry(i);

      int event_id         = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now         = b_time_now->GetLeaf("time_now")->GetValue();
      
      if(event_id != unique_event.first || time_now != unique_event.second) continue;
      
      double temp_vertex[3];
      double temp_end[3];

      int id_charge          = b_id_charge->GetLeaf("tr_id_charge")->GetValue();
      int id_energy          = b_id_energy->GetLeaf("tr_id_energy")->GetValue();
      int id_hits            = b_id_hits->GetLeaf("tr_id_hits")->GetValue();
      int n_hits             = b_n_hits->GetLeaf("tr_n_hits")->GetValue();
      temp_vertex[0]         = b_vertex->GetLeaf("tr_vertex")->GetValue(0);
      temp_vertex[1]         = b_vertex->GetLeaf("tr_vertex")->GetValue(1);
      temp_vertex[2]         = b_vertex->GetLeaf("tr_vertex")->GetValue(2);
      temp_end[0]            = b_end->GetLeaf("tr_end")->GetValue(0);
      temp_end[1]            = b_end->GetLeaf("tr_end")->GetValue(1);
      temp_end[2]            = b_end->GetLeaf("tr_end")->GetValue(2);
      float pida             = b_pida->GetLeaf("tr_pida")->GetValue();
      float chi2_mu          = b_chi2_mu->GetLeaf("tr_chi2_mu")->GetValue();
      float chi2_pi          = b_chi2_pi->GetLeaf("tr_chi2_pi")->GetValue();
      float chi2_pr          = b_chi2_pr->GetLeaf("tr_chi2_pr")->GetValue();
      float chi2_ka          = b_chi2_ka->GetLeaf("tr_chi2_ka")->GetValue();
      float length           = b_length->GetLeaf("tr_length")->GetValue();
      float kinetic_energy   = b_kinetic_energy->GetLeaf("tr_kinetic_energy")->GetValue();
      float mcs_mom_muon     = b_mcs_mom_muon->GetLeaf("tr_mcs_mom_muon")->GetValue();
      float range_mom_muon   = b_range_mom_muon->GetLeaf("tr_range_mom_muon")->GetValue();
      float range_mom_proton = b_range_mom_proton->GetLeaf("tr_range_mom_proton")->GetValue();

      TVector3 vertex(temp_vertex);
      TVector3 end(temp_end);
      
      float vertex_x = temp_vertex[0];                        
      float vertex_y = temp_vertex[1];                        
      float vertex_z = temp_vertex[2];                        
      float end_x    = temp_end[0];                        
      float end_y    = temp_end[1];                        
      float end_z    = temp_end[2];                        
                                                                                   
      // Co-ordinate offset in cm
      int sbnd_length_x = 400;
      int sbnd_length_y = 400;
      int sbnd_length_z = 500;
      
      int sbnd_offset_x = 200;
      int sbnd_offset_y = 200;
      int sbnd_offset_z = 0;

      int sbnd_border_x = 10;
      int sbnd_border_y = 20;
      int sbnd_border_z = 10;

      bool not_contained = 
        (     (vertex_x > (sbnd_length_x - sbnd_offset_x)) 
           || (vertex_x < (-sbnd_offset_x))          
           || (vertex_y > (sbnd_length_y - sbnd_offset_y)) 
           || (vertex_y < (-sbnd_offset_y + sbnd_border_y))          
           || (vertex_z > (sbnd_length_z - sbnd_offset_z)) 
           || (vertex_z < (-sbnd_offset_z))
           || (end_x    > (sbnd_length_x - sbnd_offset_x)) 
           || (end_x    < (-sbnd_offset_x))          
           || (end_y    > (sbnd_length_y - sbnd_offset_y)) 
           || (end_y    < (-sbnd_offset_y))          
           || (end_z    > (sbnd_length_z - sbnd_offset_z)) 
           || (end_z    < (-sbnd_offset_z))); 

      bool does_vtx_escape = 
        (     (vertex_x > (sbnd_length_x - sbnd_offset_x)) 
           || (vertex_x < (-sbnd_offset_x))          
           || (vertex_y > (sbnd_length_y - sbnd_offset_y)) 
           || (vertex_y < (-sbnd_offset_y + sbnd_border_y))          
           || (vertex_z > (sbnd_length_z - sbnd_offset_z)) 
           || (vertex_z < (-sbnd_offset_z)));

      bool does_end_escape = 
        (     (end_x    > (sbnd_length_x - sbnd_offset_x)) 
           || (end_x    < (-sbnd_offset_x))          
           || (end_y    > (sbnd_length_y - sbnd_offset_y)) 
           || (end_y    < (-sbnd_offset_y))          
           || (end_z    > (sbnd_length_z - sbnd_offset_z)) 
           || (end_z    < (-sbnd_offset_z))); 

      bool one_end_escapes = true;
      if(does_vtx_escape && does_end_escape) one_end_escapes   = false;
      if(!does_vtx_escape && !does_end_escape) one_end_escapes = false;

      n_dedx = b_size->GetLeaf("tr_dedx_size")->GetValue(); // Get the number of entries for the dedx & residual range branches
      std::vector<float> dedx;
      std::vector<float> residual_range;
      dedx.clear();
      residual_range.clear();
      for(unsigned int j = 0; j < n_dedx; ++j){
        b_dedx->GetEntry(j);
        b_residual_range->GetEntry(j);
        dedx.push_back(b_dedx->GetLeaf("tr_dedx")->GetValue());
        residual_range.push_back(b_residual_range->GetLeaf("tr_residual_range")->GetValue());
      }

      track_list.push_back(Track(id_charge, id_energy, id_hits, n_hits, pida, chi2_mu, chi2_pi, chi2_pr, chi2_ka, length, kinetic_energy, mcs_mom_muon, range_mom_muon, range_mom_proton,vertex, end, !not_contained, one_end_escapes, dedx, residual_range));
    
    } 
  }
  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionHelper::GetShowerList(unsigned int start, TTree *shower_tree, const std::pair<int, int> &unique_event, ShowerList &shower_list){
  
    TBranch *b_event_id   = shower_tree->GetBranch("event_id");
    TBranch *b_time_now   = shower_tree->GetBranch("time_now");
    TBranch *b_n_hits     = shower_tree->GetBranch("sh_n_hits");
    TBranch *b_vertex     = shower_tree->GetBranch("sh_start");
    TBranch *b_direction  = shower_tree->GetBranch("sh_direction");
    TBranch *b_open_angle = shower_tree->GetBranch("sh_open_angle");
    TBranch *b_length     = shower_tree->GetBranch("sh_length");
    TBranch *b_energy     = shower_tree->GetBranch("sh_energy");
    
    unsigned int n_entries = shower_tree->GetEntries();

    for(unsigned int i = 0; i < n_entries; ++i){
    
      shower_tree->GetEntry(i);

      int event_id      = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now      = b_time_now->GetLeaf("time_now")->GetValue();
      
      if(event_id != unique_event.first || time_now != unique_event.second) continue;
      
      double temp_vertex[3];
      double temp_direction[3];
      
      temp_vertex[0]    = b_vertex->GetLeaf("sh_start")->GetValue(0);
      temp_vertex[1]    = b_vertex->GetLeaf("sh_start")->GetValue(1);
      temp_vertex[2]    = b_vertex->GetLeaf("sh_start")->GetValue(2);
      temp_direction[0] = b_direction->GetLeaf("sh_direction")->GetValue(0);
      temp_direction[1] = b_direction->GetLeaf("sh_direction")->GetValue(1);
      temp_direction[2] = b_direction->GetLeaf("sh_direction")->GetValue(2);
      float open_angle  = b_open_angle->GetLeaf("sh_open_angle")->GetValue();
      float length      = b_length->GetLeaf("sh_length")->GetValue();
      float energy      = b_energy->GetLeaf("sh_energy")->GetValue();
      int n_hits        = b_n_hits->GetLeaf("sh_n_hits")->GetValue();
 
      TVector3 vertex(temp_vertex);
      TVector3 direction(temp_direction);

      shower_list.push_back(Shower(n_hits, vertex, direction, open_angle, length, energy));
    } 
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionHelper::GetMCParticleList(unsigned int start, TTree *mcparticle_tree, const std::pair<int, int> &unique_event, ParticleList &mcparticle_list){
    
    TBranch *b_event_id = mcparticle_tree->GetBranch("event_id");
    TBranch *b_time_now = mcparticle_tree->GetBranch("time_now");
    TBranch *b_id       = mcparticle_tree->GetBranch("p_id");
    TBranch *b_n_hits   = mcparticle_tree->GetBranch("p_n_hits");
    TBranch *b_pdgcode  = mcparticle_tree->GetBranch("p_pdgcode");
    TBranch *b_status   = mcparticle_tree->GetBranch("p_status");
    TBranch *b_mass     = mcparticle_tree->GetBranch("p_mass");
    TBranch *b_energy   = mcparticle_tree->GetBranch("p_energy");
    TBranch *b_vertex   = mcparticle_tree->GetBranch("p_vertex");
    TBranch *b_end      = mcparticle_tree->GetBranch("p_end");
    TBranch *b_momentum = mcparticle_tree->GetBranch("p_momentum");
    
    unsigned int n_entries = mcparticle_tree->GetEntries();

      for(unsigned int i = 0; i < n_entries; ++i){
    
      mcparticle_tree->GetEntry(i);
      
      int event_id          = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now          = b_time_now->GetLeaf("time_now")->GetValue();
      
      if(event_id != unique_event.first || time_now != unique_event.second) continue;

      double temp_vertex[3];
      double temp_end[3];
      double temp_momentum[3];
      
      int id                = b_id->GetLeaf("p_id")->GetValue();
      int pdgcode           = b_pdgcode->GetLeaf("p_pdgcode")->GetValue();
      int statuscode        = b_status->GetLeaf("p_status")->GetValue();
      int n_hits            = b_n_hits->GetLeaf("p_n_hits")->GetValue();
      float mass            = b_mass->GetLeaf("p_mass")->GetValue();
      float energy          = b_energy->GetLeaf("p_energy")->GetValue();
      temp_vertex[0]        = b_vertex->GetLeaf("p_vertex")->GetValue(0);
      temp_vertex[1]        = b_vertex->GetLeaf("p_vertex")->GetValue(1);
      temp_vertex[2]        = b_vertex->GetLeaf("p_vertex")->GetValue(2);
      temp_end[0]           = b_end->GetLeaf("p_end")->GetValue(0);
      temp_end[1]           = b_end->GetLeaf("p_end")->GetValue(1);
      temp_end[2]           = b_end->GetLeaf("p_end")->GetValue(2);
      temp_momentum[0]      = b_momentum->GetLeaf("p_momentum")->GetValue(0);
      temp_momentum[1]      = b_momentum->GetLeaf("p_momentum")->GetValue(1);
      temp_momentum[2]      = b_momentum->GetLeaf("p_momentum")->GetValue(2);
 
      TVector3 vertex(temp_vertex);
      TVector3 end(temp_end);
      TVector3 momentum(temp_momentum);

      mcparticle_list.push_back(Particle(id, pdgcode, statuscode, n_hits, mass, energy, vertex, end, momentum));
      
    }
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionHelper::GetRecoParticleFromTrack1EscapingDistanceCut(const TrackList &track_list, ParticleList &recoparticle_list){

    // Assign ridiculously short length to initiate the longest track length
    float longest_track_length      = -std::numeric_limits<float>::max();
    unsigned int longest_track_id   =  std::numeric_limits<unsigned int>::max();
   
    // Check if exactly 1 track escapes
    bool exactly_one_escapes = false;
    unsigned int n_escaping = 0;
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &trk(track_list[i]);
      if(trk.m_one_end_contained) n_escaping++;
    }
    if(n_escaping == 1) exactly_one_escapes = true;

    // Check if the track is contained and passes the distance cut to the escaping border
    bool contained_and_passes_distance_cut = false;
    if(exactly_one_escapes){
      for(unsigned int i = 0; i < track_list.size(); ++i){
        const Track &trk(track_list[i]);
        if(trk.m_one_end_contained){
          // Find out if the neutrino vertex is far enough from the escaping face
          float distance_to_intersection_point = -std::numeric_limits<float>::max();
          // Loop over the fiducial planes and find out which the escaping particle passed through
          PlaneList planes;
          EventSelectionHelper::GetSBNDFiducialPlanes(planes);
          for(const Plane &plane : planes){
            if(!EventSelectionHelper::CheckIfTrackIntersectsPlane(plane, trk)) continue;
            distance_to_intersection_point = EventSelectionHelper::GetDistanceFromTrackToPlane(plane,trk);
            if(distance_to_intersection_point > 50){
              contained_and_passes_distance_cut = true;
              break;
            }
          }
        }
      }
    }
    
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &candidate(track_list[i]);
      // Get the lengths of the tracks and find the longest track and compare to the rest of
      // the lengths
      if(candidate.m_length > longest_track_length) {
        longest_track_length = candidate.m_length;
        longest_track_id     = i;
      }
    }

    bool always_longest(true);
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      // Find out if the longest track is always 1.5x longer than all the others in the event
      const Track &track(track_list[id]);
      if(track.m_length*1.5 >= longest_track_length && id != longest_track_id) always_longest = false;
    }
  
    // Muon candidates 
    std::vector<unsigned int> mu_candidates;
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(contained_and_passes_distance_cut){
        // If one end is contained and the neutrino vertex is more than 75 cm from the escaping border
        if(track.m_one_end_contained) 
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 13, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range));
        else if(EventSelectionHelper::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range));
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range)); 
      }
      else{
        // If the Chi2 Proton hypothesis gives proton, call the track a proton
        // Otherwise, call it a muon candidate
        if(EventSelectionHelper::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range));
        else if(EventSelectionHelper::GetPdgByChi2MuonCandidate(track) == 13 || (id == longest_track_id && always_longest))
          mu_candidates.push_back(id);
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range));
      }
    }

    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi, muon.m_dedx, muon.m_residual_range));
      return;
    }
    
    // If more than one muon candidate exists
    bool foundTheMuon(false);
    unsigned int muonID = std::numeric_limits<unsigned int>::max();

    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      const Track &candidate(track_list[id]);
      if(longest_track_id == id && always_longest) {
        muonID = id;
        foundTheMuon = true;
        break;
      }
    }
    if(!foundTheMuon) {
      // Find the smallest chi^2 under the muon hypothesis
      muonID = EventSelectionHelper::GetMuonByChi2(track_list, mu_candidates);
      if(muonID != std::numeric_limits<unsigned int>::max()) foundTheMuon = true;
      else throw 10;
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi, muon.m_dedx, muon.m_residual_range));
    
    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range)); 
    } 
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionHelper::GetRecoParticleFromShower(const ShowerList &shower_list, const TVector3 &reco_vertex, ParticleList &recoparticle_list){

    // New method
    // Calculate distance of closest approach between 2 photon showers 
    // There must be at least 2 showers to try and find a pi0

    std::vector<unsigned int> used_photon;
    ShowerList unused_showers;
    std::vector<unsigned int>::iterator it;
    unused_showers.clear();
    used_photon.clear();

    if(shower_list.size() == 1) unused_showers.push_back(shower_list[0]);
    else{
      // Loop over showers
      for(unsigned int i = 0; i < shower_list.size(); ++i){

        // Vector to hold 'j' values of candidate photon to pair with current photon
        // Vector to hold to position at which the photons meet to call it the point at which
        // the pi0 decayed
        // Vector to hold the distance of closest approach, in order to minimise this
        std::vector<unsigned int> candidate_id_for_pair;
        std::vector<TVector3> decay_point;
        std::vector<float> c_distance;
        std::vector<float> pi0_energy;
        std::vector<float> pi0_inv_mass_diff;
        std::vector<int> total_hits;
        candidate_id_for_pair.clear();
        decay_point.clear();
        c_distance.clear();
        total_hits.clear();

        for( unsigned int j = i+1; j < shower_list.size(); ++j){
          // If we are only looking at a single photon, continue
          if(i==j) continue;

          // If the photon has already been assigned to a pi0
          for(unsigned int k = 0; k < used_photon.size(); ++k) if(used_photon[k] == j) continue;

          // Get the distance of closest approach of the current two showers we are looking at
          TVector3 dir_1, dir_2, start_1, start_2, link;

          dir_1   = shower_list[i].m_direction;
          dir_2   = shower_list[j].m_direction;
          start_1 = shower_list[i].m_vertex;
          start_2 = shower_list[j].m_vertex;
          link    = start_1 - start_2;

          float a = dir_1.Dot(dir_1); 
          float b = dir_1.Dot(dir_2); 
          float c = dir_2.Dot(dir_2); 
          float d = dir_1.Dot(link); 
          float e = dir_2.Dot(link);
          float denomenator = a*c - b*b;

          // Get the invariant mass of the current 2 showers
          float energy_1   = shower_list[i].m_energy;
          float energy_2   = shower_list[j].m_energy;

          float cos_theta  = (b / (std::sqrt(a)*std::sqrt(c)));  
          float inv_mass   = std::sqrt(2*energy_1*energy_2*(1-cos_theta));
          float pi0_mass   = 134.97; // MeV

          float mass_diff  = std::abs(pi0_mass/1000. - inv_mass);

          // If the lines are parallel
          if(denomenator == 0) continue;

          float m_closest = (b*e - c*d)/denomenator;
          float n_closest = (a*e - b*d)/denomenator;

          TVector3 d_closest  = link + ((b*e - c*d)*dir_1 - (a*e - b*d)*dir_2)*(1/denomenator);
          float mag_d_closest =  d_closest.Mag();

          TVector3 d_middle   = (start_1 + m_closest*dir_1) + 0.5*d_closest;

          // If the distance of closest approach is smaller than 15 cm 
          // or the invariant mass of the photons is within 20% of the pion mass
          // call it a candidate
          if(mag_d_closest < 15 && mass_diff * (1./pi0_mass) < 0.2) { 
            candidate_id_for_pair.push_back(j);
            decay_point.push_back(d_closest);
            c_distance.push_back(mag_d_closest);
            total_hits.push_back(shower_list[i].m_n_hits + shower_list[j].m_n_hits);
            pi0_inv_mass_diff.push_back(mass_diff);
            pi0_energy.push_back(energy_1 + energy_2);
          }
        } // Inner shower list

        if(candidate_id_for_pair.size() == 0) continue;
        if(candidate_id_for_pair.size() == 1){
          // If the location with respect to the neutrino vertex at which the photons 
          // were produced by the candidate pi0 is more than 15 cm continue
          if((reco_vertex - decay_point[0]).Mag() > 15) continue;

          // push back a pi0 corresponding to the two photons
          recoparticle_list.push_back(Particle(111, total_hits[0], reco_vertex, decay_point[0], pi0_energy[0]));
          used_photon.push_back(i);
          used_photon.push_back(candidate_id_for_pair[0]);

        } 
        else{
          // Find the minimum distance
          std::vector<float>::iterator min      = std::min_element(c_distance.begin(), c_distance.end());
          TVector3 best_decay_point             = decay_point[std::distance(c_distance.begin(), min)];
          int best_total_hits                   = total_hits[std::distance(c_distance.begin(), min)];
          float best_pi0_energy                 = pi0_energy[std::distance(c_distance.begin(), min)];
          recoparticle_list.push_back(Particle(111,best_total_hits, reco_vertex, best_decay_point, best_pi0_energy));
          used_photon.push_back(candidate_id_for_pair[std::distance(c_distance.begin(), min)]);
          used_photon.push_back(i);
        } // candidates
      } // Shower loop
      // Find any unused showers
      if(shower_list.size() > used_photon.size() && used_photon.size() != 0){
        for(unsigned int i = 0; i < shower_list.size(); ++i) {
          // If the current shower is not in the used photon list
          // Push it to unused showers
          it = std::find(used_photon.begin(), used_photon.end(), i);
          if(it == used_photon.end())
            unused_showers.push_back(shower_list[i]);
        }
      } // Used photons
    } // If more than 1 shower

    // For every unused shower, check the distance of the shower vertex from the neutrino vertex
    // If it is <= 15 cm, call it an electron
    // If it is > 15 cm, call it a photon 
    for(unsigned int i = 0; i < unused_showers.size(); ++i) {
      TVector3 shower_start = unused_showers[i].m_vertex;
      double conversion_length = (shower_start - reco_vertex).Mag();
      if(conversion_length <= 15) recoparticle_list.push_back(Particle(11, unused_showers[i].m_n_hits, unused_showers[i].m_vertex, (unused_showers[i].m_length * unused_showers[i].m_direction ), unused_showers[i].m_energy));
      else recoparticle_list.push_back(Particle(22, unused_showers[i].m_n_hits, unused_showers[i].m_vertex, (unused_showers[i].m_length * unused_showers[i].m_direction ), unused_showers[i].m_energy));
    } // Unused showers
  }

  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionHelper::GetPdgByChi2(const Track &track){

    // Push the chi2 values onto a vector to find the minimum
    // NOT MUON
    std::map<int, float> chi2_map;

    chi2_map.insert(std::map<int, float>::value_type(211,  track.m_chi2_pi));
    chi2_map.insert(std::map<int, float>::value_type(321,  track.m_chi2_ka));
    chi2_map.insert(std::map<int, float>::value_type(2212, track.m_chi2_pr));

    float min_chi2 = std::numeric_limits<float>::max();
    int best_pdg   = std::numeric_limits<int>::max();

    for(std::map<int, float>::const_iterator it = chi2_map.begin(); it != chi2_map.end(); ++it){
    
      if(it->second < min_chi2){
        min_chi2 = it->second; 
        best_pdg = it->first;
      }
    } 
    return best_pdg;
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionHelper::GetMuonByChi2(const TrackList &tracks, const std::vector<unsigned int> &mu_candidates){

    // Loop over muon candidates and find smallest corresponding chi2_mu
    float min_chi2_mu = std::numeric_limits<float>::max();
    unsigned int min_chi2_id = std::numeric_limits<unsigned int>::max();

    // Find the smallest chi^2 under the muon hypothesis
    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      const Track &candidate(tracks[id]);
      if(candidate.m_chi2_mu < min_chi2_mu) {
        min_chi2_mu = candidate.m_chi2_mu; 
        min_chi2_id = id;
      }
    }
    return min_chi2_id;
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionHelper::GetProtonByChi2Proton(const Track &track){

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_pr < 80) return 2212;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionHelper::GetMuonByChi2Proton(const Track &track){

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_pr > 80) return 13;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionHelper::GetPdgByChi2MuonCandidate(const Track &track){

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_mu < 16) return 13;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionHelper::GetPdgByPIDA(const Track &track){

    // Muon
    if(track.m_pida >= 5 && track.m_pida < 9) return 13;

    //Pion
    if(track.m_pida >= 9 && track.m_pida < 13) return 211;
  
    //Kaon
    if(track.m_pida >= 13 && track.m_pida < 13.5) return 321;
  
    //Proton
    if(track.m_pida >= 13) return 2212;

    return std::numeric_limits<int>::max();
  
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionHelper::GetPdgByPIDAStrict(const Track &track){

    // Muon
    if(track.m_pida >= 7.5  && track.m_pida < 8) return 13;

    //Pion
    if(track.m_pida >= 8    && track.m_pida < 9) return 211;
  
    //Kaon
    if(track.m_pida >= 12.9 && track.m_pida < 13.5) return 321;
  
    //Proton
    if(track.m_pida >= 16.9 && track.m_pida < 17.4) return 2212;

    return std::numeric_limits<int>::max();
  
  }
  
  //------------------------------------------------------------------------------------------ 
  
  EventSelectionHelper::Track::Track(const int mc_id_charge, const int mc_id_energy, const int mc_id_hits, const int n_hits, const float pida, const float chi2_mu, const float chi2_pi, const float chi2_pr, const float chi2_ka, const float length, const float kinetic_energy, const float mcs_momentum_muon, const float range_momentum_muon, const float range_momentum_proton, const TVector3 &vertex, const TVector3 &end, const bool &contained, const bool &one_end_contained, const std::vector<float> &dedx, const std::vector<float> &residual_range) :
    m_mc_id_charge(mc_id_charge),
    m_mc_id_energy(mc_id_energy),
    m_mc_id_hits(mc_id_hits),
    m_n_hits(n_hits),
    m_pida(pida),
    m_chi2_mu(chi2_mu),
    m_chi2_pi(chi2_pi),
    m_chi2_pr(chi2_pr),
    m_chi2_ka(chi2_ka),
    m_length(length),
    m_kinetic_energy(kinetic_energy),
    m_mcs_mom_muon(mcs_momentum_muon),
    m_range_mom_muon(range_momentum_muon),
    m_range_mom_proton(range_momentum_proton),
    m_vertex(vertex),
    m_end(end),
    m_contained(contained),
    m_one_end_contained(one_end_contained),
    m_dedx(dedx),
    m_residual_range(residual_range) {}

  //------------------------------------------------------------------------------------------ 
  
  EventSelectionHelper::Shower::Shower(const int n_hits, const TVector3 &vertex, const TVector3 &direction, const float open_angle, const float length, const float energy) :
    m_n_hits(n_hits),
    m_vertex(vertex),
    m_direction(direction),
    m_open_angle(open_angle),
    m_length(length),
    m_energy(energy/1000.) {}

} // namespace: selection
/*
  void EventSelectionHelper::GetRecoParticleFromTrackRaquel(const TrackList &track_list, ParticleList &recoparticle_list){

    // Assign ridiculously short length to initiate the longest track length
    float longest_track_length      = -std::numeric_limits<float>::max();
    unsigned int longest_track_id   =  std::numeric_limits<unsigned int>::max();
  
    // Get the longest track id
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &candidate(track_list[i]);
      // Get the lengths of the tracks and find the longest track and compare to the rest of
      // the lengths
      if(candidate.m_length > longest_track_length) {
        longest_track_length = candidate.m_length;
        longest_track_id     = i;
      }
    }
    bool always_longest(true);
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      // Find out if the longest track is always 1.5x longer than all the others in the event
      const Track &track(track_list[id]);
      if(track.m_length*1.5 >= longest_track_length && id != longest_track_id) always_longest = false;
    }
   
    // Muon candidates 
    std::vector<unsigned int> mu_candidates;

    // Check if exactly 1 track escapes and that it is the longest track
    bool only_longest_escapes = false;
    bool none_escape = true;
    bool longest_escapes = false;
    unsigned int n_escaping = 0;
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &trk(track_list[i]);
      if(i == longest_track_id) longest_escapes = true;
      if(trk.m_one_end_contained) n_escaping++;
    }
    if(n_escaping == 1 && longest_escapes) only_longest_escapes = true;
    if(n_escaping != 0) none_escape = false;

    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(only_longest_escapes){
        if(track.m_one_end_contained)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 13, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton, track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        else if(EventSelectionHelper::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton, track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton, track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi)); 
      }
      else if(!none_escape) continue;
      else{
        // If the Chi2 Proton hypothesis gives proton, call the track a proton
        // Otherwise, call it a muon candidate
        if(EventSelectionHelper::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        else if(EventSelectionHelper::GetMuonByChi2Proton(track) == 13 || (id == longest_track_id && always_longest) || track.m_one_end_contained)
          mu_candidates.push_back(id);
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
      }
    }

    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi));
      return;
    }
    
    // If more than one muon candidate exists
    bool foundTheMuon(false);
    unsigned int muonID = std::numeric_limits<unsigned int>::max();

    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      const Track &candidate(track_list[id]);
      if(longest_track_id == id && always_longest) {
        muonID = id;
        foundTheMuon = true;
        break;
      }
    }
    if(!foundTheMuon) {
      // Find the smallest chi^2 under the muon hypothesis
      muonID = EventSelectionHelper::GetMuonByChi2(track_list, mu_candidates);
      if(muonID != std::numeric_limits<unsigned int>::max()) foundTheMuon = true;
      else throw 10;
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi));
    
    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi)); 
    } 
  }
*/
  //------------------------------------------------------------------------------------------ 
/*
  void EventSelectionHelper::GetRecoParticleFromTrackChi2P(const TrackList &track_list, ParticleList &recoparticle_list){

    // Assign ridiculously short length to initiate the longest track length
    float longest_track_length      = -std::numeric_limits<float>::max();
    unsigned int longest_track_id   =  std::numeric_limits<unsigned int>::max();
   
    // Check if exactly 1 track escapes
    bool exactly_one_escapes = false;
    unsigned int n_escaping = 0;
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &trk(track_list[i]);
      if(trk.m_one_end_contained) n_escaping++;
    }
    if(n_escaping == 1) exactly_one_escapes = true;

    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &candidate(track_list[i]);
      // Get the lengths of the tracks and find the longest track and compare to the rest of
      // the lengths
      if(candidate.m_length > longest_track_length) {
        longest_track_length = candidate.m_length;
        longest_track_id     = i;
      }
    }

    bool always_longest(true);
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      // Find out if the longest track is always 1.5x longer than all the others in the event
      const Track &track(track_list[id]);
      if(track.m_length*1.5 >= longest_track_length && id != longest_track_id) always_longest = false;
    }
   
    // Muon candidates 
    std::vector<unsigned int> mu_candidates;
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(exactly_one_escapes){
        if(track.m_one_end_contained) 
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 13, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        else if(EventSelectionHelper::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi)); 
      }
      else{
        // If the Chi2 Proton hypothesis gives proton, call the track a proton
        // Otherwise, call it a muon candidate
        if(EventSelectionHelper::GetProtonByChi2Proton(track) == 2212){
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        }
        else if(EventSelectionHelper::GetMuonByChi2Proton(track) == 13 || (id == longest_track_id && always_longest) || track.m_one_end_contained)
          mu_candidates.push_back(id);
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
      }
    }

    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi));
      return;
    }
    
    // If more than one muon candidate exists
    bool foundTheMuon(false);
    unsigned int muonID = std::numeric_limits<unsigned int>::max();

    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      const Track &candidate(track_list[id]);
      if(longest_track_id == id && always_longest) {
        muonID = id;
        foundTheMuon = true;
        break;
      }
    }
    if(!foundTheMuon) {
      // Find the smallest chi^2 under the muon hypothesis
      muonID = EventSelectionHelper::GetMuonByChi2(track_list, mu_candidates);
      if(muonID != std::numeric_limits<unsigned int>::max()) foundTheMuon = true;
      else throw 10;
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi));
    
    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi)); 
    } 
  }
*/
  //------------------------------------------------------------------------------------------ 

/*  void EventSelectionHelper::GetRecoParticleFromTrack1Escaping(const TrackList &track_list, ParticleList &recoparticle_list){

    // Assign ridiculously short length to initiate the longest track length
    float longest_track_length      = -std::numeric_limits<float>::max();
    unsigned int longest_track_id   =  std::numeric_limits<unsigned int>::max();
   
    // Check if exactly 1 track escapes
    bool exactly_one_escapes = false;
    unsigned int n_escaping = 0;
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &trk(track_list[i]);
      if(trk.m_one_end_contained) n_escaping++;
    }
    if(n_escaping == 1) exactly_one_escapes = true;

    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &candidate(track_list[i]);
      // Get the lengths of the tracks and find the longest track and compare to the rest of
      // the lengths
      if(candidate.m_length > longest_track_length) {
        longest_track_length = candidate.m_length;
        longest_track_id     = i;
      }
    }

    bool always_longest(true);
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      // Find out if the longest track is always 1.5x longer than all the others in the event
      const Track &track(track_list[id]);
      if(track.m_length*1.5 >= longest_track_length && id != longest_track_id) always_longest = false;
    }
   
    // Muon candidates 
    std::vector<unsigned int> mu_candidates;
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(exactly_one_escapes){
        if(track.m_one_end_contained) 
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 13, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        else if(EventSelectionHelper::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi)); 
      }
      else{
        // If the Chi2 Proton hypothesis gives proton, call the track a proton
        // Otherwise, call it a muon candidate
        if(EventSelectionHelper::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
        else if(EventSelectionHelper::GetPdgByChi2MuonCandidate(track) == 13 || (id == longest_track_id && always_longest))
          mu_candidates.push_back(id);
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi));
      }
    }

    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi));
      return;
    }
    
    // If more than one muon candidate exists
    bool foundTheMuon(false);
    unsigned int muonID = std::numeric_limits<unsigned int>::max();

    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      const Track &candidate(track_list[id]);
      if(longest_track_id == id && always_longest) {
        muonID = id;
        foundTheMuon = true;
        break;
      }
    }
    if(!foundTheMuon) {
      // Find the smallest chi^2 under the muon hypothesis
      muonID = EventSelectionHelper::GetMuonByChi2(track_list, mu_candidates);
      if(muonID != std::numeric_limits<unsigned int>::max()) foundTheMuon = true;
      else throw 10;
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi));
    
    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionHelper::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi)); 
    } 
  }
*/
  //------------------------------------------------------------------------------------------ 
