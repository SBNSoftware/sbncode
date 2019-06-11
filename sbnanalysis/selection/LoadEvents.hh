#ifndef EVENT_SELECTION_TOOL_H
#define EVENT_SELECTION_TOOL_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "EventSelectionHelper.hh"
#include "../core/Event.hh"

namespace selection{
 
  /**
   * @brief  LoadEvents helper class
   */
  class LoadEvents {

    private : 
      class Track;
      class Shower;

    public : 

      typedef std::vector<Plane>               PlaneList;
      typedef std::vector<Particle>            ParticleList;
      typedef std::vector<Event>               EventList;
      typedef std::vector<Track>               TrackList;
      typedef std::vector<Shower>              ShowerList;
      
      /**
       * @brief get the pot corresponding to each individual file
       *
       * @param  subrun tree with the some of the subrun information from larsoft
       *
       * @return pot from the current file
       *
       */
      static double GetPOT(TTree *subrun);

      /**
       * @brief  load the list of events to analyse from the root file
       *
       * @param  file_name name of the root file to access
       * @param  event_list vector of events to fill
       * @param  file number of the file of the current event
       * @param  pot from the current file
       *
       */
      static void LoadEventList(const std::string &file_name, EventList &event_list, const int &file, double &pot);

    private :

      /**
       * @brief  get the list of track objects
       *
       * @param  track_tree tree to take track information from
       * @param  unique_event_list list of unique events to take track information from
       * @param  track_list vector of tracks to fill
       *
       */
      static void GetTrackList(unsigned int start, TTree *track_tree, const std::pair<int, int> &unique_event, TrackList &track_list);

      /**
       * @brief  get the list of shower objects
       *
       * @param  shower_tree tree to take shower information from
       * @param  unique_event_list list of unique events to take shower information from
       * @param  shower_list vector of showers to fill
       *
       */
      static void GetShowerList(unsigned int start, TTree *shower_tree, const std::pair<int, int> &unique_event, ShowerList &shower_list);
      
      /**
       * @brief  get the list of mc particle objects
       *
       * @param  mc particle_tree tree to take mc particle information from
       * @param  unique_event_list list of unique events to take mc particle information from
       * @param  mc particle_list vector of mc particles to fill
       *
       */
      static void GetMCParticleList(unsigned int start, TTree *mcparticle_tree, const std::pair<int, int> &unique_event, ParticleList &mcparticle_list);

      /**
       * @brief  get a list of reconstructed particles from track objects
       *
       * @param  track_list list of tracks in the event
       * @param  recoparticle_list particle list to fill
       *
       */
      static void GetRecoParticleFromTrack1EscapingDistanceCut(const TrackList &track_list, ParticleList &recoparticle_list);

      /**
       * @brief  get a list of reconstructed particles from track objects using original method
       *
       * @param  track_list list of tracks in the event
       * @param  recoparticle_list particle list to fill
       *
       */
      static void GetRecoParticleFromShower(const ShowerList &shower_list, const TVector3 &reco_vertex, ParticleList &recoparticle_list);

      /**
       * @brief  Track class 
       */
      class Track{
      
        public : 
          
          /**
           * @brief  Constructor
           *
           * @param  mc_id_charge mc TrackID corresponding to MCParticle using charge 
           * @param  mc_id_energy mc TrackID corresponding to MCParticle using energy
           * @param  mc_id_hits mc TrackID corresponding to MCParticle using hits
           * @param  pida pida value
           * @param  chi2_mu chi squared value for the muon fit of the reconstructed dEdx to the expected distribution
           * @param  chi2_pi chi squared value for the pion fit of the reconstructed dEdx to the expected distribution
           * @param  chi2_pr chi squared value for the proton fit of the reconstructed dEdx to the expected distribution
           * @param  chi2_ka chi squared value for the kaon fit of the reconstructed dEdx to the expected distribution
           * @param  length track length
           * @param  kinetic_energy track kinetic energy
           * @param  vertex vertex of the track
           * @param  end end point of the track
           * @param  contained whether or not the reconstructed track is contained within the SBND fiducial volume
           * @param  one_end_contained whether or not the reconstructed track has one end contained within the SBND fiducial volume
           *
           */
          Track(const int mc_id_charge, const int mc_id_energy, const int mc_id_hits, const int n_hits, const float pida, const float chi2_mu, const float chi2_pi, const float chi2_pr, const float chi2_ka, const float length, const float kinetic_energy, const float mcs_momentum_muon, const float range_momentum_muon, const float range_momentum_proton, const TVector3 &vertex, const TVector3 &end, const bool &contained, const bool &one_end_contained, const std::vector<float> &dedx, const std::vector<float> &residual_range);

          // Member variables
          int      m_mc_id_charge;             ///< mc TrackID corresponding to MCParticle using charge
          int      m_mc_id_energy;             ///< mc TrackID corresponding to MCParticle using energy
          int      m_mc_id_hits;               ///< mc TrackID corresponding to MCParticle using hits
          int      m_n_hits;                   ///< number of hits in the track
          float    m_pida;                     ///< pida value
          float    m_chi2_mu;                  ///< chi squared fit to the muon expected dEdx
          float    m_chi2_pi;                  ///< chi squared fit to the pion expected dEdx
          float    m_chi2_pr;                  ///< chi squared fit to the proton expected dEdx 
          float    m_chi2_ka;                  ///< chi squared fit to the kaon expected dEdx
          float    m_length;                   ///< length of the track
          float    m_kinetic_energy;           ///< kinetic energy of the track
          float    m_mcs_mom_muon;             ///< multiple coulomb scattering momentum is the particle is an escaping muon
          float    m_range_mom_muon;           ///< range momentum if the particle is a contained muon 
          float    m_range_mom_proton;         ///< range momentum if the particle is a contained proton
          TVector3 m_vertex;                   ///< vertex of the track         
          TVector3 m_end;                      ///< end of the track
          bool     m_contained;                ///< whether or not the reconstructed track is contained
          bool     m_one_end_contained;        ///< whether or not the reconstructed track has one contained end
          std::vector<float> m_dedx;           ///< vector of the dedx distribution of the reconstructed track
          std::vector<float> m_residual_range; ///< vector of the residual range distribution of the track
      
      }; // Track
      
      /**
       * @brief  Shower class 
       */
      class Shower{
      
        public : 
          
          /**
           * @brief  Constructor
           *
           * @param  vertex vertex of the shower
           * @param  direction direction of the shower
           * @param  open_angle opening angle at the vertex of the shower
           * @param  length length of the shower
           *
           */
          Shower(const int n_hits, const TVector3 &vertex, const TVector3 &direction, const float open_angle, const float length, const float energy);

          // Member variables
          int      m_n_hits;     ///< number of hits in the shower
          TVector3 m_vertex;     ///< vertex of the shower 
          TVector3 m_direction;  ///< direction of the shower
          float    m_open_angle; ///< opening angle at the vertex of the shower
          float    m_length;     ///< length of the shower
          float    m_energy;     ///< energy of the shower

      }; // Shower
  }; // LoadEvents
} // namespace: selection
#endif
