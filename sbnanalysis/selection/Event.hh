#ifndef EVENT_H
#define EVENT_H

#include "Plane.hh"
#include "Particle.hh"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include "TH1.h"
#include "TF1.h"

namespace selection{
  
  /*
   * @brief typedef for a topology defintion
   */
  typedef std::map< std::vector< int >, int > TopologyMap;
  /*
   * @brief typedef for a list of selection::Particles
   */
  typedef std::vector<Particle> ParticleList;
  /*
   * @brief typedef for a list of selection::Particles
   */
  typedef std::vector<Plane> PlaneList;

  /**
   * @brief  Event class
   */
  class Event{

    public :
    
      /**
       * @brief  Default Constructor
       */
      Event(); 

      /**
       * @brief  Constructor
       *
       * @param  mc_particles list of the MC particle objects in the event
       * @param  reco_particles list of the reconstructed particle objects in the event
       * @param  interaction the interaction type corresponding to the event
       * @param  scatter the genie scattering code of the event
       * @param  neutrino_pdg final neutrino pdgcode
       * @param  initial_pdg initial neutrino pdgcode 
       * @param  is_cc is this a charged or neutral current event
       * @param  mc_vertex Monte Carlo neutrino vertex 
       * @param  reco_vertex reconstructed neutrino vertex
       * @param  neutrino_energy energy of the neutrino
       * @param  neutrino_qsqr qsqr of the neutrino
       */
      Event(const ParticleList &mc_particles, 
            const ParticleList &reco_particles, 
            const unsigned int interaction, 
            const unsigned int scatter, 
            const int neutrino_pdg, 
            const int initial_pdg, 
            const bool is_cc, 
            const TVector3 &mc_vertex, 
            const TVector3 &reco_vertex, 
            const float neutrino_energy, 
            const float neutrino_qsqr); 

      /**
       * @brief  CountMCParticlesWithPdg
       *
       * @param  pdg pdg code to count
       *
       * @return number of Monte Carlo particles with the given pdg code
       */
      unsigned int CountMCParticlesWithPdg(const int pdg) const;

      /**
       * @brief  CountRecoParticlesWithPdg
       *
       * @param  pdg pdg code to count
       *
       * @return number of reconstructed partices with the given pdg code
       */
      unsigned int CountRecoParticlesWithPdg(const int pdg) const;

      /**
       * @brief  CheckMCTopology
       *
       * @param  topology map of the chosen topology which holds the number of each chosen type of particle to look for
       *
       * @return boolean as to whether the event has the desired Monte Carlo topology
       */
      bool CheckMCTopology(const TopologyMap &topology) const;

      /**
       * @brief  CheckRecoTopology
       *
       * @param  topology map of the chosen topology which holds the number of each chosen type of particle to look for
       *
       * @return boolean as to whether the event has the desired reconstructed topology
       */
      bool CheckRecoTopology(const TopologyMap &topology) const;

      /**
       * @brief  Get the list of MC particls for this event
       */
      ParticleList GetMCParticleList() const;

      /**
       * @brief  Get the list of reconstructed particls for this event
       */
      ParticleList GetRecoParticleList() const;
      
      /**
       * @brief  Get the file id the current event came from
       */
      int GetFileId() const;

      /**
       * @brief  Get the id of the current event
       */
      int GetId() const;

      /**
       * @brief  Get the interaction type of the event \n
       * <tt>
       * 0 : Unknown \n
       * 1 : Weak CC \n
       * 2 : Weak NC \n
       * 3 : Weak CC + NC + Interference \n
       * 4 : Nucleon decay \n
       * </tt>
       */
      int GetInteractionType() const;

      /**
       * @brief  Get the scattering code of the event: the physical process \n
       *  <tt>
       * 0  : Unknown \n
       * 1  : QE \n
       * 2  : Single kaon \n
       * 3  : DIS \n
       * 4  : Resonant \n
       * 5  : Coherent \n
       * 6  : Diffractive \n
       * 7  : \f$ \nu \f$- e elastic \n
       * 8  : Inverse \f$ \mu \f$ decay \n
       * 9  : AM \f$ \nu - \gamma \f$ \n
       * 10 : MEC \n
       * 11 : Coherent elastic \n
       * 12 : Inverse \f$ \beta \f$ decay \n
       * 13 : Glashow resonance \n
       * 14 : IMD Annihilation \n
       * </tt>
       *
       */
      int GetPhysicalProcess() const;
      
      /**
       * @brief  Get the neutrino pdg code in the event
       */
      int GetNeutrinoPdgCode() const;

      /**
       * @brief  Get the initial neutrino pdg code in the event to see if it oscillated
       */
      int GetInitialNeutrinoPdgCode() const;

      /**
       * @brief  Get whether the true neutrino interaction happened within the SBND fiducial 
       *         volume
       */
      bool IsSBNDTrueFiducial() const;
      
      /**
       * @brief  Get whether all the reconstructed tracks in an event are contained
       */
      bool AllRecoContained() const;
     
      /**
       * @brief  Get if the event is CC or NC
       */
      bool GetIsCC() const;

      /**
       * @brief  Get the Monte Carlo neutrino vertex position
       */
      TVector3 GetMCNuVertex() const;

      /**
       * @brief  Get the reconstructed neutrino vertex position
       */
      TVector3 GetRecoNuVertex() const;
      
      /**
       * @brief  Get the true neutrino energy
       */
      float GetTrueNuEnergy() const;
      
      /**
       * @brief  Get the true neutrino qsqr
       */
      float GetTrueNuQ2() const;
      
      /**
       * @brief  Get the most energetic reconstructed particle
       *
       * @return Particle most energetic reco
       */
      Particle GetMostEnergeticRecoParticle() const;

      /**
       * @brief  Get the most energetic true particle
       *
       * @return Particle most energetic true
       */
      Particle GetMostEnergeticTrueParticle() const;

      /**
       * @brief  Get the minimum x,y,z positions of the SBND fiducial volume
       *
       * @return Vector of lowest x,y,z positions
       */
      TVector3 GetMinimumFiducialDimensions() const;

      /**
       * @brief  Get the maximum x,y,z positions of the SBND fiducial volume
       *
       * @return Vector of highest x,y,z positions
       */
      TVector3 GetMaximumFiducialDimensions() const;

      /**
       * @brief  Get the number of escaping reconstructed particles
       */
      int NumberOfEscapingRecoParticles() const;
      
      /**
       * @brief  Get the number of true escaping particles
       */
      int NumberOfEscapingMCParticles() const;
      
    private : 

      /*
       * @brief  Get the number of escaping particles in the given list
       */
      int NumberOfEscapingParticles(const ParticleList &particles) const;

      /**
       * @brief  CountParticlesWithPdg
       *
       * @param  pdg pdg code to count
       *
       * @return number of partices with the given pdg code
       */
      unsigned int CountParticlesWithPdg(const int pdg, const ParticleList &particle_list) const;
      
      /**
       * @brief  CheckTopology
       *
       * @param  topology map of the chosen topology which holds the number of each chosen type of particle to look for
       *
       * @return boolean as to whether the event has the desired topology
       */
      bool CheckTopology(const TopologyMap &topology, const ParticleList &particle_list) const;

      /**
       * @brief  Get the most energetic particle
       *
       * @return Particle most energetic
       */
      Particle GetMostEnergeticParticle(const ParticleList &particle_list) const;

      // Member variables
      ParticleList       m_mc_particles;       ///< vector of Monte Carlo particles
      ParticleList       m_reco_particles;     ///< vector of reconstructed particles
      unsigned int       m_interaction;        ///< interaction type of the event
      unsigned int       m_scatter;            ///< scatter code for the event: physical process
      int                m_nu_pdg;             ///< Neutrino pdg code of the event
      int                m_init_pdg;           ///< Initial neutrino pdg code of the event
      unsigned int       m_charged_pi;         ///< Number of charged pions in the event
      unsigned int       m_neutral_pi;         ///< Number of neutral pions in the event
      bool               m_is_cc;              ///< whether the event contains and CC or NC interaction
      TVector3           m_reco_vertex;        ///< reconstructed neutrino vertex
      TVector3           m_mc_vertex;          ///< reconstructed neutrino vertex
      float              m_neutrino_energy;    ///< true neutrino energy
      float              m_neutrino_qsqr;      ///< true neutrino qsqr
      int                m_file;               ///< file id
      int                m_id;                 ///< event id
      float              m_sbnd_border_x;      ///< fiducial border in x for the sbnd detector
      float              m_sbnd_border_y;      ///< fiducial border in y for the sbnd detector
      float              m_sbnd_border_z;      ///< fiducial border in z for the sbnd detector
      float              m_sbnd_offset_x;      ///< offset in x for the sbnd detector
      float              m_sbnd_offset_y;      ///< offset in y for the sbnd detector
      float              m_sbnd_offset_z;      ///< offset in z for the sbnd detector
      float              m_sbnd_half_length_x; ///< detector half length in x
      float              m_sbnd_half_length_y; ///< detector half length in y
      float              m_sbnd_half_length_z; ///< detector half length in z

      float              m_unfilled = -99999.;

  }; // Event
} // selection

#endif
