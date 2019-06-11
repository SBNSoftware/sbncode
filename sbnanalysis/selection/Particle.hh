#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <iostream>
#include <vector>
#include "TVector3.h"

namespace selection{
  
  /**
   * @brief  Particle class
   */
  class Particle{

    public : 

      /**
       * @brief  Constructor for MC particles 
       *
       * @param  id of the particle
       * @param  pdg of the particle
       * @param  status code of the particle
       * @param  n_hits number of hits the particle has
       * @param  mass mass of the particle
       * @param  energy total energy of the particle
       * @param  vertex start point of the track
       * @param  end end of the track
       * @param  momentum momentum of the track
       *
       */
      Particle(const int mc_id, 
               const int pdg, 
               const int status, 
               const int n_hits,
               const float mass, 
               const float energy, 
               const TVector3 &vertex, 
               const TVector3 &end, 
               const TVector3 &momentum); 

      /**
       * @brief  Constructor for reconstructed tracks 
       *
       * @param  mc_id_charge MC ID using the charge method
       * @param  mc_id_energy MC ID using the energy method
       * @param  mc_id_hits MC ID using the hits method
       * @param  pdg of the particle
       * @param  n_hits number of hits the particle has
       * @param  kinetic_energy kinetic/visible energy of the particle
       * @param  missing_energy missing energy of the particle
       * @param  mcs_momentum_muon momentum for an escaping muon
       * @param  range_momentum_muon momentum for a contained muon
       * @param  range_momentum_proton momentum for a contained proton
       * @param  length length of the particle
       * @param  range the range of the particle
       * @param  vertex start point of the track
       * @param  end end of the track
       * @param  chi2p chi2 under the proton hypothesis
       * @param  chi2mu chi2 under the muon hypothesis
       * @param  chi2pi chi2 under the pion hypothesis
       * @param  dedx dEdx of the particle
       * @param  residual_range Residual range of the particle
       * @param  pitch the pitch of the track 
       *
       */
      Particle(const int mc_id_charge, 
               const int mc_id_energy, 
               const int mc_id_hits, 
               const int pdg, 
               const int n_hits, 
               const float kinetic_energy, 
               const float missing_energy, 
               const float mcs_momentum_muon, 
               const float range_momentum_muon, 
               const float range_momentum_proton, 
               const float length, 
               const float range, 
               const TVector3 &vertex, 
               const TVector3 &end, 
               const float &chi2p, 
               const float &chi2mu, 
               const float &chi2pi, 
               const std::vector<float> &dedx, 
               const std::vector<float> &residual_range, 
               const std::vector<float> &pitch);

      /**
       * @brief  Constructor for reconstructed showers 
       *
       * @param  pdg of the particle
       * @param  n_hits number of hits the particle has
       * @param  vertex start point of the shower
       * @param  end end of the shower
       * @param  energy energy of the shower
       *
       */
      Particle(const int pdg, 
               const int n_hits, 
               const TVector3 &vertex, 
               const TVector3 &end, 
               const float &energy);

      /**
       * @brief  Get the mass from the pdg code
       *
       * @param  pdg pdg of the particle
       *
       */
      float GetMassFromPdg(const int pdg) const;

      /**
       * @brief  Get the pdg code
       */
      int GetPdgCode() const;

      /**
       * @brief Status code of the particle \n
       *
       * <tt>
       * -1  : Undefined \n
       *  0  : Initial state \n
       *  1  : Stable final state \n
       *  2  : Intermediate state \n
       *  3  : Decayed state \n
       * 10 : Correlated nucleon \n
       * 11 : Nucleon target \n
       * 12 : DIS pre-fragmented hadronic final state \n
       * 13 : Pre-decay resonant state \n
       * 14 : Hadron in nucleus \n
       * 15 : Final state nuclear remnant \n
       * 16 : Nucleon cluster target \n
       * </tt> 
       */
      int GetStatusCode() const;

      /**
       * @brief  Get the number of hits
       */
      int GetNumberOfHits() const;
      /**
       * @brief  Get the mass
       */
      float GetMass() const;

      /**
       * @brief  Get the energy
       */
      float GetEnergy() const;

      /**
       * @brief  Get the chi2 proton
       */
      float GetChi2P() const;

      /**
       * @brief  Get the chi2 muon
       */
      float GetChi2Mu() const;

      /**
       * @brief  Get the chi2 pion
       */
      float GetChi2Pi() const;

      /**
       * @brief  Get the kinetic energy
       */
      float GetKineticEnergy() const;

      /**
       * @brief  Get the length
       */
      float GetLength() const;

      /**
       * @brief  Get the vertex
       */
      TVector3 GetVertex() const;

      /**
       * @brief  Get the end
       */
      TVector3 GetEnd() const;

      /**
       * @brief  Flip the particle's vertex and end point if needed 
       */
      void FlipTrack();

      /**
       * @brief  Get the momentum
       */
      TVector3 GetMomentum() const;

      /**
       * @brief  Get the momentum module
       */
      float GetModulusMomentum() const;

      /**
       * @brief get momentum if escaping muon
       */
      TVector3 GetMCSMomentumMuon() const;

      /**
       * @brief get momentum if contained muon
       */
      TVector3 GetRangeMomentumMuon() const;

      /**
       * @brief get momentum if contained proton
       */
      TVector3 GetRangeMomentumProton() const;
      
      /**
       * @brief  Get the MCParticle id
       */
      int GetMCId() const;

      /**
       * @brief  Get the MCParticle id corresponding to a reco track particle using charge
       */
      int GetMCParticleIdCharge() const;

      /**
       * @brief  Get the MCParticle id corresponding to a reco track particle using energy
       */
      int GetMCParticleIdEnergy() const;

      /**
       * @brief  Get the MCParticle id corresponding to a reco track particle using hits
       */
      int GetMCParticleIdHits() const;

      /**
       * @brief  Get whether the particle has calorimetry
       */
      bool GetHasCalorimetry() const;

      /**
       * @brief  Get whether the particle is from a reconstructed track
       */
      bool GetFromRecoTrack() const;

      /**
       * @brief  Get the cos(theta) of the particle regarding the z direction
       */
      float GetCosTheta() const;

      /**
       * @brief Get whether a track is within the SBND fiducial volume
       */
      bool GetTrackContained() const;

      /**
       * @brief Get whether only 1 end of a track escapes the SBND fiducial volume
       */
      bool GetOneEndTrackContained() const;

    private :

      int      m_mc_id_charge;             ///< mc TrackID corresponding to MCParticle using charge
      int      m_mc_id_energy;             ///< mc TrackID corresponding to MCParticle using energy
      int      m_mc_id_hits;               ///< mc TrackID corresponding to MCParticle using hits
      int      m_mc_id;                    ///< mc TrackID 
      int      m_n_hits;                   ///< number of hits 
      int      m_pdg;                      ///< pdg code
      int      m_status;                   ///< status code
      float    m_mass;                     ///< mass of the particle
      float    m_kinetic_energy;           ///< kinetic energy of the particle
      float    m_missing_energy;           ///< missing energy of the particle
      float    m_energy;                   ///< energy of the particle
      float    m_length;                   ///< length of the particle track
      float    m_range;                    ///< range of the particle track
      float    m_costheta;                 ///< cos(theta) of the particle
      float    m_chi2p;                    ///< chi2 under the proton hypothesis
      float    m_chi2mu;                   ///< chi2 under the muon hypothesis
      float    m_chi2pi;                   ///< chi2 under the pion hypothesis
      float    m_mcs_mom_muon;             ///< multiple coulomb scattering momentum is the particle is an escaping muon
      float    m_range_mom_muon;           ///< range momentum if the particle is a contained muon 
      float    m_range_mom_proton;         ///< range momentum if the particle is a contained proton
      bool     m_has_calorimetry;          ///< whether or not the particle has calorimetry
      bool     m_from_reco_track;          ///< whether the particle is from a reconstructed track
      std::vector<float> m_dedx;           ///< vector of the dedx distribution of the reconstructed track
      std::vector<float> m_residual_range; ///< vector of the residual range distribution of the track
      std::vector<float> m_pitch;          ///< vector of the track pitch distribution of the reconstructed track
      TVector3 m_vertex;                   ///< particle start position
      TVector3 m_end;                      ///< particle end position
      TVector3 m_momentum;                 ///< particle momentum
      float    m_sbnd_border_x;            ///< fiducial border in x for the sbnd detector
      float    m_sbnd_border_y;            ///< fiducial border in y for the sbnd detector
      float    m_sbnd_border_z;            ///< fiducial border in z for the sbnd detector
      float    m_sbnd_offset_x;            ///< offset in x for the sbnd detector
      float    m_sbnd_offset_y;            ///< offset in y for the sbnd detector
      float    m_sbnd_offset_z;            ///< offset in z for the sbnd detector
      float    m_sbnd_length_x;            ///< detector half length in x
      float    m_sbnd_length_y;            ///< detector half length in y
      float    m_sbnd_length_z;            ///< detector half length in z


  }; // Particle
} // Selection
#endif
