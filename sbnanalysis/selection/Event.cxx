#include "Event.hh"
#include "LoadEvents.hh"
#include "EventSelectionHelper.hh"
namespace selection{
  
  Event::Event() : 
    m_mc_particles(m_unfilled),
    m_reco_particles(m_unfilled),
    m_interaction(m_unfilled),
    m_scatter(m_unfilled),
    m_nu_pdg(m_unfilled),
    m_init_pdg(m_unfilled),
    m_is_cc(m_unfilled),
    m_mc_vertex(m_unfilled),
    m_reco_vertex(m_unfilled),
    m_neutrino_energy(m_unfilled),
    m_neutrino_qsqr(m_unfilled), 
    m_sbnd_half_length_x(m_unfilled),
    m_sbnd_half_length_y(m_unfilled),
    m_sbnd_half_length_z(m_unfilled),
    m_sbnd_offset_x(m_unfilled),
    m_sbnd_offset_y(m_unfilled),
    m_sbnd_offset_z(m_unfilled),
    m_sbnd_border_x(m_unfilled),
    m_sbnd_border_y(m_unfilled),
    m_sbnd_border_z(m_unfilled){}

  //------------------------------------------------------------------------------------------ 

  Event::Event(const ParticleList &mc_particles, 
               const ParticleList &reco_particles, 
               const unsigned int interaction, 
               const unsigned int scatter, 
               const int neutrino_pdg, 
               const int initial_pdg, 
               const bool is_cc, 
               const TVector3 &mc_vertex, 
               const TVector3 &reco_vertex, 
               const float neutrino_energy, 
               const float neutrino_qsqr) :
    m_mc_particles(mc_particles),
    m_reco_particles(reco_particles),
    m_interaction(interaction),
    m_scatter(scatter),
    m_nu_pdg(neutrino_pdg),
    m_init_pdg(initial_pdg),
    m_is_cc(is_cc),
    m_mc_vertex(mc_vertex),
    m_reco_vertex(reco_vertex),
    m_neutrino_energy(neutrino_energy),
    m_neutrino_qsqr(neutrino_qsqr){

      // Co-ordinate offset in cm
      m_sbnd_half_length_x = 400;
      m_sbnd_half_length_y = 400;
      m_sbnd_half_length_z = 500;

      m_sbnd_offset_x = 200;
      m_sbnd_offset_y = 200;
      m_sbnd_offset_z = 0;

      m_sbnd_border_x = 10;
      m_sbnd_border_y = 20;
      m_sbnd_border_z = 10;
    }

  //------------------------------------------------------------------------------------------ 
    
  unsigned int Event::CountMCParticlesWithPdg(const int pdg) const{
 
    return this->CountParticlesWithPdg(pdg, m_mc_particles);

  }

  //------------------------------------------------------------------------------------------ 

  unsigned int Event::CountRecoParticlesWithPdg(const int pdg) const{
    
    return this->CountParticlesWithPdg(pdg, m_reco_particles);

  }

  //------------------------------------------------------------------------------------------ 

  bool Event::CheckMCTopology(const TopologyMap &topology) const{
  
    return this->CheckTopology(topology, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------ 

  bool Event::CheckRecoTopology(const TopologyMap &topology) const{
  
    return this->CheckTopology(topology, m_reco_particles);

  }
  
  //------------------------------------------------------------------------------------------ 

  Particle Event::GetMostEnergeticRecoParticle() const{
 
    return this->GetMostEnergeticParticle(m_reco_particles);

  }
  
  //------------------------------------------------------------------------------------------ 

  Particle Event::GetMostEnergeticTrueParticle() const{
 
    return this->GetMostEnergeticParticle(m_mc_particles);

  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetFileId() const{
    return m_file;
  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetId() const{
    return m_id;
  }

  //------------------------------------------------------------------------------------------ 

  ParticleList Event::GetMCParticleList() const{
  
    return m_mc_particles;

  }

  //------------------------------------------------------------------------------------------ 

  ParticleList Event::GetRecoParticleList() const{
  
    return m_reco_particles;

  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetInteractionType() const{
  
    return m_interaction;
  
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int Event::GetPhysicalProcess() const{
  
    return m_scatter;
  
  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetNeutrinoPdgCode() const{
  
    return m_nu_pdg;

  }

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Event::GetMinimumFiducialDimensions() const{
    return TVector3((-m_sbnd_offset_x + m_sbnd_border_x), 
                    (-m_sbnd_offset_y + m_sbnd_border_y), 
                    (-m_sbnd_offset_z + m_sbnd_border_z));
  }

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Event::GetMaximumFiducialDimensions() const{
    return TVector3((m_sbnd_half_length_x - m_sbnd_offset_x - m_sbnd_border_x), 
                    (m_sbnd_half_length_y - m_sbnd_offset_y - m_sbnd_border_y), 
                    (m_sbnd_half_length_z - m_sbnd_offset_z - m_sbnd_border_z));
  }
  
  //------------------------------------------------------------------------------------------ 
  
  bool Event::IsSBNDTrueFiducial() const{
       
    // Check the neutrino interaction vertex is within the fiducial volume      
     float nu_vertex_x = m_mc_vertex[0];                        
     float nu_vertex_y = m_mc_vertex[1];                        
     float nu_vertex_z = m_mc_vertex[2];                
     float min_fid_x = Event::GetMinimumFiducialDimensions()[0];
     float min_fid_y = Event::GetMinimumFiducialDimensions()[1];
     float min_fid_z = Event::GetMinimumFiducialDimensions()[2];
     float max_fid_x = Event::GetMaximumFiducialDimensions()[0];
     float max_fid_y = Event::GetMaximumFiducialDimensions()[1];
     float max_fid_z = Event::GetMaximumFiducialDimensions()[2];
                                                                                 
     if (    (nu_vertex_x > max_fid_x)  
          || (nu_vertex_x < min_fid_x)
          || (nu_vertex_y > max_fid_y)
          || (nu_vertex_y < min_fid_y)
          || (nu_vertex_z > max_fid_z)
          || (nu_vertex_z < min_fid_z)) return false;

     return true;
  }

  //------------------------------------------------------------------------------------------ 
 
  bool Event::AllRecoContained() const{
    for(const Particle &p : m_reco_particles){
      if(p.GetFromRecoTrack() && !p.GetTrackContained()) return false;
    }
    return true;
  }
  
  //------------------------------------------------------------------------------------------ 

  int Event::NumberOfEscapingRecoParticles() const{
    return this->Event::NumberOfEscapingParticles(this->Event::GetRecoParticleList());
  }
  
  //------------------------------------------------------------------------------------------ 
      
  int Event::NumberOfEscapingMCParticles() const{
    return this->Event::NumberOfEscapingParticles(this->Event::GetMCParticleList());
  }
  
  //------------------------------------------------------------------------------------------ 
      
  int Event::NumberOfEscapingParticles(const ParticleList &particles) const{
    int escaping = 0;
    for(const Particle &p : particles){
      if(!p.GetFromRecoTrack() || !p.GetHasCalorimetry()) continue;
      if(p.GetOneEndTrackContained()) escaping++;
    }
    return escaping;
  }
  
  //------------------------------------------------------------------------------------------ 

  bool Event::GetIsCC() const{
  
    return m_is_cc;
  
  }

  //------------------------------------------------------------------------------------------ 

  TVector3 Event::GetMCNuVertex() const{
  
    return m_mc_vertex;
  
  }

  //------------------------------------------------------------------------------------------ 

  TVector3 Event::GetRecoNuVertex() const{
  
    return m_reco_vertex;
  
  }

  //------------------------------------------------------------------------------------------ 
  
  float Event::GetTrueNuEnergy() const{
  
    return m_neutrino_energy;

  }
  float Event::GetTrueNuQ2() const{
  
    return m_neutrino_qsqr;

  }
  
  //------------------------------------------------------------------------------------------ 

  unsigned int Event::CountParticlesWithPdg(const int pdg, const ParticleList &particle_list) const{
    unsigned int particle_counter = 0;
    for(unsigned int i = 0; i < particle_list.size(); ++i) if(particle_list[i].GetPdgCode() == pdg && particle_list[i].GetNumberOfHits() >= 5) particle_counter++;
    return particle_counter;
  }

  //------------------------------------------------------------------------------------------ 
  
  bool Event::CheckTopology(const TopologyMap &topology, const ParticleList &particle_list) const{
    // Loop over the map
    for( TopologyMap::const_iterator it = topology.begin(); it != topology.end(); ++it ){
      // Define temporary variables for the current map element
      std::vector< int > pdg_codes = it->first; 
      int n_total                  = it->second;
      // Count the number of particles in the current event with the same PDG codes 
      // as given by the chosen topology
      int counter = 0;
      // Loop over particles in current event
      for(unsigned int i = 0; i < pdg_codes.size(); ++i){
        counter += this->CountParticlesWithPdg(pdg_codes[i], particle_list);
      }
      if(counter != n_total) return false;
    }
    return true;
  }
  
  //------------------------------------------------------------------------------------------
 
  Particle Event::GetMostEnergeticParticle(const ParticleList &particle_list) const{
    float highest_energy   = -std::numeric_limits<float>::max();
    unsigned int energy_id = std::numeric_limits<unsigned int >::max();
    for(unsigned int i = 0; i < particle_list.size(); ++i){
      if(!particle_list[i].GetHasCalorimetry()) continue;
      if(particle_list[i].GetEnergy() > highest_energy) energy_id = i;
    }
    return particle_list[energy_id];
  }
  
} // selection
