#include "../include/CC0piAnalysisHelper.h"

namespace selection{

  //------------------------------------------------------------------------------------------

  float CC0piAnalysisHelper::GetMCCC0piNeutrinoEnergy(const Event &e) {
    return GetCC0piNeutrinoEnergy(e.GetMCParticleList()); 
  }

  //------------------------------------------------------------------------------------------
 
  float CC0piAnalysisHelper::GetRecoCC0piNeutrinoEnergy(const Event &e) {
    return GetCC0piNeutrinoEnergy(e.GetRecoParticleList()); 
  }
  
  //------------------------------------------------------------------------------------------
 
  float CC0piAnalysisHelper::GetCC0piNeutrinoEnergy(const ParticleList &particle_list) {
    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float V     = 0.02950;   // Nucleon removal energy, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float reco=0, e, p, cth;   // track variables
    
    // Vector of z direction
    TVector3 z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if( particle_list[i].GetPdgCode() == 13 ){
        // Get the values needed
        e    = particle_list[i].GetEnergy();
        p    = particle_list[i].GetMomentum().Mag();
        cth  = (1/p) * (particle_list[i].GetMomentum()).Dot(z);
        
        reco = (1/(m_n - V - e + p*cth))*((m_n - V)*e - (m_mu*m_mu*0.5) + m_n*V - (V*V*0.5) + (m_p*m_p - m_n*m_n)*0.5);
        return reco;
      }
    }
    return reco;
  }

  //------------------------------------------------------------------------------------------
  
} // selection
