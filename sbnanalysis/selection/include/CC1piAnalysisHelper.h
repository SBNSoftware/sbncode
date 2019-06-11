#ifndef CC1PI_ANALYSIS_HELPER_H
#define CC1PI_ANALYSIS_HELPER_H

#include <string>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "GeneralAnalysisHelper.h"
#include "EventSelectionTool.h"
#include "Event.h"
#include "Particle.h"

namespace selection{
  
  /**
   * @brief  CC1piAnalysisHelper helper class
   */
  class CC1piAnalysisHelper {

    public : 

      typedef std::vector<Particle> ParticleList;
      typedef std::vector<Event>    EventList;
      
      /**                                                              
       * @brief  Gives the number of times the Muon is the longest track and the pion or proton is the second longest track.                                     
       *
       * @param  e current event
       * @param  signal_map_topology the chosen topology 
       * @param  count_longest count the number of times the muon is the longest
       * @param  count_second_longest count the number of times the proton or pion is the second longest
       *
       * @return Particle matrix with statistics based on length and pdg
       */
      static ParticleMatrix LengthBasedStatistics(const Event &e, const TopologyMap signal_map_topology, ParticleMatrix &count_longest, ParticleMatrix &count_second_longest);
      
      /**                                                              
       * @brief Returns the MC momentum transfer for cc1pi     
       * @param pdg                        
       */
      static float GetMCQ2WithPdg(const Event &e, const int pdg);
      /**                                                              
       * @brief Returns the Reco momentum transfer for cc1pi     
       * @param pdg                        
       */
      static float GetRecoQ2WithPdg(const Event &e, const int pdg);

      /**                                                              
       * @brief Returns MC Energy of the particle with longest track for cc1p                         
       */
      static float GetMCEnergyLongest(const Event &e);
      /**                                                              
       * @brief Returns Reco Energy  of the particle with longest track for cc1p    
       */
      static float GetRecoEnergyLongest(const Event &e);
      /**                                                              
       * @brief Returns Reco Energy  of the particle with the second longest track for cc1p    
       */
      static float GetMCEnergySecondLongest(const Event &e);
      /**                                                              
       * @brief Returns Reco Energy  of the particle with the second longest track for cc1p    
       */
      static float GetRecoEnergySecondLongest(const Event &e);
      /**                                                              
       * @brief Returns MC Energy  of the particle with longest track for cc1p    
       */
      static float GetMCKineticEnergyLongest(const Event &e);
      /**                                                              
       * @brief Returns Reco Energy  of the particle with longest track for cc1p    
       */
      static float GetRecoKineticEnergyLongest(const Event &e);
      /**                                                              
       * @brief Returns MC kinetic Energy  of the particle with the second longest track for cc1p    
       */
      static float GetMCKineticEnergySecondLongest(const Event &e);
      /**                                                              
       * @brief Returns Reco kinetic Energy  of the particle with the second longest track for cc1p    
       */
      static float GetRecoKineticEnergySecondLongest(const Event &e);
      /**                                                              
       * @brief Returns MC momentum module  of the particle with longest track for cc1p    
       */
      static float GetMCModulusMomentumLongest(const Event &e);
      /**                                                              
       * @brief Returns Reco momentum module  of the particle with longest track for cc1p    
       */
      static float GetRecoModulusMomentumLongest(const Event &e);
      /**                                                              
       * @brief Returns momentum module  of the particle with the second longest track for cc1p    
       */
      static float GetMCModulusMomentumSecondLongest(const Event &e);
      /**                                                              
       * @brief Returns momentum module of the particle with the second longest track for cc1p    
       */
      static float GetRecoModulusMomentumSecondLongest(const Event &e);
      /**                                                              
       * @brief Returns the energy of the delta particle produced in a resonance                                                     
       */
      static float GetMCDeltaEnergy(const Event &e);
      
      /**
       * @brief  Get the true neutrino energy for CC 1pi interactions
       *
       * @param  track muon track to use in the calculation
       *
       * @return float reconstructed neutrino energy 
       */
      static float GetMCCC1piNeutrinoEnergy(const Event &e);
      
      /**
       * @brief  Get the reconstructed neutrino energy for CC 1pi interactions
       *
       * @param  event
       *
       * @return float reconstructed neutrino energy 
       */
      static float GetRecoCC1piNeutrinoEnergy(const Event &e);
      
      /**
       * @brief  Get the reconstructed neutrino energy for CC 1pi interactions (METHOD 2)
       *
       * @param  event
       *
       * @return float reconstructed neutrino energy 
       */
      static float GetRecoCC1piNeutrinoEnergyMethod2(const Event &e);

    private : 

      /**                                                              
       * @brief Returns the energy of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      static float GetEnergyLongest(const Event &e, const ParticleList &particle_list);
      /**                                                              
       * @brief Returns the energy of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      static float GetEnergySecondLongest(const Event &e, const ParticleList &particle_list);
      /**                                                              
       * @brief Returns the kinetic energy of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      static float GetKineticEnergyLongest(const Event &e, const ParticleList &particle_list);
      /**                                                              
       * @brief Returns the kinetic energy of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      static float GetKineticEnergySecondLongest(const Event &e, const ParticleList &particle_list);
      /**                                                              
       * @brief Returns the momentum module of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      static float GetModulusMomentumLongest(const Event &e, const ParticleList &particle_list);
      /**                                                              
       * @brief Returns the momentum module of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      static float GetModulusMomentumSecondLongest(const Event &e, const ParticleList &particle_list);
      
      /**                                                              
       * @brief  Get the energy of the delta particle produced in a resonance
       *
       * @param  event
       * @param  particle_list                                     
       *
       * @return delta energy
       */
      static float GetDeltaEnergy(const Event &e, const ParticleList &particle_list);
      
      /**                                                              
       * @brief  Get the energy of the delta particle produced in a resonance with a proton in the final state
       *
       * @param  event
       * @param  particle_list                                     
       *
       * @return delta energy
       */
      static float GetDeltaEnergy_p(const ParticleList &particle_list);
      
      /**
       * @brief  Get the neutrino energy for CC1pi
       *
       * @param  particle list
       *
       * @return neutrino energy
       */

      static float GetCC1piNeutrinoEnergy(const ParticleList &particle_list);
      /**
       * @brief  Get the neutrino energy for CC1pi (METHOD 2)
       *
       * @param  particle list
       *
       * @return neutrino energy
       */
      static float GetCC1piNeutrinoEnergyMethod2(const ParticleList &particle_list);

  }; // CC1piAnalysisHelper
} // namespace: selection
#endif
