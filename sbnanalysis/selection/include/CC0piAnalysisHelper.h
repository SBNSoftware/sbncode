#ifndef CC0PI_ANALYSIS_HELPER_H
#define CC0PI_ANALYSIS_HELPER_H

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
   * @brief  CC0piAnalysisHelper helper class
   */
  class CC0piAnalysisHelper {

    public : 

      typedef std::vector<Particle> ParticleList;
      typedef std::vector<Event>    EventList;
      
      /**
       * @brief  Get the true neutrino energy for CC 0pi interactions
       *
       * @param  track muon track to use in the calculation
       *
       * @return float reconstructed neutrino energy 
       */
      static float GetMCCC0piNeutrinoEnergy(const Event &e);

      /**
       * @brief  Get the reconstructed neutrino energy for CC 0pi interactions
       *
       * @param  track muon track to use in the calculation
       *
       * @return float reconstructed neutrino energy 
       */
      static float GetRecoCC0piNeutrinoEnergy(const Event &e);

      

    private : 
      
      /**
       * @brief  Get the neutrino energy for CC0pi
       */      
      static float GetCC0piNeutrinoEnergy(const ParticleList &particle_list);



  }; // CC0piAnalysisHelper
} // namespace: selection
#endif
