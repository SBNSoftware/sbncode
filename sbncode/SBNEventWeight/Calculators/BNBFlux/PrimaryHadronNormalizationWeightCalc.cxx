//Built from PrimaryHadronNormalizationWeightCalc.cxx in ubcode/EventWeight/Calculator/Flux/BNBPrimaryHadron directory

#include "FluxCalcPrep.h"


namespace sbn {
  namespace evwgh {

    std::pair<bool, double> FluxWeightCalc::PHNWeightCalc(simb::MCFlux flux, float rand){
      // We need to guard against unphysical parameters 
      bool parameters_pass = false;

      double weight = rand+1;

      if(weight > 0) parameters_pass = true;
      else{parameters_pass = false;}

      std::pair<bool, double> output(parameters_pass, weight);

      return output; 

    }

  }  // namespace evwgh
}  // namespace sbn



