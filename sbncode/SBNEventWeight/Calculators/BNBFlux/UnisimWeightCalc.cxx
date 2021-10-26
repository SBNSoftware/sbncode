//Built from FluxUnisimWeightCalc.cxx in ubcode/EventWeight/Calculator/Flux/ directory

#include "FluxCalcPrep.h"


namespace sbn {
  namespace evwgh {

    double FluxWeightCalc::UnisimWeightCalc(double enu, int ptype, int ntype, double randomN, bool noNeg)
    {//same copy from the FluxWeightCalc.cxx in ubsim/EventWeight/
      //Keng Lin June 2021

      double weight = 1;

      int bin = int(enu/0.05); //convert energy (in GeV) into 50 MeV bin


      //  pseudocode:
      //   Scaled Reweighting = ScaleFactor * Reweighting + ( 1 - ScaleFactor) * Central Value
      //
      double scaled_pos = fScalePos*fRWpos[ptype][ntype][bin] + 
        (1-fScalePos)*fCV[ptype][ntype][bin];

      double scaled_neg = fScaleNeg*fRWneg[ptype][ntype][bin] + 
        (1-fScaleNeg)*fCV[ptype][ntype][bin];

      //  pseudocode:
      //   Check value of Random Number for this universe  such that:
      //   If random# > 0
      //       Weight =  1 + ( random# * [ { Scaled Reweighting[pos] / Central Value } - 1 ] )
      //   If random# < 0
      //       Weight =  1 - ( random# * [ { Scaled Reweighting[neg] / Central Value } - 1 ] )
      //
      //   if there is only one systematic histogram offered sub this:
      //   If random# < 0
      //       Weight =  1 - ( random# * [ < 2 - { Scaled Reweighting[pos] / Central Value } > - 1 ] )
      //

      if(randomN > 0){      
        double syst = randomN*((scaled_pos/fCV[ptype][ntype][bin])-1);
        weight = 1 + (syst);

        if(scaled_pos == 0) weight = 1;
      }

      else if(noNeg == true){      
        double syst = randomN*( (2 - (scaled_pos/fCV[ptype][ntype][bin])) - 1);           
        weight = 1 - (syst);      

        if(scaled_pos == 0) weight = 1;

      }
      else{
        double syst = randomN*((scaled_neg/fCV[ptype][ntype][bin])-1);
        weight = 1 - (syst);    

        if(scaled_neg == 0) weight = 1;

      }

      if(fabs(fCV[ptype][ntype][bin]) < 1.e-12) weight = 1;

      if(weight < 0) weight = 1; 
      if(weight > 30) weight = 30; 
      if(!(std::isfinite(weight))){
        std::cout << "UniSim " <<  fParameterSet.fName << " : Failed to get a finite weight" << std::endl;      
        weight = 30;}

      if( (ntype == 0 || ntype == 1) && ptype == 1) weight = 1;//nue/nuebar from pion
      if( (ntype == 1 || ntype == 3) && ptype == 3) weight = 1;//nuebar/numubar from pi or charged kaon

      //      std::cout<<weight<<" ("<<randomN<<") ";
      return weight;

    }

  }  // namespace evwgh
}  // namespace sbn



