//Built from PrimaryHadronFeynmanScaling.cxx in ubcode/EventWeight/Calculator/Flux/BNBPrimaryHadron directory

#include "FluxCalcPrep.h"


namespace sbn {
  namespace evwgh {

    std::pair<bool, double> FluxWeightCalc::PHFSWeightCalc(simb::MCFlux flux, std::vector<float> rand){

      // 
      //  Largely built off the MiniBooNE code 
      //  but this is intended to expand beyond it
      //
      //  Edits from MiniBooNE code:
      //
      //  JZ (6/2017) : Remove max weight, set to max_limit of a double
      //  JZ (6/2017) : Fixed typo in the E_cm calculation per Mike S. suggestions 
      //   

      // We need to guard against unphysical parameters 
      bool parameters_pass = false;

      // Get Neutrino Parent Kinimatics       
      double HadronMass = 0.4937;

//      if(fabs(fprimaryHad) == 321) HadronMass = 0.4937; //Charged Kaon
//      else{ 
//        throw art::Exception(art::errors::StdException)
//          << "Feynman Scaling is only configured for Charged Kaons ";
//      }

      TLorentzVector HadronVec; 
      double HadronPx = flux.ftpx;
      double HadronPy = flux.ftpy;
      double HadronPz = flux.ftpz;
      double HadronE  = sqrt(HadronPx*HadronPx + 
          HadronPy*HadronPy + 
          HadronPz*HadronPz + 
          HadronMass*HadronMass);
      HadronVec.SetPxPyPzE(HadronPx,HadronPy,HadronPz,HadronE);

      // Get Initial Proton Kinitmatics 
      //   CURRENTLY GSimple flux files drop information about 
      //   the initial state proton, but that this 
      TLorentzVector ProtonVec;
      double ProtonMass = 0.9382720;
      double ProtonPx = 0;
      double ProtonPy = 0;
      double ProtonPz = 8.89; //GeV
      double ProtonE  = sqrt(ProtonPx*ProtonPx +
          ProtonPy*ProtonPy +
          ProtonPz*ProtonPz +
          ProtonMass*ProtonMass);
      ProtonVec.SetPxPyPzE(ProtonPx,ProtonPy,ProtonPz,ProtonE);

      //
      //  From Mike S. :
      //  The fitting was done without this error so I think that you need to use:
      //           ecm = sqrt(2.*mp**2+2.*mp*eb)
      //
      double E_cm = sqrt(2*ProtonMass*ProtonMass + 2*ProtonMass*ProtonVec.E());

      double gamma = (ProtonE + ProtonMass) / E_cm; 
      double gammaBeta = ProtonVec.P()/E_cm;

      double MaxThreshold = 2.053;
//      if(fabs(fprimaryHad) == 321){
//        MaxThreshold = 2.053; // Mike S. Comment: K+: lambda + proton
//      }
//      else{
//        throw art::Exception(art::errors::StdException)
//          << "Feynman Scaling is only configured for Charged Kaons ";  
//      }

      double E_cm_max = (E_cm*E_cm
          +HadronMass*HadronMass
          -MaxThreshold*MaxThreshold)/(2*E_cm);
      double p_cm_max = sqrt(E_cm_max*E_cm_max-HadronMass*HadronMass);

      double HadronPT = HadronVec.Pt();
      double HadronPparallel = HadronVec.P()*cos(HadronVec.Theta());
      double ParallelCMp = gamma*HadronPparallel - gammaBeta*HadronE;

      double xF;

      xF = ParallelCMp / p_cm_max; 

      // Define the scale central value cross section

      // Pull out the parameters for the Feynman Scaling 
      double c1 = FitVal.at(0);     
      double c2 = FitVal.at(1);     
      double c3 = FitVal.at(2);     
      double c4 = FitVal.at(3);     
      double c5 = FitVal.at(4);     
      double c6 = FitVal.at(5);     
      double c7 = FitVal.at(6);     

      double CV = c1*(HadronVec.P()*HadronVec.P()/HadronE)*exp(-1.*c3*pow(fabs(xF),c4) 
          - c7*pow(fabs(HadronPT*xF),c6)
          - c2*HadronPT
          - c5*HadronPT*HadronPT);

      if(CV < 0) CV = 0;
      if(fabs(xF) > 1) CV = 0;
      // Define the variations around that cross section
      //    To do this we will call the a function from WeightCalc
      //    that will generate a set of smeared parameters based 
      //    on the covariance matrix that we imported
      std::vector< double > FSKPlusFitSmeared = MultiGaussianSmearing(FitVal, FitCov, rand);//Defined in sbncode/SBNEventWeight/Base/SmearingUtils.h

      // Pull out the parameters for the Feynman Scaling 
      double smeared_c1 = FSKPlusFitSmeared.at(0);     
      double smeared_c2 = FSKPlusFitSmeared.at(1);     
      double smeared_c3 = FSKPlusFitSmeared.at(2);     
      double smeared_c4 = FSKPlusFitSmeared.at(3);     
      double smeared_c5 = FSKPlusFitSmeared.at(4);     
      double smeared_c6 = FSKPlusFitSmeared.at(5);     
      double smeared_c7 = FSKPlusFitSmeared.at(6);     

      if(smeared_c1 > 0 && 
          smeared_c2 > 0 && 
          smeared_c3 > 0 && 
          smeared_c5 > 0 && 
          smeared_c7 > 0){
        parameters_pass = true;  
      }

      double RW = smeared_c1*(HadronVec.P()*HadronVec.P()/HadronE)*exp(-1.*smeared_c3*pow(fabs(xF),smeared_c4) 
          - smeared_c7*pow(fabs(HadronPT*xF),smeared_c6)
          - smeared_c2*HadronPT
          - smeared_c5*HadronPT*HadronPT);

      if(RW < 0) RW = 0;
      if(fabs(xF) > 1) RW = 0;

      double weight = 1; 

      if(RW < 0 || CV < 0){
        weight = 1;
      }
      else if(CV < 1.e-12){
        weight = 1;
      }
      else{
        weight *= RW/CV;
      }

      if(weight < 0) weight = 0;
      if(weight > 30) weight = 30;
      if(!(std::isfinite(weight))){
        std::cout << "FS : Failed to get a finite weight" << std::endl; 
        weight = 30;}

      std::pair<bool, double> output(parameters_pass, weight);



      return output; 

    }



  }  // namespace evwgh
}  // namespace sbn



