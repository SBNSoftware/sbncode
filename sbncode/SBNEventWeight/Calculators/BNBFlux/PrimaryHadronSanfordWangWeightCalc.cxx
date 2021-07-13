//Built from PrimaryHadronSanfordWang.cxx in ubcode/EventWeight/Calculator/Flux/BNBPrimaryHadron directory

#include "FluxCalcPrep.h"


namespace sbn {
	namespace evwgh {

  std::pair<bool, double> FluxWeightCalc::PHSWWeightCalc(simb::MCFlux flux, std::vector<double> rand){

    // 
    //  Largely built off the MiniBooNE code 
    //  but this is intended to expand beyond it
    //
    //  Edits from MiniBooNE code:
    //
    //  JZ (6/2017) : Remove max weight, set to max_limit of a double
    //  JZ (6/2017) : Changed guards on c9 to be double point precision 
    //   

    // We need to guard against unphysical parameters 
    bool parameters_pass = true;

    // Get Neutrino Parent Kinimatics       
    double HadronMass = 0.4976;

    
//    if(std::find(fprimaryHad.begin(), fprimaryHad.end(), 130) != fprimaryHad.end() ||
//       std::find(fprimaryHad.begin(), fprimaryHad.end(), 310) != fprimaryHad.end() ||
//       std::find(fprimaryHad.begin(), fprimaryHad.end(), 311) != fprimaryHad.end())
//      { 
//	HadronMass = 0.4976; //Neutral Kaon
//      }
//      else{ 
//	throw art::Exception(art::errors::StdException)
//	  << "Sanford-Wang is only configured for Netrual Kaons";
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
      
      // Define the central value cross section
      // Pull out the parameters for the Sanford-Wang Fit
      double c1 = FitVal.at(0);     
      double c2 = FitVal.at(1);     
      double c3 = FitVal.at(2);     
      double c4 = FitVal.at(3);     
      double c5 = FitVal.at(4);     
      double c6 = FitVal.at(5);     
      double c7 = FitVal.at(6);     
      double c8 = FitVal.at(7);     
      double c9 = FitVal.at(8);     

      //Sanford-Wang Parameterization 
      //  Eq 11 from PhysRevD.79.072002
      
      double CV = c1 * pow(HadronVec.P(), c2) * 
	(1. - HadronVec.P()/(ProtonVec.P() - c9)) *
	exp(-1. * c3 * pow(HadronVec.P(), c4) / pow(ProtonVec.P(), c5)) *
	exp(-1. * c6 * HadronVec.Theta() *(HadronVec.P() - c7 * ProtonVec.P() * pow(cos(HadronVec.Theta()), c8)));
      
      // Check taken from MiniBooNE code
      if((HadronVec.P()) > ((ProtonVec.P()) - (c9))){
	CV = 0;
      } 
           
      // Define the variations around that cross section
      //    To do this we will call the a function from WeightCalc
      //    that will generate a set of smeared parameters based 
      //    on the covariance matrix that we imported
      std::vector< double > SWK0FitSmeared = MultiGaussianSmearing(FitVal, FitCov, rand);       
     
      // Pull out the smeared parameters for the Sanford-Wang Fit
      double smeared_c1 = SWK0FitSmeared.at(0);     
      double smeared_c2 = SWK0FitSmeared.at(1);     
      double smeared_c3 = SWK0FitSmeared.at(2);     
      double smeared_c4 = SWK0FitSmeared.at(3);     
      double smeared_c5 = SWK0FitSmeared.at(4);     
      double smeared_c6 = SWK0FitSmeared.at(5);     
      double smeared_c7 = SWK0FitSmeared.at(6);     
      double smeared_c8 = SWK0FitSmeared.at(7);     
      double smeared_c9 = SWK0FitSmeared.at(8);     

      /// Perform the MiniBooNE 
      if(smeared_c1 < 0 || 
	 smeared_c3 < 0 || 
	 smeared_c6 < 0){
	parameters_pass = false;	
      }

      double RW = smeared_c1 * pow(HadronVec.P(), smeared_c2) * 
	(1. - HadronVec.P()/(ProtonVec.P() - smeared_c9)) *
	exp(-1. * smeared_c3 * pow(HadronVec.P(), smeared_c4) / pow(ProtonVec.P(), smeared_c5)) *
	exp(-1. * smeared_c6 * HadronVec.Theta() *(HadronVec.P() - smeared_c7 * ProtonVec.P() * pow(cos(HadronVec.Theta()), smeared_c8)));

      // Check taken from MiniBooNE code
      if((HadronVec.P()) > ((ProtonVec.P()) - (smeared_c9))){
	RW = 0;
      } 
           
      double weight = 1; 

      if(RW < 0 || CV < 0){
	weight = 1;
      }
      else if(fabs(CV) < 1.e-12){
	weight = 1;
      }
      else{
	weight *= RW/CV;
      }

      if(weight < 0) weight = 0; 
      if(weight > 30) weight = 30; 
      if(!(std::isfinite(weight))){
	std::cout << "SW : Failed to get a finite weight" << std::endl;      
	weight = 30;
      }

      std::pair<bool, double> output(parameters_pass, weight);

      return output; 


  }


	}  // namespace evwgh
}  // namespace sbn



