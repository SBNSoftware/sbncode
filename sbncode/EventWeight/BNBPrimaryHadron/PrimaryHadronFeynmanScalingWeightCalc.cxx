//Feynmann scaling xsec fit for kaon+
//
//    This code is intended to generate weights for neutrinos originating from K+ decays
//    utilizing Feynman Scaling coming for world data on p+Be 
//    
//    A this code is adopted from the MiniBooNE flux paper and the MiniBooNE reweighting framework 
//    along with code written by Raquel Castillo, Zarko, MiniBooNE (Steve Brice and Mike S), and Athula 
//
//    Current person adding comments and functions is Joseph Zennamo (jaz8600@fnal.gov)
//

#include "../WeightCalcCreator.h"
#include "../WeightCalc.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include <vector>
#include "TH1.h"
#include "TArrayD.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include <TChain.h>

using namespace std;

namespace sbncode {
namespace evwgh {
  class PrimaryHadronFeynmanScalingWeightCalc : public WeightCalc
  {
  public:
    PrimaryHadronFeynmanScalingWeightCalc();
    void Configure(fhicl::ParameterSet const& p);
    std::pair< bool, double > MiniBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand);
    std::pair< bool, double > MicroBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand);
    virtual std::vector<std::vector<double> > GetWeight(art::Event & e);
    

  private:
    CLHEP::RandGaussQ *fGaussRandom;
    std::vector<double> ConvertToVector(TArrayD const* array);
    std::string fGenieModuleLabel;
    std::vector<std::string> fParameter_list;
    float fParameter_sigma;
    int fNmultisims;
    int fprimaryHad;
    std::string fWeightCalc;
    std::string ExternalDataInput;
    double fScaleFactor;
    TFile* file;
    std::vector< std::vector< double > > fWeightArray; 
    std::string fMode;
    std::vector<double> FSKPlusFitVal;
    TMatrixD* FSKPlusFitCov;    


     DECLARE_WEIGHTCALC(PrimaryHadronFeynmanScalingWeightCalc)
  };
  PrimaryHadronFeynmanScalingWeightCalc::PrimaryHadronFeynmanScalingWeightCalc()
  {
  }

  void PrimaryHadronFeynmanScalingWeightCalc::Configure(fhicl::ParameterSet const& p)
  {

    // Here we do all our fhicl file configureation
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");
    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    std::cout << pset.to_string() << std::endl;

    fParameter_list		=   pset.get<std::vector<std::string> >("parameter_list");
    fParameter_sigma		=   pset.get<float>("parameter_sigma");
    fNmultisims			=   pset.get<int>("number_of_multisims");
    fprimaryHad			=   pset.get<int>("PrimaryHadronGeantCode");
    std::string dataInput       =   pset.get< std::string >("ExternalData");
    fWeightCalc                 =   pset.get<std::string>("weight_calculator");
    fMode                       =   pset.get<std::string>("mode");
    fScaleFactor                =   pset.get<double>("scale_factor");
    // Getting External Data:
    //   Now let's get the external data that we will need 
    //   to assess what the central value of the p+Be -> K+
    //   cross section is and then also what are reasonable 
    //   variations around that central value.
    //   First: Pull in the file
    cet::search_path sp("FW_SEARCH_PATH");
    std::string ExternalDataInput = sp.find_file(dataInput);
    file = new TFile(Form("%s",ExternalDataInput.c_str()));
    //  Second: Define what we want to gather from that file which is two fold
    std::vector< std::string > pname; // these are what we will extract from the file
    pname.resize(2); // there are two items

    // We want the parameterization of the Feynman scaling and the covariance matrix
    // which characterizes the uncertainties and how they correlate across the parameterization
    pname[0] = "FS/KPlus/FSKPlusFitVal";
    pname[1] = "FS/KPlus/FSKPlusFitCov";
    TArrayD* FSKPlusFitValArray = (TArrayD*) file->Get(pname[0].c_str());
    //TArrayD is the most annoying format I have ever experienced so let's convert it to a vector
    FSKPlusFitVal = PrimaryHadronFeynmanScalingWeightCalc::ConvertToVector(FSKPlusFitValArray);    
    FSKPlusFitCov = (TMatrixD*) file->Get(pname[1].c_str());
    *(FSKPlusFitCov) *= fScaleFactor*fScaleFactor;
    std::cout << "Scale Factor being applied : " << fScaleFactor << std::endl; 

    // This is done but it is important to note that the parameterization maps such that 
    // FSKPlusFitVal[0] corresponds to the diagonal element FSKPlusFitCov[0][0]
    
    //Prepare random generator
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));
    
    //
    //  This part is very important. You will need more than a single random number
    //  per event. In fact you will need per universe as many random numbers as there
    //  are rows in the covariance matrix. 
    //
    //  To properly perform Cholesky decomposition you'll use a correlated set of random
    //  numbers when creating the smeared cross section methods. Please see documentation
    //  in WeightCalc to help with this.
    //
    fWeightArray.resize(fNmultisims);
        
    for (unsigned int i=0;i<fWeightArray.size();i++) {
      fWeightArray[i].resize(FSKPlusFitCov->GetNcols());      
      if (fMode.find("multisim") != std::string::npos ){
	for(unsigned int j = 0; j < fWeightArray[i].size(); j++){
	  fWeightArray[i][j]=fGaussRandom->shoot(&rng->getEngine(GetName()),0,1.);
	}
      }
      else{
	std::fill(fWeightArray[i].begin(), fWeightArray[i].end(), 1.);
      }
    }
  }
  std::vector<std::vector<double> > PrimaryHadronFeynmanScalingWeightCalc::GetWeight(art::Event & e)
  {

    //Collect the event's Flux information
    //     This specifically deals with the neutrino type and parentage
    art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
    e.getByLabel(fGenieModuleLabel,mcFluxHandle);
    std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;

    
    //Collect event's MC truth information
    //  This specifically deals with the neutrino energy and 
    //  counting how many interactions there are per event 
    //  (neutrino counting is CRITICALLY important for applying the 
    //   correct weights and not ending up with unphysical values)
    art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
    e.getByLabel(fGenieModuleLabel,mctruthHandle);
    std::vector<simb::MCTruth> const& mclist = *mctruthHandle;
    
    //Create a vector of weights for each neutrino 
    std::vector< std::vector<double> > weight;
    weight.resize(mclist.size());
       
    //If there aren't any neutrinos in the event then just return an empty weight vector
    if(mclist.size() == 0 || fluxlist.size() == 0){return weight;}

    // Let's start by iterating through each of the neutrino interactions 
    for(unsigned int inu = 0; inu < fluxlist.size(); inu++){

      // First let's check that the parent of the neutrino we are looking for is 
      //  the particle we intended it to be, if not set all weights to 1
      if (fluxlist[inu].ftptype != fprimaryHad){
	weight[inu].resize(fNmultisims);
	std::fill(weight[inu].begin(), weight[inu].end(), 1);
	continue; //now move on to the next neutrino
      }// Hadronic parent check
          
      //Let's make a weights based on the calculator you have requested 
      if(fMode.find("multisim") != std::string::npos){       
	for (unsigned int i = 0; i < fWeightArray.size() && int(weight[inu].size()) <= fNmultisims; i++) {
	  if(fWeightCalc.find("MicroBooNE") != std::string::npos){
	    
	    //
	    //This way we only have to call the WeightCalc once
	    // 
	    std::pair<bool, double> test_weight =
	      MicroBooNEWeightCalc(fluxlist[inu], fWeightArray[i]);
	    
	    if(test_weight.first){
	      weight[inu].push_back(test_weight.second);
	    }
	  }
	  if(fWeightCalc.find("MiniBooNE") != std::string::npos){

	    //
	    //This way we only have to call the WeightCalc once
	    // 
	    std::pair<bool, double> test_weight =
	      MiniBooNEWeightCalc(fluxlist[inu], fWeightArray[i]);
	    
	    if(test_weight.first){
	      weight[inu].push_back(test_weight.second);
	    }
	  }
	}//Iterate through the number of universes      
      } // make sure we are multisiming

        

    }//Iterating through each neutrino 

    

    return weight;
  }

  //////////////////////////////
  ////////       Auxilary Functions
  //////////////////////////////

  //// 
  //   Use the MiniBooNE Implementation to determine the weight 
  ////
  std::pair< bool, double> PrimaryHadronFeynmanScalingWeightCalc::MiniBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand){
    
  // We now know that we have the hadronic parent that we want so we need to
    // prepare the kinimatic information about the iteraction to go into the Feynman Scaling 
    // calculation
    
    // Computations are based on MiniBooNE code that can be found here:
    // 
    //   http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/~checkout~/...
    //                            miniboone/AnalysisFramework/Utilities/src/PhysicsFunctions.F
    //
    //  Attributed author in MiniBooNE code is Mike Shaevitz 
    //
    //   Useful references:
    //       C.Mariani, G.Cheng, J.Conrad, and M.Shaevitz, PhysRevD.84.114021
    //       MiniBooNE Collaboration, PhysRevD.79.072002
    //
    
    // We need to guard against unphysical parameters 
    bool parameters_pass = false;

    // Get Neutrino Parent Kinimatics       
    double HadronMass;
      
      if(fabs(fprimaryHad) == 321) HadronMass = 0.4937; //Charged Kaon
      else{ 
	throw art::Exception(art::errors::StdException)
	  << "Feynman Scaling is only configured for Charged Kaons ";
      }

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
      

      double E_cm = sqrt(2*ProtonMass + 2*ProtonMass*ProtonVec.P()); // Is this right?
      double gamma = (ProtonE + ProtonMass) / E_cm; 
      double gammaBeta = ProtonVec.P()/E_cm;

      double MaxThreshold;
      if(fabs(fprimaryHad) == 321){
	MaxThreshold = 2.053; // Mike S. Comment: K+: lambda + proton
      }
      else{
	throw art::Exception(art::errors::StdException)
	  << "Feynman Scaling is only configured for Charged Kaons ";	
      }

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
      double c1 = FSKPlusFitVal.at(0);     
      double c2 = FSKPlusFitVal.at(1);     
      double c3 = FSKPlusFitVal.at(2);     
      double c4 = FSKPlusFitVal.at(3);     
      double c5 = FSKPlusFitVal.at(4);     
      double c6 = FSKPlusFitVal.at(5);     
      double c7 = FSKPlusFitVal.at(6);     
      
      double CV = c1*(HadronVec.P()*HadronVec.P()/HadronE)*exp(-1.*c3*pow(fabs(xF),c4) 
							       - c7*pow(fabs(HadronPT*xF),c6)
							       - c2*HadronPT
							       - c5*HadronPT*HadronPT);

      if(CV < 0) CV = 0;

      // Define the variations around that cross section
      //    To do this we will call the a function from WeightCalc
      //    that will generate a set of smeared parameters based 
      //    on the covariance matrix that we imported
      std::vector< double > FSKPlusFitSmeared = WeightCalc::MultiGaussianSmearing(FSKPlusFitVal, FSKPlusFitCov, rand); 
      
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

      //      if(RW < 0) RW = 0;
      
      double weight = 1; 

      if(RW == 0 || CV == 0){
	weight = 1;
      }
      else if(CV < 1.e-12){
	weight = 1;
      }
      else if(fabs(xF) > 1){
	weight = 1;
      }
      else{
	weight *= RW/CV;
      }

      if(weight < 0) weight = 1;
      if(weight >30) weight = 30;
      if(weight != weight) weight = 30; //From MiniBooNE; In Fortran Nan > 30 

      std::pair<bool, double> output(parameters_pass, weight);

      return output; 

  }// Done with the MiniBooNE function


  //// 
  //   Use the MicroBooNE Implementation to determine the weight 
  ////
  std::pair<bool, double> PrimaryHadronFeynmanScalingWeightCalc::MicroBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand){
    
    double weight = 1;
    bool parameters_pass = true;

    std::pair<bool, double> output(parameters_pass, weight);
    
    return output; 

  }// Done with the MicroBooNE function


  //// 
  //  This converts TArrayD to std::vector< double > 
  ///
  std::vector<double> PrimaryHadronFeynmanScalingWeightCalc::ConvertToVector(TArrayD const* array) {
    std::vector<double> v(array->GetSize());
    std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
	      v.begin());
    return v;
  } // ConvertToVector()


  REGISTER_WEIGHTCALC(PrimaryHadronFeynmanScalingWeightCalc)
}
}
