//Sanford-Wang xsec fit for neutral kaon
//
//    This code is intended to generate weights for neutrinos originating from K0 decays
//    utilizing a Sanford-Wang fit of world data  
//    
//    A this code is adopted from the MiniBooNE flux paper and the MiniBooNE reweighting framework 
//    along with code written by Raquel Castillo, Zarko, MiniBooNE (Steve Brice and Mike S), and Athula 
//
//    Current person adding comments and functions is Joseph Zennamo (jaz8600@fnal.gov)
//

#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"

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
  class PrimaryHadronNormalizationWeightCalc : public WeightCalc
  {
  public:
    PrimaryHadronNormalizationWeightCalc();
    void Configure(fhicl::ParameterSet const& p);
    std::pair< bool, double > MiniBooNEWeightCalc(simb::MCFlux flux, double rand);
    std::pair< bool, double > MicroBooNEWeightCalc(simb::MCFlux flux, double rand);
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
    double fScaleFactor;
    TFile* file;
    std::vector< double > fWeightArray; 
    std::string fMode;

     DECLARE_WEIGHTCALC(PrimaryHadronNormalizationWeightCalc)
  };
  PrimaryHadronNormalizationWeightCalc::PrimaryHadronNormalizationWeightCalc()
  {
  }

  void PrimaryHadronNormalizationWeightCalc::Configure(fhicl::ParameterSet const& p)
  {

    // Here we do all our fhicl file configureation
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");
    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    std::cout << pset.to_string() << std::endl;

    fParameter_list		=   pset.get<std::vector<std::string> >("parameter_list");
    fParameter_sigma		=   pset.get<float>("parameter_sigma");
    fNmultisims			=   pset.get<int>("number_of_multisims");
    fprimaryHad			=   pset.get< int >("PrimaryHadronGeantCode");
    fWeightCalc                 =   pset.get<std::string>("weight_calculator");
    fMode                       =   pset.get<std::string>("mode");
    fScaleFactor                =   pset.get<double>("scale_factor");
    
    //Prepare random generator
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));
    
    //
    // This is the meat of this, select random numbers that will be the 
    //   the reweighting
    //
    fWeightArray.resize(2*fNmultisims);
           
    for (unsigned int i=0;i<fWeightArray.size();i++) {
      if (fMode.find("multisim") != std::string::npos ){
	fWeightArray[i]=fGaussRandom->shoot(&rng->getEngine(GetName()),1,1.);
      }
      
      else{
	fWeightArray[i] = 1.;
      }
    }
    
  }// End Configure

  std::vector<std::vector<double> > PrimaryHadronNormalizationWeightCalc::GetWeight(art::Event & e)
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
      // 
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
  std::pair< bool, double> PrimaryHadronNormalizationWeightCalc::MiniBooNEWeightCalc(simb::MCFlux flux, double rand){
    
    // We need to guard against unphysical parameters 
    bool parameters_pass = false;

    double weight = rand; 

    if(weight > 0) parameters_pass = true;
    else{parameters_pass = false;}
    
    
      std::pair<bool, double> output(parameters_pass, weight);
      
      return output; 

  }// Done with the MiniBooNE function


  //// 
  //   Use the MicroBooNE Implementation to determine the weight 
  ////
  std::pair<bool, double> PrimaryHadronNormalizationWeightCalc::MicroBooNEWeightCalc(simb::MCFlux flux, double rand){
    
    double weight = 1;
    bool parameters_pass = true;

    std::pair<bool, double> output(parameters_pass, weight);
    
    return output; 

  }// Done with the MicroBooNE function


  //// 
  //  This converts TArrayD to std::vector< double > 
  ///
  std::vector<double> PrimaryHadronNormalizationWeightCalc::ConvertToVector(TArrayD const* array) {
    std::vector<double> v(array->GetSize());
    std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
	      v.begin());
    return v;
  } // ConvertToVector()


  REGISTER_WEIGHTCALC(PrimaryHadronNormalizationWeightCalc)
}
}

