//Sanford-Wang Central Value, Splined HARP data for Pi+
//
//    This code is intended to generate weights for neutrinos originating from Pi+ decays
//    utilizing Sanford-Wang and Spline fits of the HARP data on p+Be->Pi+ 
//    
//    This code is adopted from the MiniBooNE flux paper and the MiniBooNE reweighting framework 
//    along with code written by Raquel Castillo, Zarko, MiniBooNE (Steve Brice and Mike S), and Athula 
//
//    Current person adding comments and functions is Joseph Zennamo (jaz8600@fnal.gov)
//        [this code was designed follow the style of PrimaryHadronFeynmanScalingWeightCalc] 
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
#include "TSpline.h"
#include "TMatrixD.h"
#include "TDecompChol.h"
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
  class PrimaryHadronSWCentralSplineVariationWeightCalc : public WeightCalc
  {
  public:
    PrimaryHadronSWCentralSplineVariationWeightCalc();
    void Configure(fhicl::ParameterSet const& p);
    std::pair< bool, double > MiniBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand);
    std::pair< bool, double > MicroBooNEWeightCalc(simb::MCFlux flux,std::vector<double> rand);
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
    std::string ExternalFitInput;
    double fScaleFactor;
    TFile* file;
    TFile* Fitfile;
    std::vector< std::vector< double > > fWeightArray; 
    std::string fMode;

    TMatrixD* HARPCov;    
    TMatrixD* HARPLowerTriangluarCov;    
    TMatrixD* HARPXSec;    

    std::vector<double> HARPmomentumBounds;
    std::vector<double> HARPthetaBounds;
    std::vector<double> SWParam;

    bool fIsDecomposed;
    


     DECLARE_WEIGHTCALC(PrimaryHadronSWCentralSplineVariationWeightCalc)
  };
  PrimaryHadronSWCentralSplineVariationWeightCalc::PrimaryHadronSWCentralSplineVariationWeightCalc()
  {
  }

  void PrimaryHadronSWCentralSplineVariationWeightCalc::Configure(fhicl::ParameterSet const& p)
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
    std::string fitInput        =   pset.get< std::string >("ExternalFit");
    fWeightCalc                 =   pset.get<std::string>("weight_calculator");
    fMode                       =   pset.get<std::string>("mode");
    fScaleFactor                =   pset.get<double>("scale_factor");
    
    ////////////////////////
    //
    //  We will configure the down stream code to run on the primary 
    //   hadron that we want to be assigning an uncertainty to 
    //
    ///////////////////////
    std::string HadronName; 
    std::string HadronAbriviation;
    if(fprimaryHad == 211){
      HadronName = "PiPlus";
      HadronAbriviation = "PP";
    }
    else if(fprimaryHad == -211){
      HadronName = "PiMinus";
      HadronAbriviation = "PM";
    }
    else{ 
	throw art::Exception(art::errors::StdException)
	  << "sanford-wang is only configured for charged pions ";
      }

    ///////////////////
    //
    // Get HARP Data:
    //   Now let's get the external data that we will need 
    //   to assess what the variations we can make around our
    //   central value cross section 
    //
    //   First: Pull in the file
    ///////////////////
    cet::search_path sp("FW_SEARCH_PATH");
    std::string ExternalDataInput = sp.find_file(dataInput);
    file = new TFile(Form("%s",ExternalDataInput.c_str()));
    //  Second: Define what we want to gather from that file which is two fold
    std::vector< std::string > pname; // these are what we will extract from the file
    pname.resize(4); // there are 4 items

    // Define what we will gather from HARP data
    pname[0] = Form("HARPData/%s/%sCrossSection",HadronName.c_str(),HadronAbriviation.c_str()); // Cross Section
    pname[1] = Form("HARPData/%s/%scovarianceMatrix",HadronName.c_str(),HadronAbriviation.c_str()); // Covariance Matrix
    pname[2] = Form("HARPData/%s/%smomentumBoundsArray",HadronName.c_str(),HadronAbriviation.c_str()); // Momentum Bounds
    pname[3] = Form("HARPData/%s/%sthetaBoundsArray",HadronName.c_str(),HadronAbriviation.c_str()); // Theta Bounds

    HARPXSec = (TMatrixD*) file->Get(pname[0].c_str());
    HARPCov  = (TMatrixD*) file->Get(pname[1].c_str());

    TArrayD* HARPmomentumBoundsArray = (TArrayD*) file->Get(pname[2].c_str());
    HARPmomentumBounds = 
      PrimaryHadronSWCentralSplineVariationWeightCalc::ConvertToVector(HARPmomentumBoundsArray);    

    TArrayD* HARPthetaBoundsArray = (TArrayD*) file->Get(pname[3].c_str());
    HARPthetaBounds = 
      PrimaryHadronSWCentralSplineVariationWeightCalc::ConvertToVector(HARPthetaBoundsArray);    
    
    /////////////////
    //
    //   Extract the Sanford-Wang Fit Parmeters
    //
    ////////////////

    
    std::string ExternalFitInput = sp.find_file(fitInput);
    Fitfile = new TFile(Form("%s",ExternalFitInput.c_str()));

    std::string fitname; // these are what we will extract from the file

    fitname = Form("SW/%s/SW%sFitVal",HadronName.c_str(),HadronName.c_str()); // Sanford-Wang Fit Parameters

    TArrayD* SWParamArray = (TArrayD*) Fitfile->Get(fitname.c_str());
    SWParam = 
      PrimaryHadronSWCentralSplineVariationWeightCalc::ConvertToVector(SWParamArray);    




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
    fWeightArray.resize(2*fNmultisims);
        
    for (unsigned int i=0;i<fWeightArray.size();i++) {
      fWeightArray[i].resize(HARPCov->GetNcols());
      if (fMode.find("multisim") != std::string::npos ){
	for(unsigned int j = 0; j < fWeightArray[i].size(); j++){
	  fWeightArray[i][j]=fGaussRandom->shoot(&rng->getEngine(GetName()),0,1.);
	}//Iterate over the covariance matrix size
      }
      else{
	std::fill(fWeightArray[i].begin(), fWeightArray[i].end(), 1.);
      }
    }//iterate over the number of universes 

    //////
    // Decompose the Covariance Matrix here that way we don't have to do it many times
    /////

    //perform Choleskey Decomposition
    TDecompChol dc = TDecompChol(*(HARPCov));
    if(!dc.Decompose())
      {
    	throw art::Exception(art::errors::StdException)
	  << "Cannot decompose covariance matrix to begin smearing.";
      }

    //Get upper triangular matrix. This maintains the relations in the
    //covariance matrix, but simplifies the structure.
    fIsDecomposed = true;
    HARPLowerTriangluarCov = new TMatrixD(dc.GetU());  
       
  }//End Configure

  std::vector<std::vector<double> > PrimaryHadronSWCentralSplineVariationWeightCalc::GetWeight(art::Event & e)
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
	for (unsigned int i = 0; int(weight[inu].size()) < fNmultisims; i++) {

	  if(fWeightCalc.find("MicroBooNE") != std::string::npos){
	    
	    if(MicroBooNEWeightCalc(fluxlist[inu], fWeightArray[i]).first){
	      weight[inu].push_back(MicroBooNEWeightCalc(fluxlist[inu], fWeightArray[i]).second);
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
  std::pair< bool, double> PrimaryHadronSWCentralSplineVariationWeightCalc::MiniBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand){
    
    bool parameters_pass = true;
    
    double weight = 1;

    ///////////////////////////////
    //
    // We now know that we have the hadronic parent that we want so we need to
    // prepare the kinematic information about the iteraction to go into the Sanford-Wang
    // calculation
    //
    // Computations are based on MiniBooNE code that can be found here:
    // 
    //   http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/~checkout~/...
    //                            miniboone/AnalysisFramework/Utilities/src/PhysicsFunctions.F
    //
    //  >>> Based on MiniBooNE function MultisimMatrix_RawMesonProd
    //
    //
    //
    //   Useful references:
    //       MiniBooNE Collaboration, PhysRevD.79.072002
    //       D.Schmitz Thesis; http://lss.fnal.gov/archive/thesis/2000/fermilab-thesis-2008-26.pdf
    //
    //////////////////////////////

    double c1 = SWParam[0];
    double c2 = SWParam[1];
    double c3 = SWParam[2];
    double c4 = SWParam[3];
    double c5 = SWParam[4];
    double c6 = SWParam[5];
    double c7 = SWParam[6];
    double c8 = SWParam[7];
    double c9 = 1.0; // This isn't in the table but it is described in the text
    
    //  Lay out the event kinimatics 
    double HadronMass;
      
      if(fabs(fprimaryHad) == 211) HadronMass = 0.13957010; //Charged Pion
      else{ 
	throw art::Exception(art::errors::StdException)
	  << "sanford-wang is only configured for charged pions ";
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

      ////////
      //  
      //   Based on MiniBooNE code to evaluate the theta value within below the maximum 
      //    HARP theta coverage. This helps keep the splines well formed and constrains  
      //    the uncertainties at very low neutrino energy
      //
      ////////

      double ThetaOfInterest;
      if(HadronVec.Theta() > 0.195){
	ThetaOfInterest = 0.195;
      }
      else{
	ThetaOfInterest = HadronVec.Theta();
      }


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
      
      //Sanford-Wang Parameterization 
      //  Eq 11 from PhysRevD.79.072002

      double CV = c1 * pow(HadronVec.P(), c2) * 
	(1. - HadronVec.P()/(ProtonVec.P() - c9)) *
	exp(-1. * c3 * pow(HadronVec.P(), c4) / pow(ProtonVec.P(), c5)) *
	exp(-1. * c6 * ThetaOfInterest *(HadronVec.P() - c7 * ProtonVec.P() * pow(cos(ThetaOfInterest), c8)));

      // Check taken from MiniBooNE code
      if(HadronVec.P() > (ProtonVec.P() - c9)){CV = 0;} 

      //////
      //
      // Now that we have our central value based on the Sanford-Wang parameterization we 
      // need to create variations around this value. MiniBooNE did this by performing 
      // spline fits to the HARP data. This is a reproduction of that code:
      // 
      ////

      double* HARPmomentumBins = HARPmomentumBounds.data(); // Convert std::vector to array
      double* HARPthetaBins = HARPthetaBounds.data();       // Convert std::vector to array

      int Ntbins = int(HARPthetaBounds.size()) - 1;
      int Npbins = int(HARPmomentumBounds.size()) - 1;
      
      //
      //  Using the HARP cross section and covariance matrices 
      //  we will now create variations around the measured HARP
      //  meson production cross sections. 
      //
      
      // Important notes of the HARP data: 
      //  The cross section matrix has 13 momentum bins and 6 theta bins 
      //  The covariance matrix contains the correlated uncertainties across
      //  all 78 cross section measurements. 
      //
      //  The covariance matrix encodes the uncertainty on a given HARP 
      //  ANALYSIS BIN instead of being a 3D matrix. 
      //
      //  A HARP analysis bin is defined as 
      //  bin = momentum[theta[]] meaning that:
      //      analysis bin 0 is the zeroth theta and zeroth momentum bin
      //  but analysis bin 27 is the 3rd theta and 5th momentum bin
      // 
      // The first thing to do is convert our cross section matrix into 
      // an std::vector for the analysis bins
      //
      std::vector< double > HARPCrossSectionAnalysisBins;
      HARPCrossSectionAnalysisBins.resize(int(Ntbins*Npbins));

      int anaBin = 0;
      for(int pbin = 0; pbin < Npbins; pbin++){
	for(int tbin = 0; tbin < Ntbins; tbin++){	
	  HARPCrossSectionAnalysisBins[anaBin] = HARPXSec[0][pbin][tbin];
	  anaBin++;
	}
      }
      
      //
      // Now using this we can vary the cross section based on the HARP covariance matrix 
      // this will allow us to reweigh each cross section measurement based on
      // the multigaussian smearing of this matrix. 
      //

      std::vector< double > smearedHARPCrossSectionAnalysisBins = 
	WeightCalc::MultiGaussianSmearing(HARPCrossSectionAnalysisBins, HARPLowerTriangluarCov, fIsDecomposed,rand); 
      HARPCrossSectionAnalysisBins.clear();

      //
      //   Check all the smeared cross sections, if any come out to be negative then 
      //   we will not pass this given parameter set.
      //

      for(int check = 0; check < int(smearedHARPCrossSectionAnalysisBins.size()); check++){
	if(smearedHARPCrossSectionAnalysisBins[check] < 0){ parameters_pass = false;}
      }

      //  
      // With a checked set of smeared cross sections we can convert our analysis bin
      // std::vector back into a TMatrixD, though this is a bit annoying since TMatrices 
      // are an annoying format. 
      //
      TMatrixD*  smearedHARPXSec = new TMatrixD(Npbins,Ntbins);
      
      anaBin = 0;
      for(int pbin = 0; pbin < Npbins; pbin++){
	  for(int tbin = 0; tbin < Ntbins; tbin++){
	  smearedHARPXSec[0][pbin][tbin] = smearedHARPCrossSectionAnalysisBins[anaBin];	  
	  anaBin++;
	}
      }

      //
      // We want to now make a vector of histograms that are 
      //  projections across the cross section matrix, 
      //  this means we want to have vectors where one histogram 
      //  is made for each bin on variable X where the histogram
      //  then contains bins for variable Y. 
      //
      //  This will allow us to do a 2D interpolation across the
      //  full parameter space using cubic spline fits 
      //

      std::vector< TH1F > MomentumBins;
      MomentumBins.resize(Ntbins);
      for(int bin = 0; bin < Ntbins; bin++){
	MomentumBins[bin] = TH1F("HARPp",";;;", Npbins, HARPmomentumBins);
	MomentumBins[bin].SetBit(kCanDelete);
      } 
      
      //Setup vectors of splines for fitting
      std::vector< TSpline3 > SplinesVsMomentum; 
      SplinesVsMomentum.resize(Ntbins);

      for(int tbin = 0; tbin < Ntbins; tbin++){
	  for(int pbin = 0; pbin < Npbins; pbin++){
	
	    //
	    // We want to create spline at fixed theta values, that way we can 
	    // study the cross section at the given meson momentum
	    MomentumBins[tbin].SetBinContent(pbin+1, smearedHARPXSec[0][pbin][tbin]);		    
	    // This will give us the cross section in discreet slices of theta (tbin)
	    // but each histogram in the vector will be the cross section in bins of
	    // momentum that can then be cublic splined  
	    //
	  }
	
	  //
	  // in a given theta slice we can now spline over the cross section entries
	  // in bins of momentum
	  //
	  // Important note about Splines, MiniBooNE used constraints on the
	  // second derivative of the first and final knot point and required 
	  // that they both equal to zero, this is not naturally the case in 
	  // in TSpline3 but it is default in DCSPLC (the Fortran CERN library spline function)
	  // this helps to control the smoothness of the spline and minimizes the variation bin to bin
	  // For TSpine3 this is controlled by:
	  //
	  //  b1 = constrain first knot 1st derivative 
	  //  b2 = constrain first knot 2nd derivative 
	  //  e1 = constrain final knot 1st derivative 
	  //  e2 = constrain final knot 2nd derivative 
	  //
	  //  the numbers that follow are the values you constrain those conditions to
	  //
	  SplinesVsMomentum[tbin] = TSpline3(&(MomentumBins[tbin]),"b2e2",0,0); 

      }
      MomentumBins.clear();
      delete smearedHARPXSec;
      //
      // Now that we have the 1D Splines over momentum we want to know
      // what the cross section is for the meson's momentum in bins of
      // theta that way we can extract the cross section at the meson's
      // theta. 
      //
      TH1F* ThetaBins = new TH1F("HARPt",";;;", Ntbins, HARPthetaBins);
      ThetaBins->SetBit(kCanDelete);
      for(int tbin = 0; tbin < Ntbins; tbin++){
	ThetaBins->SetBinContent(tbin+1, SplinesVsMomentum[tbin].Eval(HadronVec.P()));
      } 

      ////      
      // Now we can spline across the theta bins and we can extract the exact value of the 
      // cross section for the meson's theta value
      //
      //  Options described above.
      /////

      TSpline3* FinalSpline = new TSpline3(ThetaBins,"b2e2",0,0);

      double RW = FinalSpline->Eval(ThetaOfInterest);

      SplinesVsMomentum.clear();
      delete ThetaBins;
      delete FinalSpline;

      weight = RW/CV;

      //
      // These guards are inherited from MiniBooNE code
      //
      // These were defined here: 
      // 
      //   cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/~checkout~/ ... 
      //              miniboone/AnalysisFramework/MultisimMatrix/src/MultisimMatrix.inc    
      //
      //  This forces any negative spline fit to be 1     

      ////////
      //  Possible Bug.
      ////////
      if(weight < 0) weight = 1;  // This seems to be a feature in the MiniBooNE code
                                  //  It looks like the intension is to set this to zero
                                  //  but it is set to 1 before it is set to zero 
      //  This locks in a max weight for all universes 
      if(weight > 30) weight = 30;
      if(weight != weight) weight = 30; //From MiniBooNE; in Fortran Nan > 30
      
     
      std::pair<bool, double> output(parameters_pass, weight);

      return output; 

  }// Done with the MiniBooNE function


  //// 
  //   Use the MicroBooNE Implementation to determine the weight 
  ////
  std::pair<bool, double> PrimaryHadronSWCentralSplineVariationWeightCalc::MicroBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand){
    
    double weight = 1;
    bool parameters_pass = true;

    std::pair<bool, double> output(parameters_pass, weight);
    
    return output; 

  }// Done with the MicroBooNE function


  //// 
  //  This converts TArrayD to std::vector< double > 
  ///
  std::vector<double> PrimaryHadronSWCentralSplineVariationWeightCalc::ConvertToVector(TArrayD const* array) {
    std::vector<double> v(array->GetSize());
    std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
	      v.begin());
    return v;
  } // ConvertToVector()


  REGISTER_WEIGHTCALC(PrimaryHadronSWCentralSplineVariationWeightCalc)
}
}

