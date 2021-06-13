//Feynmann scaling xsec fit for kaon+
//
//	 This code is ported from ubcode to sbncode at June 2021
//			-- Keng Lin
//
//    This code is intended to generate weights for neutrinos originating from K+ decays
//    utilizing Feynman Scaling coming for world data on p+Be 
//    
//    A this code is adopted from the MiniBooNE flux paper and the MiniBooNE reweighting framework 
//    along with code written by Raquel Castillo, Zarko, MiniBooNE (Steve Brice and Mike S), and Athula 
//
//    Current person adding comments and functions is Joseph Zennamo (jaz8600@fnal.gov)
//

#include "sbncode/SBNEventWeight/Base/WeightCalcCreator.h"
#include "sbncode/SBNEventWeight/Base/WeightCalc.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Persistency/Provenance/ModuleContext.h"
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

namespace sbn{
	namespace evwgh {
		class PrimaryHadronFeynmanScalingWeightCalc : public WeightCalc
		{
			public:
				PrimaryHadronFeynmanScalingWeightCalc() = default;
				void Configure(fhicl::ParameterSet const& p, CLHEP::HepRandomEngine& engine) override;
				//    std::pair< bool, double > MiniBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand);
				std::pair< bool, double > MicroBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand);
//				std::vector<std::vector<double> > GetWeight(art::Event & e) override;
  std::vector<float> GetWeight(art::Event& e, size_t inu) override;

			private:
				std::vector<double> ConvertToVector(TArrayD const* array);

				std::string fGenieModuleLabel{};
				std::vector<std::string> fParameter_list{};
				float fParameter_sigma{};
				int fNmultisims{};
				int fprimaryHad{};
				std::string fWeightCalc{};
				std::string ExternalDataInput{};
				double fScaleFactor{};
				TFile* file{nullptr};
				std::vector< std::vector< double > > fWeightArray{};
				std::string fMode{};
				std::vector<double> FSKPlusFitVal{};
				TMatrixD* FSKPlusFitCov{nullptr};
//				double fSeed{};
//				bool fUseMBRands{false};

				DECLARE_WEIGHTCALC(PrimaryHadronFeynmanScalingWeightCalc)
		};

		void PrimaryHadronFeynmanScalingWeightCalc::Configure(fhicl::ParameterSet const& p,
				CLHEP::HepRandomEngine& engine)
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
					fWeightArray[i].resize(FSKPlusFitCov->GetNcols());      
					if (fMode.find("multisim") != std::string::npos ){
						for(unsigned int j = 0; j < fWeightArray[i].size(); j++){
							fWeightArray[i][j] = CLHEP::RandGaussQ::shoot(&engine, 0, 1.);
						}
					}
					else{
						std::fill(fWeightArray[i].begin(), fWeightArray[i].end(), 1.);
					}
				}

		}

std::vector<float> PrimaryHadronFeynmanScalingWeightCalc::GetWeight(art::Event& e, size_t inu) {

	std::vector<float> weights(fWeightArray.size(), 1);
	return weights;
	}


		std::pair<bool, double> PrimaryHadronFeynmanScalingWeightCalc::MicroBooNEWeightCalc(simb::MCFlux flux, std::vector<double> rand){

			std::pair<bool, double> output(true, 0);

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
