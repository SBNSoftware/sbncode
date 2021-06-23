// FluxUnisimWeightCalc.cxx
//
// Handles event weights for GENIE systematics studies
//
// Updated by Marco Del Tutto on Feb 18 2017
//
// Ported from uboonecode to larsim on Feb 14 2017
//   by Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>
//
// Ported to/adapted for SBNCode, Dec 2020, A. Mastbaum
//
// ----
// Ported from ubcode; objects used here are based on
//		SBNEventweight/Calculator/CrossSections/GenieWeightCalc.cxx
//
//  - GetWeight() function is re-introduced but it has bug to be fixed;
//		by Keng Lin Jun. 2021



#include "art/Framework/Principal/Event.h"
#include "nugen/NuReweight/art/NuReweight.h"

#include "nusimdata/SimulationBase/MCFlux.h" //new for flux
#include "nusimdata/SimulationBase/MCTruth.h"
//#include "nusimdata/SimulationBase/GTruth.h"
//#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"

#include "sbncode/SBNEventWeight/Base/WeightCalc.h"
#include "sbncode/SBNEventWeight/Base/WeightCalcCreator.h"


#include "TH1F.h"//need these to read histograms
#include "TFile.h"

namespace sbn {
	namespace evwgh {

		class FluxUnisimWeightCalc : public WeightCalc {
			public:
				FluxUnisimWeightCalc() : WeightCalc() {}

				//Read FHiCL and store the settings for the reweighting environment and the calculator.
				void Configure(fhicl::ParameterSet const& pset,
						CLHEP::HepRandomEngine& engine) override;

				//GetWeight() returns the final weights as a vector
				//	each weight calculation is carried out by WeightCalc() function;
				//inu indicates the ith parameter.
				std::vector<float> GetWeight(art::Event& e, size_t inu) override;

				//Function for evaluating a specific weight
				//enu - neutrino energy from simb::MCTruth; 
				//ptype-  parent particles label: pi, k, k0, mu 
				//ntype - neutrino flavor label: numu, numubar, nue, nuebar 
				//uni - the nth universe
				//noNeg - formulas for weights depending on input histograms.
				double MicroBooNEWeightCalc(double enu, int ptype, int ntype, int uni, bool noNeg);

			private:
				//std::vector<rwgt::NuReweight> rwVector;//reweighter? Delete this, if this is for Genie weights

				std::string fGenieModuleLabel;

				std::string fMode;//should be under fParameterSet.fRWType, tb fixed CHECK
				double fScalePos{}; 
				double fScaleNeg{}; 

				std::vector<double> fWeightArray{};//a vector of random numbers
				//replaced by std::map<EventWeightParameter, std::vector<float> > fParameterSet.fParameterMap
				//CHECK, cannot handld the cituation that there are two sets of fWerightArray; i.e. from AddParameter()

				//contents of the input histograms;
				//[mu/pi/k-/k] x [nue,anue,numu,anumu] x [bin#]
				double fCV[4][4][200];
				double fRWpos[4][4][200];
				double fRWneg[4][4][200];
				bool PosOnly{false};


				DECLARE_WEIGHTCALC(FluxUnisimWeightCalc)
		};



		// Implement the above functions.

		//Configure() need to:
		//- Configure reweighting environment (should be done somewhere else):
		//	- Set fMode;
		//- Configure calculator setting: 
		//	- Load histograms and save them;
		//
		//To do list
		//[ ] locate root files (search_path does not work).

		// FHiCL configuration reads as
		//module_type: "SBNEventWeight"//taken care in WeightManager.h
		//genie_module_label: generator//taken care in WeightManager.h
		//horncurrent: {
		//	<<  Reweighting  Environment >>
		//		  type: FluxUnisim //taken care in EventWeightParameterSet.cxx
		//		  random_seed: 7   //taken care in WeightManager.h
		//		  parameter_list: ["horncurrent"] //processed in EventWeightParameterSet.cxx
		//		  mode: multisim                  //processed in EventWeightParameterSet.cxx
		//		  number_of_multisims: 1000       //processed in EventWeightParameterSet.cxx
		//
		//	<<  Calculator Settings >>
		//		  CentralValue_hist_file: "*.root"
		//		  PositiveSystematicVariation_hist_file: "*.root"
		//		  NegativeSystematicVariation_hist_file: ".root"
		//		  scale_factor_pos: 1
		//		  scale_factor_neg: 1

		//	 << maybe outdated >>
		//		  weight_calculator: "MicroBooNE" //All MicroBooNE should
		//
		//	<< Outdated >>
		//		  use_MiniBooNE_random_numbers: false
		//			 }
		//
		//Notation: taken care - no need to set up here; can be used directly

		void FluxUnisimWeightCalc::Configure(fhicl::ParameterSet const& p,
				CLHEP::HepRandomEngine& engine) {
			//Collect FHiCL parameters:
			//0. Global config
			fGenieModuleLabel = p.get<std::string>("genie_module_label");//use this label to get MC*Handle
			const fhicl::ParameterSet& pset = p.get<fhicl::ParameterSet>(GetName());

			//1. << Reweighting Environment >>
			auto const& pars = pset.get<std::vector<std::string> >("parameter_list");
			for (size_t i=0; i<pars.size(); i++) {
				//Check, no sigma is used in this script.
				fParameterSet.AddParameter(pars[i], 1);//parsigmas[i]);
			}

			fMode = pset.get<std::string>("mode");//3 types: multisim/pmNsigma/fixed

			int number_of_multisims = pset.get<int>("number_of_multisims", 1);

			//CHECK mode is not recognized.
			std::cout<<"\nfPSet configure "<<GetFullName()<<","<<fMode<<","<<number_of_multisims<<std::endl;
			fParameterSet.Configure(GetFullName(), fMode, number_of_multisims);//set fParameterSet members
			//second input is ReweightType rwtype; how does ReweightType work?
			std::cout<<"Outcome: "<<fParameterSet.fName<<","<<fParameterSet.fRWType<<","<<fParameterSet.fNuniverses<<std::endl;
			fParameterSet.Sample(engine);//random_seed is loaded

			//--New version--
			//fWeightArray = random#*1000
			//random# taken care in  sbnobj/sbnobj/Common/SBNEventWeight/EventWeightParameterSet.h 
			
			for (auto const& it : fParameterSet.fParameterMap) {

				std::cout << GetFullName() << ": Setting up " << it.first.fName << std::endl;
				fWeightArray.resize(fParameterSet.fNuniverses);

				for (size_t i=0; i<fParameterSet.fNuniverses; i++) {
					fWeightArray[i] = it.second[i];
				}
			}
			//---
			//generate random numbers;
//			fWeightArray.resize(number_of_multisims);
//			if (fMode == "multisim" )
//				for (double& weight : fWeightArray) weight = CLHEP::RandGaussQ::shoot(&engine, 0, 1.);
//			else
//				for (double& weight : fWeightArray) weight=1.;



			//2. << Calculator Settings >>
			cet::search_path sp("FW_SEARCH_PATH");

			/// Grab the histogram related to the CV
			std::string dataInput1	= pset.get< std::string >("CentralValue_hist_file");
			std::string cvfile		= sp.find_file(dataInput1);
			TFile fcv(Form("%s",cvfile.c_str()));

			/// Grab the histogram related to the variation 
			std::string dataInput2pos	= pset.get< std::string >("PositiveSystematicVariation_hist_file");
			std::string rwfilepos		= sp.find_file(dataInput2pos);
			TFile frwpos(Form("%s", rwfilepos.c_str()));

			std::string dataInput2neg	= pset.get< std::string >("NegativeSystematicVariation_hist_file");
			std::string rwfileneg		= sp.find_file(dataInput2neg);
			TFile frwneg(Form("%s", rwfileneg.c_str()));

			if(dataInput2pos == dataInput2neg) PosOnly = true;//true - for skin depth, use only one variations
    

			int ptype[4] = {1,2,3,4}; //mu, pi, k0, k
			int ntype[4] = {1,2,3,4}; //nue, anue, numu, anumu

			for (int iptyp=0;iptyp<4;iptyp++) {
				for (int intyp=0;intyp<4;intyp++) {
					for (int ibin=0;ibin<200;ibin++) { //Grab events from ibin+1 
						fCV[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (fcv.Get(Form("h5%d%d",ptype[iptyp],ntype[intyp]))))->GetBinContent(ibin+1);
						fRWpos[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (frwpos.Get(Form("h5%d%d",ptype[iptyp],ntype[intyp]))))->GetBinContent(ibin+1);
						fRWneg[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (frwneg.Get(Form("h5%d%d",ptype[iptyp],ntype[intyp]))))->GetBinContent(ibin+1);
					}// energy bin
				}//   type of neutrinos
			}//     type of hadron parent 


			fcv.Close();
			frwpos.Close();
			frwneg.Close(); 

			fScalePos 	= pset.get<double>("scale_factor_pos");
			fScaleNeg 	= pset.get<double>("scale_factor_neg");

		}

		//K: New version requires:
		// GetWeight() needs:
		// - art::Handle< std::vector<simb::MCFlux> > 
		// - art::Handle< std::vector<simb::MCTruth> >
		// - boolean (true - identical in histograms?)

		std::vector<float> FluxUnisimWeightCalc::GetWeight(art::Event& e, size_t inu) {


			//--- Copy over from ubcode
			//
			//Collect the event's Flux information
			//     This specifically deals with the neutrino type and parentage
			//
			// It is better to have MCFlux only;
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
			//	std::vector< std::vector<double> > weight;
			//	weight.resize(mclist.size());

			// No neutrinos in this event
			std::cout<<"CHECK mclist size is "<<mclist.size()<<std::endl;
			std::vector<float> weights(fWeightArray.size(), 0);
			if(mclist.size() == 0) return weights;

			//Iterate through each neutrino in the event

			//Resize vector to the number of universes you want to generate
			//				weight[inu].resize(fNuni);

			//containers for the parent and neutrino type information
			int ptype = std::numeric_limits<int>::max(); 
			int ntype = std::numeric_limits<int>::max();

			// Discover the neutrino parent type
			//     This contains the neutrino's parentage information
			if (      fluxlist[inu].fptype==13  || fluxlist[inu].fptype==-13  ) ptype = 0;//mu
			else if ( fluxlist[inu].fptype==211 || fluxlist[inu].fptype==-211 ) ptype = 1;//pi
			else if ( fluxlist[inu].fptype==130                               ) ptype = 2;//K0
			else if ( fluxlist[inu].fptype==321 || fluxlist[inu].fptype==-321 ) ptype = 3;//K
			else {
				throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<fluxlist[0].fptype<< std::endl;
			}

			// Discover the neutrino type
			//     This contains the neutrino's flavor information
			std::cout<<"CHECK current fluxlist gives fntype = "<<fluxlist[inu].fntype<<std::endl;
			if (      fluxlist[inu].fntype== 12  ) ntype=0;//nue
			else if ( fluxlist[inu].fntype==-12  ) ntype=1;//nuebar
			else if ( fluxlist[inu].fntype== 14  ) ntype=2;//numu
			else if ( fluxlist[inu].fntype==-14  ) ntype=3;//numubar
			else {
				throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ntype "<<fluxlist[inu].fntype<< std::endl;
			}

			// Collect neutrino energy
			//CHECK cant get mclist info; same error as fluxlist
			double enu= mclist[inu].GetNeutrino().Nu().E();
//			std::cout<<"CHECK we got enu = "<<enu<<std::endl;
//			std::cout<<"CHECK we got ntype = "<<ntype<<std::endl;
//			std::cout<<"CHECK we got ptype = "<<ptype<<std::endl;

			//Let's make a weights based on the calculator you have requested 
			if(fMode=="multisim"){//continue calculation with multisim
				std::cout<<"fMode in"<<std::endl;
				for (size_t i=0;i<weights.size();i++) {
					weights[i]=MicroBooNEWeightCalc(enu, ptype, ntype, i, PosOnly);
				}//Iterate through the number of universes      
			}

			std::cout<<"Working on inu:"<<inu<<" with weight "<<weights[inu]<<" size "<<weights.size()<<std::endl;
			std::cout<<"\n The DUMMY VERSION WORKS!********\n\n\n"<<std::endl;

			//--- Copy over from ubcode

			return weights;
		}

		double FluxUnisimWeightCalc::MicroBooNEWeightCalc(double enu, int ptype, int ntype, int uni, bool noNeg)
		{//same copy from the FluxUnisimWeightCalc.cxx in ubsim/EventWeight/
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
			//   Check value of Random Number array [RAND] for this universe [uni] such that:
			//   If RAND[uni] > 0
			//       Weight =  1 + ( RAND[uni] * [ { Scaled Reweighting[pos] / Central Value } - 1 ] )
			//   If RAND[uni] < 0
			//       Weight =  1 - ( RAND[uni] * [ { Scaled Reweighting[neg] / Central Value } - 1 ] )
			//
			//   if there is only one systematic histogram offered sub this:
			//   If RAND[uni] < 0
			//       Weight =  1 - ( RAND[uni] * [ < 2 - { Scaled Reweighting[pos] / Central Value } > - 1 ] )
			//

			if(fWeightArray[uni] > 0){      
				double syst = fWeightArray[uni]*((scaled_pos/fCV[ptype][ntype][bin])-1);
				weight = 1 + (syst);

				if(scaled_pos == 0) weight = 1;
			}

			else if(noNeg == true){      
				double syst = fWeightArray[uni]*( (2 - (scaled_pos/fCV[ptype][ntype][bin])) - 1);           
				weight = 1 - (syst);      

				if(scaled_pos == 0) weight = 1;

			}
			else{
				double syst = fWeightArray[uni]*((scaled_neg/fCV[ptype][ntype][bin])-1);
				weight = 1 - (syst);    

				if(scaled_neg == 0) weight = 1;

			}

			if(fabs(fCV[ptype][ntype][bin]) < 1.e-12) weight = 1;

			if(weight < 0) weight = 1; 
			if(weight > 30) weight = 30; 
			if(!(std::isfinite(weight))){
				std::cout << "UniSim " <<  fParameterSet.fName << " : Failed to get a finite weight" << std::endl;      
				weight = 30;}

			if( (ntype == 0 || ntype == 1) && ptype == 1) weight = 1;//nue from mu or pi
			if( (ntype == 1 || ntype == 3) && ptype == 3) weight = 1;//numu from pi or k

			std::cout<<weight<<" ";
			return weight;

		}

		REGISTER_WEIGHTCALC(FluxUnisimWeightCalc)

	}  // namespace evwgh
}  // namespace sbn



