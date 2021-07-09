#include "FluxCalcPrep.h"


namespace sbn {
	namespace evwgh {

		//Configure everything here! 
		void FluxWeightCalc::Configure(fhicl::ParameterSet const& p,
				CLHEP::HepRandomEngine& engine) {

			fGenieModuleLabel = p.get<std::string>("genie_module_label");//use this label to get MC*Handle
			const fhicl::ParameterSet& pset = p.get<fhicl::ParameterSet>(GetName());

			//0. << Reweighting Environment >>
			
			//skip "random_seed"
			auto const& pars = pset.get<std::vector<std::string> >("parameter_list");
			std::vector< float > parsigmas(pars.size(), 1.0);
//			std::cout<<"Default "<<parsigmas[0]<<std::endl;
//CHECK		std::cout<<"after get_if_present "<<parsigmas[0]<<std::endl;
//CHECK		if(pset.has_key("parameter_sigma") )	std::cout<<" Valid"<<std::endl;
			if (pars.size() != parsigmas.size()) {
				throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
					<< "parameter_list and parameter_sigma length mismatch."
					<< std::endl;
			}

			if(!pset.get_if_present("parameter_sigma", parsigmas)){
				std::cout<<" `parameter_sigma` was not set; now it is set as 1"<<std::endl;
			}

			std::string fMode = pset.get<std::string>("mode");//3 types: multisim/pmNsigma/fixed
			int number_of_multisims = pset.get<int>("number_of_multisims", 1);

			for (size_t i=0; i<pars.size(); i++) {
				//Check, no sigma is used in this script.
				fParameterSet.AddParameter(pars[i], parsigmas[i]);//parsigmas[i]);
			}


			//1. << Calculator Settings >>
			std::string calctype = pset.get< std::string >("calc_type");//Unisim,PrimaryHadronSWCentralSplineVariation,PrimaryHadronFeynmanScaling,PrimaryHadronSanfordWang,PrimaryHadronNormalization
			std::cout<<"Calculator type: "<<calctype<<std::endl;

			fScalePos 	= pset.get<double>("scale_factor_pos");
			if(!pset.get_if_present("scale_factor_neg",fScaleNeg)){
				std::cout<<" `scale_factor_neg` is set to 1."<<std::endl;
				}

			cet::search_path sp("FW_SEARCH_PATH");
			bool validC = true;
			
			//----------------------------
			//-- Non hadrons production --
			//----------------------------
			if( calctype == "Unisim"){//Unisim Calculator
				fParameterSet.Configure(GetFullName(), fMode, number_of_multisims);//set fParameterSet members
				fParameterSet.Sample(engine);//random_seed is loaded at sbncode/Base/WeightManager.h
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
				//Those May07*.root use the following number for histograms
				int cptype[4] = {1,2,3,4}; //mu, pi, k0, k
				int cntype[4] = {1,2,3,4}; //nue, anue, numu, anumu

				for (int iptyp=0;iptyp<4;iptyp++) {
					for (int intyp=0;intyp<4;intyp++) {
						for (int ibin=0;ibin<200;ibin++) { //Grab events from ibin+1 
							fCV[iptyp][intyp][ibin]=(dynamic_cast<TH1F*>	(fcv.Get(Form("h5%d%d",cptype[iptyp],cntype[intyp]))))->GetBinContent(ibin+1);
							fRWpos[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (frwpos.Get(Form("h5%d%d",cptype[iptyp],cntype[intyp]))))->GetBinContent(ibin+1);
							fRWneg[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (frwneg.Get(Form("h5%d%d",cptype[iptyp],cntype[intyp]))))->GetBinContent(ibin+1);
						}// energy bin
					}//   type of neutrinos
				}//type of hadron parent 
				fcv.Close();
				frwpos.Close();
				frwneg.Close(); 

			//-------------
			//-- Hadrons --
			//-------------
			} else if( calctype.compare(0, 13,"PrimaryHadron") == 0){//Hadron Calculators
				std::cout<<"CHECK  Find a Hadron"<<std::endl;
				fParameterSet.Configure(GetFullName(), fMode, 2*number_of_multisims);//double the numbers of universes
				fParameterSet.Sample(engine);

				fprimaryHad	=   pset.get< std::vector<int>>("PrimaryHadronGeantCode");//for Feynman Scaling
				if( calctype == "PrimaryHadronNormalization" ){
					//Nothing to configure

				} else{//Hadron Calculator, slightly complicated

					if( calctype == "PrimaryHadronFeynmanScaling" ){

					}else if( calctype == "PrimaryHadronSanfordWang" ){

					}else if( calctype == "PrimaryHadronSWCentralSplineVariation" ){

					}else	validC = false;
				}

			} else			validC = false;

			if (!validC){
				throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
					<<" calculator "+calctype + "is invalid"
					<<std::endl;
			}


		exit(0);
//
//
//
//			cet::search_path sp("FW_SEARCH_PATH");
//			//need specific files for a calculator.
//			std::string invalidC="";
//			if( calctype == "Unisim"){//Unisim Calculator
//			} else if( calctype.compare(0, 13,"PrimaryHadron") == 0){//Hadron Calculator
//
//
//
//				std::vector< std::string > pname; // these are what we will extract from the file
////				pname.resize(2); // there are two items
//
//				// We want the parameterization of the Feynman scaling and the covariance matrix
//				// which characterizes the uncertainties and how they correlate across the parameterization
//
//				std::string HadronName; 
//				std::string HadronAbriviation;
//				//Note, only look at the first pdg code to determine what histograms to load
//				switch( fprimaryHad[0]){
//					case -321://Kaon-
//						break;
//					case 321://Kaon+
//						pname.push_back("FS/KPlus/FSKPlusFitVal");
//						pname.push_back("FS/KPlus/FSKPlusFitCov");
//						break;
//					case 130://K0L	
//					case 310://K0S
//					case 311://K0
//						pname.push_back("SW/K0s/SWK0sFitVal");
//						pname.push_back("SW/K0s/SWK0sFitCov");
//						break;
//					case -211://pi-
//						HadronName = "PiMinus";
//						HadronAbriviation = "PM";
//						break;
//					case 211://pi+
//						HadronName = "PiPlus";
//						HadronAbriviation = "PP";
//						break;
//					default:
//						throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
//							<<"PrimaryHadronGeantCode "<<fprimaryHad[0]<<" has no assigned calculator"
//							<<std::endl;
//				}
//				if(fabs(fprimaryHad[0]) == 211){//fitInput also configured here
//					pname.push_back( Form("HARPData/%s/%sCrossSection",HadronName.c_str(),HadronAbriviation.c_str())); // Cross Section
//					pname.push_back( Form("HARPData/%s/%scovarianceMatrix",HadronName.c_str(),HadronAbriviation.c_str())); // Covariance Matrix
//					pname.push_back( Form("HARPData/%s/%smomentumBoundsArray",HadronName.c_str(),HadronAbriviation.c_str())); // Momentum Bounds
//					pname.push_back( Form("HARPData/%s/%sthetaBoundsArray",HadronName.c_str(),HadronAbriviation.c_str())); // Theta Bounds
//
//					std::string fitname = Form("SW/%s/SW%sFitVal",HadronName.c_str(),HadronName.c_str()); // Sanford-Wang Fit Parameters
//					TArrayD* SWParamArray = (TArrayD*) Fitfile->Get(fitname.c_str());
//					//CHECK , not finished
//
//				}
//				
//				//  Second: Define what we want to gather from that file which is two fold
//				TArrayD* FitValArray = (TArrayD*) file->Get(pname[0].c_str());
//				//TArrayD is the most annoying format I have ever experienced so let's convert it to a vector
//				FitVal = FluxWeightCalc::ConvertToVector(FitValArray);    
//				FitCov = (TMatrixD*) file->Get(pname[1].c_str());
//				*(FitCov) *= fScalePos*fScalePos;
//				std::cout << "Scale Factor being applied : " << fScalePos << std::endl; 
//
//
//			} 

		}//End of Configure() function


		std::vector<float> FluxWeightCalc::GetWeight(art::Event& e, size_t inu) {

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

			//CHECK, load random# here
			//--New version--
			//fWeightArray has 1000 random number
			//random# taken care in  sbnobj/sbnobj/Common/SBNEventWeight/EventWeightParameterSet.h 

				std::vector<double> fWeightArray{};//a vector of random numbers
			for (auto const& it : fParameterSet.fParameterMap) {

				std::cout << GetFullName() << ": Load random numbers in " << it.first.fName << std::endl;
				fWeightArray.resize(fParameterSet.fNuniverses);

				for (size_t i=0; i<fParameterSet.fNuniverses; i++) {
					fWeightArray[i] = it.second[i];
				}
			}
			// No neutrinos in this event
			std::vector<float> weights(fWeightArray.size(), 0);
			if(mclist.size() == 0) return weights;

			//Iterate through each neutrino in the event

			//Resize vector to the number of universes you want to generate
			//				weight[inu].resize(fNuni);

			//Unisim specific
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
			if (      fluxlist[inu].fntype== 12  ) ntype=0;//nue
			else if ( fluxlist[inu].fntype==-12  ) ntype=1;//nuebar
			else if ( fluxlist[inu].fntype== 14  ) ntype=2;//numu
			else if ( fluxlist[inu].fntype==-14  ) ntype=3;//numubar
			else {
				throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ntype "<<fluxlist[inu].fntype<< std::endl;
			}

			// Collect neutrino energy
			double enu= mclist[inu].GetNeutrino().Nu().E();

			//Let's make a weights based on the calculator you have requested 

			if(fParameterSet.fRWType == EventWeightParameterSet::kMultisim){
//				std::cout<<"<<--- CHECK Filling Weights: "<<std::endl;
				for (size_t i=0;i<weights.size();i++) {
//					if(validate_code > 4) break;
					weights[i]=UnisimWeightCalc(enu, ptype, ntype, fWeightArray[i], PosOnly);
//					if(weights[i]>1) validate_code++;
				}//Iterate through the number of universes      
			}

			std::cout<<"Working on inu:"<<inu<<" with first weight "<<weights[0]<<" size "<<weights.size()<<std::endl;
			std::cout<<"\n The DUMMY VERSION WORKS!********\n\n\n"<<std::endl;

			return weights;
		}

		std::vector<double> FluxWeightCalc::ConvertToVector(TArrayD const* array) {
			std::vector<double> v(array->GetSize());
			std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
					v.begin());
			return v;
		} // ConvertToVector()

		REGISTER_WEIGHTCALC(FluxWeightCalc)
	}  // namespace evwgh
}  // namespace sbn

