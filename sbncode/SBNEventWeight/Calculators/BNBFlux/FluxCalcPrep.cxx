#include "FluxCalcPrep.h"


namespace sbn {
	namespace evwgh {

		//Configure everything here! 
		void FluxWeightCalc::Configure(fhicl::ParameterSet const& p,
				CLHEP::HepRandomEngine& engine) {
			std::cout<<"SBNEventWeight : configure "<< GetName() << std::endl;

			fGenieModuleLabel = p.get<std::string>("genie_module_label");//use this label to get MC*Handle
			const fhicl::ParameterSet& pset = p.get<fhicl::ParameterSet>(GetName());

			//0. << Reweighting Environment >>
			
			//skip "random_seed"
			auto const& pars = pset.get<std::vector<std::string> >("parameter_list");
			std::vector< float > parsigmas(pars.size(), 1.0);
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

			fParameterSet.Configure(GetFullName(), fMode, number_of_multisims);
			for (size_t i=0; i<pars.size(); i++) {
				//Check, no sigma is used in this script.
				fParameterSet.AddParameter(pars[i], parsigmas[i]);//parsigmas[i]);
			}


			//1. << Calculator Settings >>
			CalcType = pset.get< std::string >("calc_type");//Unisim,PrimaryHadronSWCentralSplineVariation,PrimaryHadronFeynmanScaling,PrimaryHadronSanfordWang,PrimaryHadronNormalization
			std::cout<<" Calculator type: "<<CalcType<<std::endl;

			fScalePos 	= pset.get<double>("scale_factor_pos");
			if(!pset.get_if_present("scale_factor_neg",fScaleNeg)){
				std::cout<<"  auto-asignment: scale_factor_neg = 1."<<std::endl;
				}

			cet::search_path sp("FW_SEARCH_PATH");
			bool validC = true;
			
			//----------------------------
			//-- Non hadrons production --
			//----------------------------
			
//	CHECK				std::string dataInput       =   pset.get< std::string >("ExternalData");
//					std::string ExternalDataInput = sp.find_file(dataInput);
//					//std::cout<<__LINE__<<" CHECK find "<<ExternalDataInput<<std::endl;
//					TFile* file = new TFile(Form("%s",ExternalDataInput.c_str()));

			fParameterSet.Sample(engine);//random_seed is loaded at sbncode/Base/WeightManager.h
			if( CalcType == "Unisim"){//Unisim Calculator

				std::string dataInput1	= pset.get< std::string >("CentralValue_hist_file");
//				std::cout<<__LINE__<<" CHECK find "<<dataInput1<<std::endl;
				std::string cvfile_path	= sp.find_file(dataInput1);
//				std::cout<<__LINE__<<" CHECK "<<cvfile_path<<std::endl;
				TFile fcv(Form("%s",cvfile_path.c_str()));

				/// Grab the histogram related to the variation 
				std::string dataInput2pos	= pset.get< std::string >("PositiveSystematicVariation_hist_file");
				std::string rwfilepos		= sp.find_file(dataInput2pos);
				TFile frwpos(Form("%s", rwfilepos.c_str()));

				std::string dataInput2neg	= pset.get< std::string >("NegativeSystematicVariation_hist_file");
				std::string rwfileneg		= sp.find_file(dataInput2neg);
				TFile frwneg(Form("%s", rwfileneg.c_str()));


				if(dataInput2pos == dataInput2neg) PosOnly = true;//true - for skin depth, use only one variations
				//Those May07*.root use the following convention to name histograms
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
			} else if( CalcType.compare(0, 13,"PrimaryHadron") == 0){//Hadron Calculators
//				std::cout<<"CHECK  Find a Hadron"<<std::endl;

				fprimaryHad	=   pset.get< std::vector<int>>("PrimaryHadronGeantCode");//for Feynman Scaling


				if( CalcType == "PrimaryHadronNormalization" ){//k-

					fParameterSet.Sample(engine);//Yes.. load it twice;
					
					//		CHECK, use old method to generate random#;
					//		fParameterSet.Sample(engine);
					//Nothing else to configure

				} else{//Hadron Calculator, slightly complicated

					std::string dataInput       =   pset.get< std::string >("ExternalData");
					std::string ExternalDataInput = sp.find_file(dataInput);
					//std::cout<<__LINE__<<" CHECK find "<<ExternalDataInput<<std::endl;
					TFile* file = new TFile(Form("%s",ExternalDataInput.c_str()));

					std::vector< std::string > pname; // these are what we will extract from the file
					if( CalcType == "PrimaryHadronFeynmanScaling" ){//k+


						pname.push_back("FS/KPlus/FSKPlusFitVal");
						pname.push_back("FS/KPlus/FSKPlusFitCov");

						TArrayD* FSKPlusFitValArray = (TArrayD*) file->Get(pname[0].c_str());
						FitVal = FluxWeightCalc::ConvertToVector(FSKPlusFitValArray);    
						FitCov = (TMatrixD*) file->Get(pname[1].c_str());
						*(FitCov) *= fScalePos*fScalePos;

					}else if( CalcType == "PrimaryHadronSanfordWang" ){//k0
						pname.push_back("SW/K0s/SWK0sFitVal");
						pname.push_back("SW/K0s/SWK0sFitCov");
						TArrayD* SWK0FitValArray = (TArrayD*) file->Get(pname[0].c_str());
						//TArrayD is the most annoying format I have ever experienced so let's convert it to a vector
						FitVal = FluxWeightCalc::ConvertToVector(SWK0FitValArray);    
						FitCov = (TMatrixD*) file->Get(pname[1].c_str());
						*(FitCov) *= fScalePos*fScalePos;


					}else if( CalcType == "PrimaryHadronSWCentralSplineVariation" ){//pi+-

						std::string fitInput = pset.get< std::string >("ExternalFit");
						std::string HadronName; 
						std::string HadronAbriviation;
						if(fprimaryHad[0] == 211){
							HadronName = "PiPlus";
							HadronAbriviation = "PP";
						}
						else if(fprimaryHad[0] == -211){
							HadronName = "PiMinus";
							HadronAbriviation = "PM";
						}
						else{ 
							throw art::Exception(art::errors::StdException)
								<< "sanford-wang is only configured for charged pions ";
						}
						pname.push_back( Form("HARPData/%s/%sCrossSection",HadronName.c_str(),HadronAbriviation.c_str()) ); // Cross Section
						pname.push_back( Form("HARPData/%s/%scovarianceMatrix",HadronName.c_str(),HadronAbriviation.c_str()) ); // Covariance Matrix
						pname.push_back( Form("HARPData/%s/%smomentumBoundsArray",HadronName.c_str(),HadronAbriviation.c_str()) ); // Momentum Bounds
						pname.push_back( Form("HARPData/%s/%sthetaBoundsArray",HadronName.c_str(),HadronAbriviation.c_str()) ); // Theta Bounds

						HARPXSec = (TMatrixD*) file->Get(pname[0].c_str());
						TMatrixD* HARPCov  = (TMatrixD*) file->Get(pname[1].c_str());
						//perform Choleskey Decomposition
						TDecompChol dc = TDecompChol(*(HARPCov));
						if(!dc.Decompose()){
							throw art::Exception(art::errors::StdException)
								<< "Cannot decompose covariance matrix to begin smearing.";
						}
						//Get upper triangular matrix. This maintains the relations in the
						//covariance matrix, but simplifies the structure.
						fIsDecomposed = true;
						FitCov = new TMatrixD(dc.GetU());  
						//HARPLowerTriangluarCov = new TMatrixD(dc.GetU());  


						TArrayD* HARPmomentumBoundsArray = (TArrayD*) file->Get(pname[2].c_str());
						HARPmomentumBounds = FluxWeightCalc::ConvertToVector(HARPmomentumBoundsArray);

						TArrayD* HARPthetaBoundsArray = (TArrayD*) file->Get(pname[3].c_str());
						HARPthetaBounds = FluxWeightCalc::ConvertToVector(HARPthetaBoundsArray);

						/////////////////
						//
						//   Extract the Sanford-Wang Fit Parmeters
						//
						////////////////


						std::string ExternalFitInput = sp.find_file(fitInput);
						TFile* Fitfile = new TFile(Form("%s",ExternalFitInput.c_str()));

						std::string fitname; // these are what we will extract from the file
						fitname = Form("SW/%s/SW%sFitVal",HadronName.c_str(),HadronName.c_str()); // Sanford-Wang Fit Parameters

						TArrayD* SWParamArray = (TArrayD*) Fitfile->Get(fitname.c_str());
						SWParam = FluxWeightCalc::ConvertToVector(SWParamArray);


					}else	validC = false;//the calculator name might have typos



				}//end of special Hadron calculator configurations
			} else	validC = false; //the calculator name is way too off.

			if(validC){//load random numbers based on FitCov
				//{2*multisim vector}
				//{{Ncols() elements}} <-- feed these amount everytime;
				for( int index = 1; index < 2*(FitCov->GetNcols()); index ++){
					fParameterSet.Sample(engine);//load it 2*<number_of_multisims> times
				}
				//					if(validC){
				//						std::cout << "Scale Factor being applied : " << fScalePos << std::endl; 
				//
				//						//Generate 2d Random Numbers here
				//						fWeightArray.resize(2*number_of_multisims);
				//
				//						for (unsigned int i=0;i<fWeightArray.size();i++) {
				//							fWeightArray[i].resize(FitCov->GetNcols());
				//							if(fParameterSet.fRWType == EventWeightParameterSet::kMultisim){
				//								for(unsigned int j = 0; j < fWeightArray[i].size(); j++){
				//									fWeightArray[i][j] = CLHEP::RandGaussQ::shoot(&engine, 0, 1.);
				//								}//Iterate over the covariance matrix size
				//							}
				//						}
				//					}



			}else {
				throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
					<<" calculator "+CalcType + "is invalid"
					<<std::endl;
			}
			std::cout<<"SBNEventWeight : finish configuration."<<std::endl;
		}//End of Configure() function


		std::vector<float> FluxWeightCalc::GetWeight(art::Event& e, size_t inu) {
			std::cout<<"SBNEventWeight : getweight for the "<<inu<<" th particles of an event"<< std::endl;
			//MCFlux & MCTruth
			art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
			e.getByLabel(fGenieModuleLabel,mcFluxHandle);
			std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;

			auto const& mclist = *e.getValidHandle<std::vector<simb::MCTruth>>(fGenieModuleLabel);
			//one line definition above.
			//art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
			//e.getByLabel(fGenieModuleLabel,mctruthHandle);
			//std::vector<simb::MCTruth> const& mclist = *mctruthHandle;

			// No neutrinos in this event gives 0;
			int NUni = fParameterSet.fNuniverses;
			std::vector<float> weights;//( mclist.size(), 0);
			if(mclist.size() == 0){ 
				std::cout<<"EMPTY WEIGHTS"<<std::endl;
				return weights;
			}

			//Iterate through each neutrino in the event
//CHECK no reseting			
//			std::cout<<"CHECK Reset weight counters for "<<CalcType<<std::endl;
//
//			wc = 0;
//			wcn = 0;
//			wc0 = 0;
//			wc1 = 0;
//			wc30 = 0;

//			std::cout<<__LINE__<<" Check type "<<CalcType<<std::endl;
			if( CalcType == "Unisim"){//Unisim Calculator
//			std::cout<<__LINE__<<" Got a Unisim  "<<inu<<" "<<CalcType<<std::endl;
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
						//					weights[i]=UnisimWeightCalc(enu, ptype, ntype, fWeightArray[i], PosOnly);
						double randomN = (fParameterSet.fParameterMap.begin())->second[i];
						weights[i]=UnisimWeightCalc(enu, ptype, ntype, randomN , PosOnly);//AddParameter result does not work here;

						if(weights[i]<0){
							wcn++;
						}else if((weights[i]-0)<1e-30){
							wc0++;
						}else if(fabs(weights[i]-30) < 1e-30){
							wc30++;
						} else if(fabs(weights[i]-1)<1e-30){
							wc1++;
						} else {
							wc++;
						}


						//					if(weights[i]>1) validate_code++;
					}//Iterate through the number of universes      
				}
			} else{//then this must be PrimaryHadron


				// First let's check that the parent of the neutrino we are looking for is 
				//  the particle we intended it to be, if not set all weights to 1
			if (std::find(fprimaryHad.begin(), fprimaryHad.end(),(fluxlist[inu].ftptype)) == fprimaryHad.end() ){//if it does not contain any particles we need get 1
//			std::cout<<__LINE__<<" Check ftptype "<<fluxlist[inu].ftptype<<" and requested "<<fprimaryHad[0]<<std::endl;
					weights.resize( NUni);
					std::fill(weights.begin(), weights.end(), 1);
					std::cout<<"  "<<fluxlist[inu].ftptype<<" not found. All weights are 1 "<<std::endl;
					return weights;//done, all 1
				}// Hadronic parent check

				if(fParameterSet.fRWType == EventWeightParameterSet::kMultisim){
//				std::cout<<"Check 	Multisim weights size "<<weights.size()<<" Universe "<<NUni<<std::endl;

					for (unsigned int i = 0; int(weights.size()) < NUni; i++) {//if all weights are 1, no need to calculate weights;
//						std::cout<<"CHECK Universe: "<<NUni<<std::endl;
						std::pair<bool, double> test_weight;
						
						std::vector< float > Vrandom = (fParameterSet.fParameterMap.begin())->second;//vector of random #
						std::vector< float > subVrandom;//sub-vector of random numbers;
						if( CalcType == "PrimaryHadronNormalization"){//Normalization
							test_weight = PHNWeightCalc(fluxlist[inu], Vrandom[i]);
						} else {
							subVrandom = {Vrandom.begin()+i*FitCov->GetNcols(), Vrandom.begin()+(i+1)*FitCov->GetNcols()};
							if( CalcType == "PrimaryHadronFeynmanScaling"){//FeynmanScaling
								test_weight = PHFSWeightCalc(fluxlist[inu], subVrandom);
							} else if( CalcType == "PrimaryHadronSanfordWang"){//SanfordWang
								test_weight = PHSWWeightCalc(fluxlist[inu], subVrandom);
							} else if( CalcType == "PrimaryHadronSWCentralSplineVariation"){//SWCentaralSplineVariation
								//						std::cout<<__LINE__<<" Get this calculator "<<std::endl;
								test_weight = PHSWCSVWeightCalc(fluxlist[inu], subVrandom);

							} else throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": this shouldnt happen.."<<std::endl;
						}

						if(test_weight.first){	
							weights.push_back(test_weight.second);
							double tmp_weight = test_weight.second;
							std::cout<<i<<" ("<<tmp_weight<<") ";
							if(tmp_weight<0){
								wcn++;
							}else if((tmp_weight-0)<1e-30){
								wc0++;
							}else if(fabs(tmp_weight-30) < 1e-30){
								wc30++;
							} else if(fabs(tmp_weight-1)<1e-30){
								wc1++;
							} else {
								wc++;
							}
						
						};

					}//Iterate through the number of universes      
				}//Yes, Multisim
			}


			std::cout<<"Weights counter: (normal, <0, ==0, 1, 30)= ";
			std::cout<<wc<<", "<<wcn<<", "<<wc0<<", "<<wc1<<", "<<wc30<<std::endl;
			std::cout<<"Next partile/event\n"<<std::endl;

			return weights;
		}//GetWeight()


		std::vector<double> FluxWeightCalc::ConvertToVector(TArrayD const* array) {
			std::vector<double> v(array->GetSize());
			std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
					v.begin());
			return v;
		} // ConvertToVector()

		REGISTER_WEIGHTCALC(FluxWeightCalc)
	}  // namespace evwgh
}  // namespace sbn
