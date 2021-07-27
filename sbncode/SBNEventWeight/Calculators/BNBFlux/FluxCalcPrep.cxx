#include "FluxCalcPrep.h"


namespace sbn {
  namespace evwgh {

    //Configure everything here! 
    void FluxWeightCalc::Configure(fhicl::ParameterSet const& p,
        CLHEP::HepRandomEngine& engine) {
      std::cout<<"SBNEventWeight Flux: configure "<< GetName() << std::endl;

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
        std::cout<<"SBNEventWeight Flux: `parameter_sigma` was not set; now it is set as 1"<<std::endl;
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
      std::cout<<"SBNEventWeight Flux: Calculator type: "<<CalcType<<std::endl;

      fScalePos   = pset.get<double>("scale_factor_pos");
      if(!pset.get_if_present("scale_factor_neg",fScaleNeg)){
        std::cout<<"SBNEventWeight Flux: auto-asignment: scale_factor_neg = 1."<<std::endl;
        }

      cet::search_path sp("FW_SEARCH_PATH");
      bool validC = true;
      
      //----------------------------
      //-- Non hadrons production --
      //----------------------------
      

      fParameterSet.Sample(engine);//random_seed is loaded at sbncode/Base/WeightManager.h

      if( CalcType == "Unisim"){//Unisim Calculator

        std::string dataInput1  = pset.get< std::string >("CentralValue_hist_file");
        std::string cvfile_path  = sp.find_file(dataInput1);
        TFile fcv(Form("%s",cvfile_path.c_str()));

        std::string dataInput2pos  = pset.get< std::string >("PositiveSystematicVariation_hist_file");
        std::string rwfilepos    = sp.find_file(dataInput2pos);
        TFile frwpos(Form("%s", rwfilepos.c_str()));

        std::string dataInput2neg  = pset.get< std::string >("NegativeSystematicVariation_hist_file");
        std::string rwfileneg    = sp.find_file(dataInput2neg);
        TFile frwneg(Form("%s", rwfileneg.c_str()));

        if(dataInput2pos == dataInput2neg) PosOnly = true;//true - for skin depth
        //Those May07*.root use the following convention to name histograms
        int cptype[4] = {1,2,3,4}; //mu, pi, k0, k
        int cntype[4] = {1,2,3,4}; //nue, anue, numu, anumu

        for (int iptyp=0;iptyp<4;iptyp++) {
          for (int intyp=0;intyp<4;intyp++) {
            for (int ibin=0;ibin<200;ibin++) { //Grab events from ibin+1 
              fCV[iptyp][intyp][ibin]=(dynamic_cast<TH1F*>  (fcv.Get(Form("h5%d%d",cptype[iptyp],cntype[intyp]))))->GetBinContent(ibin+1);
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

        fprimaryHad  =   pset.get< std::vector<int>>("PrimaryHadronGeantCode");

        if( CalcType == "PrimaryHadronNormalization" ){//k-

          fParameterSet.Sample(engine);//Yes.. load it again;

        } else{//Other Hadron Calculators

          std::string dataInput       =   pset.get< std::string >("ExternalData");
          std::string ExternalDataInput = sp.find_file(dataInput);
          TFile* file = new TFile(Form("%s",ExternalDataInput.c_str()));

          std::vector< std::string > pname; // these are keys to histograms
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

            FitVal = FluxWeightCalc::ConvertToVector(SWK0FitValArray);//TArrayD--> vector
            FitCov = (TMatrixD*) file->Get(pname[1].c_str());
            *(FitCov) *= fScalePos*fScalePos;

          }else if( CalcType == "PrimaryHadronSWCentralSplineVariation" ){//pi+-

            std::string fitInput = pset.get< std::string >("ExternalFit");
            std::string HadronName; 
            std::string HadronAbriviation;
            if(fprimaryHad[0] == 211){
              HadronName = "PiPlus";
              HadronAbriviation = "PP";
            } else if(fprimaryHad[0] == -211){
              HadronName = "PiMinus";
              HadronAbriviation = "PM";
            } else{ 
              throw art::Exception(art::errors::StdException)
                << "sanford-wang is only configured for charged pions ";
            }
            pname.push_back( Form("HARPData/%s/%sCrossSection",HadronName.c_str(),HadronAbriviation.c_str()) ); // Cross Section
            pname.push_back( Form("HARPData/%s/%scovarianceMatrix",HadronName.c_str(),HadronAbriviation.c_str()) ); // Covariance Matrix
            pname.push_back( Form("HARPData/%s/%smomentumBoundsArray",HadronName.c_str(),HadronAbriviation.c_str()) ); // Momentum Bounds
            pname.push_back( Form("HARPData/%s/%sthetaBoundsArray",HadronName.c_str(),HadronAbriviation.c_str()) ); // Theta Bounds

            HARPXSec = (TMatrixD*) file->Get(pname[0].c_str());
            TMatrixD* HARPCov  = (TMatrixD*) file->Get(pname[1].c_str());

            TDecompChol dc = TDecompChol(*(HARPCov));//perform Choleskey Decomposition
            if(!dc.Decompose()){
              throw art::Exception(art::errors::StdException)
                << "Cannot decompose covariance matrix to begin smearing.";
            }

            //Get upper triangular matrix. This maintains the relations in the
            //  covariance matrix, but simplifies the structure.
            fIsDecomposed = true;
            FitCov = new TMatrixD(dc.GetU());  


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


          }else  validC = false;//slightly incorrect calculator name in *fcl

          if(validC){//load random numbers based on FitCov
            //{2*multisim vector}
            //{{Ncols() elements}} <-- feed these amount everytime;
            for( int index = 1; index < 2*(FitCov->GetNcols()); index ++){
              fParameterSet.Sample(engine);//load it 2*<number_of_multisims> times
            }
  
          }else {
            throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": "
              <<" calculator "+CalcType + "is invalid"
              <<std::endl;
          }


        }//end of special Hadron calculator configurations

      } else  validC = false; //the calculator name is way too off.
      std::cout<<"SBNEventWeight : finish configuration."<<std::endl;
    }//End of Configure() function


    std::vector<float> FluxWeightCalc::GetWeight(art::Event& e, size_t inu) {
//      std::cout<<"SBNEventWeight : getweight for the "<<inu<<" th particles of an event"<< std::endl;
      //MCFlux & MCTruth
      art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
      e.getByLabel(fGenieModuleLabel,mcFluxHandle);
      std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;
      //or do the above 3 lines in one line
      auto const& mclist = *e.getValidHandle<std::vector<simb::MCTruth>>(fGenieModuleLabel);

      // If mo neutrinos in this event, gives 0 weight;
      int NUni = fParameterSet.fNuniverses;
      std::vector<float> weights;//( mclist.size(), 0);
      if(mclist.size() == 0){ 
        std::cout<<"SBNEventWeight Flux: EMPTY WEIGHTS"<<std::endl;
        return weights;
      }

      //Iterate through each neutrino in the event
      if( CalcType == "Unisim"){//Unisim Calculator
        weights.resize(NUni);
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

          for (size_t i=0;i<weights.size();i++) {
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
          }//Iterate through the number of universes      
        }
      } else{//then this must be PrimaryHadron


        // First let's check that the parent of the neutrino we are looking for is 
        //  the particle we intended it to be, if not set all weights to 1
      if (std::find(fprimaryHad.begin(), fprimaryHad.end(),(fluxlist[inu].ftptype)) == fprimaryHad.end() ){//if it does not contain any particles we need get 1
          weights.resize( NUni);
          std::fill(weights.begin(), weights.end(), 1);
          return weights;//done, all 1
        }// Hadronic parent check

        if(fParameterSet.fRWType == EventWeightParameterSet::kMultisim){

          for (unsigned int i = 0; int(weights.size()) < NUni; i++) {//if all weights are 1, no need to calculate weights;
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
                test_weight = PHSWCSVWeightCalc(fluxlist[inu], subVrandom);

              } else throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": this shouldnt happen.."<<std::endl;
            }

            if(test_weight.first){  
              weights.push_back(test_weight.second);
              double tmp_weight = test_weight.second;
//              std::cout<<i<<" ("<<tmp_weight<<") ";
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

      std::cout<<"SBNEventWeight Flux: Weights counter: (normal, <0, ==0, 1, 30)= ";
      std::cout<<wc<<", "<<wcn<<", "<<wc0<<", "<<wc1<<", "<<wc30<<std::endl;
//      std::cout<<"Next partile/event\n"<<std::endl;

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

