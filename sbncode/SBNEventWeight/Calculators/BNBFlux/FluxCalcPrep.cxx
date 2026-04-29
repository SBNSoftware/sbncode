#include "FluxCalcPrep.h"
//Need to include headers explicitly after v9_33_00

#include "fhiclcpp/ParameterSet.h"


namespace sbn {
  namespace evwgh {

    //Configure everything here! 
    void FluxWeightCalc::Configure(fhicl::ParameterSet const& p,
				   CLHEP::HepRandomEngine& engine) {
      std::cout<<"SBNEventWeight Flux: configure "<< GetName() << std::endl;

      fGeneratorModuleLabel = p.get<std::string>("generator_module_label");//use this label to get MC*Handle
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
        //Note: no sigma is used in this script.
        fParameterSet.AddParameter(pars[i], parsigmas[i]);//parsigmas[i]);
      }


      //1. << Calculator Settings >>
      CalcType = pset.get< std::string >("calc_type");//Unisim,PrimaryHadronSWCentralSplineVariation,PrimaryHadronFeynmanScaling,PrimaryHadronSanfordWang,PrimaryHadronNormalization
      std::cout<<"SBNEventWeight Flux: Calculator type: "<<CalcType<<std::endl;

      // For FluxHist
      fUseFluxHist = pset.get<bool>("use_flux_hist", false);

      fScalePos   = pset.get<double>("scale_factor_pos");
      if(!pset.get_if_present("scale_factor_neg",fScaleNeg)){
	std::cout<<"SBNEventWeight Flux: auto-assignment: scale_factor_neg = 1."<<std::endl;
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

      }

      // G4BNB to BooNE


      else if (CalcType == "FluxHist") {
	
        // For debugging you can keep hardcoded paths for now.
        // Later you can switch back to FHiCL-controlled file names.
	//	std::string cv_file =
          //"/exp/sbnd/app/users/rohanr/larsoft_v10_14_02/srcs/sbncode/sbncode/SBNEventWeight/jobs/flux/nuMom_immParent_merged.root";
	//	  "/pnfs/sbnd/scratch/users/rohanr/beamweight/final_hists/nuMom_immParent_merged_200.root";
	//	std::string rw_file =
          //"/exp/sbnd/app/users/rohanr/larsoft_v10_14_02/srcs/sbncode/sbncode/SBNEventWeight/jobs/flux/nuMom_immParent_oldFiles_correctedName.root";
	//	  "/pnfs/sbnd/scratch/users/rohanr/beamweight/final_hists/nuMom_immParent_oldFiles_200_merged.root";

	std::string dataInputCV = pset.get<std::string>("cv_hist_file");
	std::string dataInputRW = pset.get<std::string>("rw_hist_file");

	cet::search_path sp("FW_SEARCH_PATH");

	std::string cv_file = sp.find_file(dataInputCV);
	std::string rw_file = sp.find_file(dataInputRW);

        TFile fcv(cv_file.c_str());
        TFile frw(rw_file.c_str());

        if (fcv.IsZombie() || frw.IsZombie()) {
          throw cet::exception("FluxHist")
            << "Could not open FluxHist ROOT file(s):\n"
            << "  CV = " << cv_file << "\n"
            << "  RW = " << rw_file << "\n";
        }

	std::vector<std::tuple<int,std::string>> validCombinations = {
	  { 211,  "numu"    },
	  { 211,  "nue"     },
	  {-211,  "numubar" },
	  {-211,  "nuebar"  },
	  { 321,  "numu"    },
	  { 321,  "nue"     },
	  {-321,  "numubar" },
	  {-321,  "nuebar"  },
	  { 130,  "numu"     },
	  { 130,  "numubar"  },
	  { 130,  "nue"     },
	  { 130,  "nuebar"  },
	  {  13,  "numu"    },
	  {  13,  "nuebar"  },
	  { -13,  "numubar" },
	  { -13,  "nue"     }
	};

	for (size_t i = 0; i < validCombinations.size(); ++i) {
	  int pdg = std::get<0>(validCombinations[i]);
	  std::string nut = std::get<1>(validCombinations[i]);

	  std::string hname = Form("hSBND_%d_%s_3mom", pdg, nut.c_str());

	  TH3D* hcv = dynamic_cast<TH3D*>(fcv.Get(hname.c_str()));
	  TH3D* hrw = dynamic_cast<TH3D*>(frw.Get(hname.c_str()));

	  if (!hcv || !hrw) {
	    throw cet::exception("FluxHist")
	      << "Missing TH3D histogram: " << hname << "\n";
	  }

	  fCVHist3D[i] = dynamic_cast<TH3D*>(hcv->Clone(Form("%s_cv_clone", hname.c_str())));
	  fRWHist3D[i] = dynamic_cast<TH3D*>(hrw->Clone(Form("%s_rw_clone", hname.c_str())));

	  if (!fCVHist3D[i] || !fRWHist3D[i]) {
	    throw cet::exception("FluxHist")
	      << "Failed to clone TH3D histogram: " << hname << "\n";
	  }

	  fCVHist3D[i]->SetDirectory(nullptr);
	  fRWHist3D[i]->SetDirectory(nullptr);

	  std::cout << "Loaded 3D histogram: " << hname
		    << " CV entries=" << fCVHist3D[i]->GetEntries()
		    << " RW entries=" << fRWHist3D[i]->GetEntries()
		    << std::endl;
	}

	/*
	std::string compNames[3] = {"px", "py", "pz"};

        for (size_t i = 0; i < validCombinations.size(); ++i) {
          int pdg = std::get<0>(validCombinations[i]);
	  std::string nut = std::get<1>(validCombinations[i]);

	  for (int icomp = 0; icomp < 3; ++icomp) {
	    std::string hname = Form("hSBND_%d_%s_%s", pdg, nut.c_str(), compNames[icomp].c_str());

	    TH1* hcv = dynamic_cast<TH1*>(fcv.Get(hname.c_str()));
	    TH1* hrw = dynamic_cast<TH1*>(frw.Get(hname.c_str()));

	    if (!hcv || !hrw) {
	      throw cet::exception("FluxHist")
		<< "Missing histogram: " << hname << "\n";
	    }

	  // std::cout << "Loaded histogram: " << hname
	     << "  CV entries=" << hcv->GetEntries()
		       << "  RW entries=" << hrw->GetEntries()
		       << std::endl;
	    //
	    if (hcv->GetNbinsX() < 200 || hrw->GetNbinsX() < 200) {
	      throw cet::exception("FluxHist")
		<< "Histogram " << hname
              << " has fewer than 200 bins.\n";
	    }

	    fCVHist[i][icomp] =
              dynamic_cast<TH1*>(hcv->Clone(Form("%s_cv_clone", hname.c_str())));
            fRWHist[i][icomp] =
              dynamic_cast<TH1*>(hrw->Clone(Form("%s_rw_clone", hname.c_str())));

            if (!fCVHist[i][icomp] || !fRWHist[i][icomp]) {
              throw cet::exception("FluxHist")
                << "Failed to clone histogram: " << hname << "\n";
            }

            fCVHist[i][icomp]->SetDirectory(nullptr);
            fRWHist[i][icomp]->SetDirectory(nullptr);

	    std::cout << "Loaded histogram: " << hname
                      << "  CV entries=" << fCVHist[i][icomp]->GetEntries()
                      << "  RW entries=" << fRWHist[i][icomp]->GetEntries()
                      << std::endl;

	  }
        }
*/
        fcv.Close();
        frw.Close();
      }


      //-------------
      //-- Hadrons --
      //-------------
      else if( CalcType.compare(0, 13,"PrimaryHadron") == 0){//Hadron Calculators

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
      //      std::cout<<"SBNEventWeight : finish configuration."<<std::endl;
    }//End of Configure() function


    std::vector<float> FluxWeightCalc::GetWeight(art::Event& e, size_t inu) {
      bool count_weights = false;
      //      std::cout<<"SBNEventWeight : getweight for the "<<inu<<" th particles of an event"<< std::endl;
      //MCFlux & MCTruth
      art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
      e.getByLabel(fGeneratorModuleLabel, mcFluxHandle);
      std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;

      art::Handle< std::vector<bsim::Dk2Nu> > dk2nuHandle;
      std::vector<bsim::Dk2Nu> const* dk2nu_v = nullptr;
      e.getByLabel(fGeneratorModuleLabel, dk2nuHandle);
      if (dk2nuHandle.isValid()) {
        dk2nu_v = dk2nuHandle.product(); 
      }

      //or do the above 3 lines in one line
      auto const& mclist = *e.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorModuleLabel);

      // If no neutrinos in this event, gives 0 weight;
      int NUni = fParameterSet.fNuniverses;
      std::vector<float> weights;//( mclist.size(), 0);
      if(fluxlist.size() == 0){ 
	std::cout<<"SBNEventWeight Flux: EMPTY WEIGHTS"<<std::endl;
        return weights;
      }

      //Iterate through each neutrino in the event
      std::cout<<"SBNEventWeight Flux: Get calculator "<<CalcType<<std::endl;
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

        // Collect neutrino energy; mclist is replaced with fluxlist.
	//        double enu= mclist[inu].GetNeutrino().Nu().E();
        double enu= fluxlist[inu].fnenergyn;

        //Let's make a weights based on the calculator you have requested 

        if(fParameterSet.fRWType == EventWeightParameterSet::kMultiSim){

          for (size_t i=0;i<weights.size();i++) {
            double randomN = (fParameterSet.fParameterMap.begin())->second[i];
            weights[i]=UnisimWeightCalc(enu, ptype, ntype, randomN , PosOnly);//AddParameter result does not work here;

	    if(count_weights){
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
	    }
          }//Iterate through the number of universes      
        }
      }

      else if (CalcType == "FluxHist") {

	weights.resize(NUni, 1.0);

	if (inu >= fluxlist.size()) {
	  throw cet::exception("FluxHist")
	    << "inu out of range for fluxlist: inu=" << inu
	    << " fluxlist.size()=" << fluxlist.size() << "\n";
	}

	if (inu >= mclist.size()) {
	  throw cet::exception("FluxHist")
	    << "inu out of range for mclist: inu=" << inu
	    << " mclist.size()=" << mclist.size() << "\n";
	}

	int parent_pdg = fluxlist[inu].fptype;
	int nu_pdg     = fluxlist[inu].fntype;

	double px = mclist[inu].GetNeutrino().Nu().Px();
	double py = mclist[inu].GetNeutrino().Nu().Py();
	double pz = mclist[inu].GetNeutrino().Nu().Pz();

	double w = FluxHistWeightCalc(parent_pdg, nu_pdg, px, py, pz);

	for (int i = 0; i < NUni; ++i) {
	  weights[i] = std::isfinite(w) ? w : 1.0;
	}
      }

      /*
      else if (CalcType == "FluxHist") {

        weights.resize(NUni, 1.0);

        if (inu >= fluxlist.size()) {
          throw cet::exception("FluxHist")
            << "inu out of range for fluxlist: inu=" << inu
            << " fluxlist.size()=" << fluxlist.size() << "\n";
        }

        if (inu >= mclist.size()) {
          throw cet::exception("FluxHist")
            << "inu out of range for mclist: inu=" << inu
            << " mclist.size()=" << mclist.size() << "\n";
        }

        int parent_pdg = fluxlist[inu].fptype;
        int nu_pdg     = fluxlist[inu].fntype;

        int idx = GetFluxHistIndex(parent_pdg, nu_pdg);
        if (idx < 0 || idx >= 14) {
          return weights;
        }

        double px = mclist[inu].GetNeutrino().Nu().Px();
        double py = mclist[inu].GetNeutrino().Nu().Py();
        double pz = mclist[inu].GetNeutrino().Nu().Pz();

        double w = FluxHistWeightCalc(parent_pdg, nu_pdg, px, py, pz);

        // Deterministic event-level correction: same value in every universe slot
        for (int i = 0; i < NUni; ++i) {
          weights[i] = std::isfinite(w) ? w : 1.0;
        }
      }

      */
      /*
      else if (CalcType == "FluxHist") {

	std::cout << "Entering FluxHist GetWeight for inu = " << inu << std::endl;

        weights.resize(NUni, 1.0);

        if (inu >= fluxlist.size()) {
          throw cet::exception("FluxHist")
            << "inu out of range for fluxlist: inu=" << inu
            << " fluxlist.size()=" << fluxlist.size() << "\n";
        }

        if (inu >= mclist.size()) {
          throw cet::exception("FluxHist")
            << "inu out of range for mclist: inu=" << inu
            << " mclist.size()=" << mclist.size() << "\n";
        }

        int parent_pdg = fluxlist[inu].fptype;
        int nu_pdg     = fluxlist[inu].fntype;


        int idx = GetFluxHistIndex(parent_pdg, nu_pdg);

	//std::cout << "  FluxHist index = " << idx << std::endl;

        if (idx < 0 || idx >= 14) {
	  std::cout << "  Unsupported parent/flavor combination, returning unit weights"
                    << std::endl;
          return weights;
        }

        double pz = mclist[inu].GetNeutrino().Nu().Pz();

	std::cout << "  neutrino pz = " << pz << std::endl;

        int bin = static_cast<int>(pz / 0.05);

        if (bin < 0)   bin = 0;
        if (bin > 199) bin = 199;

	std::cout << "  histogram bin = " << bin << std::endl;

        if (fParameterSet.fRWType == EventWeightParameterSet::kMultiSim) {

          auto const& randomVec = (fParameterSet.fParameterMap.begin())->second;

	  std::cout << "  NUni = " << NUni
                    << ", randomVec.size() = " << randomVec.size() << std::endl;

          for (int i = 0; i < NUni; ++i) {

            if (i >= (int)randomVec.size()) {
              throw cet::exception("FluxHist")
                << "Random vector too short: i=" << i
                << " size=" << randomVec.size() << "\n";
            }

	    //            double randomN = randomVec[i];

            double cv = fCVHist[idx][bin];
            double rw = fRWHist[idx][bin];


            double w = 1.0;

            if (cv > 0.0) {
              double ratio = rw / cv;
              if (std::isfinite(ratio)) {
                w = ratio;//1.0 - (1.0 - ratio) * randomN;
              }
            }

            weights[i] = std::isfinite(w) ? w : 1.0;
          }
        }

	std::cout << "Leaving FluxHist GetWeight" << std::endl;
      }

      */
      else{//then this must be PrimaryHadron


        // First let's check that the parent of the neutrino we are looking for is 
        //  the particle we intended it to be, if not set all weights to 1

	simb::MCFlux flux;
	flux.ftptype = fluxlist[inu].ftptype;
	flux.ftpx = fluxlist[inu].ftpx;
	flux.ftpy = fluxlist[inu].ftpy;
	flux.ftpz = fluxlist[inu].ftpz;

	// If Dk2Nu flux, use ancestors to evaluate tptype
	if (fluxlist[inu].fFluxType == simb::kDk2Nu) {

	  for( const bsim::Ancestor & ancestor : dk2nu_v->at(inu).ancestor ) {
	    std::string aproc = ancestor.proc;
	    if( (aproc.find("BooNEHadronInelastic:BooNEpBeInteraction") != std::string::npos) && (aproc.find("QEBooNE") == std::string::npos) ) {
	      flux.ftptype = ancestor.pdg;
	      flux.ftpx = ancestor.startpx;
	      flux.ftpy = ancestor.startpy;
	      flux.ftpz = ancestor.startpz;
	    } // found it
	  } // find first inelastic
	}

	if (std::find(fprimaryHad.begin(), fprimaryHad.end(),(flux.ftptype)) == fprimaryHad.end() ){//if it does not contain any particles we need get 1
          weights.resize( NUni);
	  std::fill(weights.begin(), weights.end(), 1);
          return weights;//done, all 1
        }// Hadronic parent check

        if(fParameterSet.fRWType == EventWeightParameterSet::kMultiSim){

          for (unsigned int i = 0; int(weights.size()) < NUni; i++) {//if all weights are 1, no need to calculate weights;
	    std::pair<bool, double> test_weight;

	    std::vector< float > Vrandom = (fParameterSet.fParameterMap.begin())->second;//vector of random #
	    std::vector< float > subVrandom;//sub-vector of random numbers;
            if( CalcType == "PrimaryHadronNormalization"){//Normalization
              test_weight = PHNWeightCalc(flux, Vrandom[i]);

            } else {
              subVrandom = {Vrandom.begin()+i*FitCov->GetNcols(), Vrandom.begin()+(i+1)*FitCov->GetNcols()};
              if( CalcType == "PrimaryHadronFeynmanScaling"){//FeynmanScaling
                test_weight = PHFSWeightCalc(flux, subVrandom);

              } else if( CalcType == "PrimaryHadronSanfordWang"){//SanfordWang
                test_weight = PHSWWeightCalc(flux, subVrandom);

              } else if( CalcType == "PrimaryHadronSWCentralSplineVariation"){//SWCentaralSplineVariation
                test_weight = PHSWCSVWeightCalc(flux, subVrandom);
              } else throw cet::exception(__PRETTY_FUNCTION__) << GetName() << ": this shouldnt happen.."<<std::endl;
            }

            if(test_weight.first){  
              weights.push_back(test_weight.second);
              double tmp_weight = test_weight.second;
	      //              std::cout<<i<<" ("<<tmp_weight<<") ";
	      if(count_weights){
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
	      }

            };

          }//Iterate through the number of universes      
        }//Yes, MultiSim
      }

      if(count_weights){
	std::cout<<"SBNEventWeight Flux: Weights counter: (normal, <0, ==0, 1, 30)= ";
	std::cout<<wc<<", "<<wcn<<", "<<wc0<<", "<<wc1<<", "<<wc30<<std::endl;
      }
//      std::cout<<"Next partile/event\n"<<std::endl;

      return weights;
    }//GetWeight()


    std::vector<double> FluxWeightCalc::ConvertToVector(TArrayD const* array) {
      std::vector<double> v(array->GetSize());
      std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
          v.begin());
      return v;
    } // ConvertToVector()

    double FluxWeightCalc::FluxHistWeightCalc(int parent_pdg, int nu_pdg,
					      double px, double py, double pz)
    {
      static int debug_counter = 0;

      int idx = GetFluxHistIndex(parent_pdg, nu_pdg);
      if (idx < 0) return 1.0;

      TH3D* hcv = fCVHist3D[idx];
      TH3D* hrw = fRWHist3D[idx];

      if (!hcv || !hrw) return 1.0;

      int binx = hcv->GetXaxis()->FindBin(px);
      int biny = hcv->GetYaxis()->FindBin(py);
      int binz = hcv->GetZaxis()->FindBin(pz);

      if (binx < 1) binx = 1;
      if (binx > hcv->GetNbinsX()) binx = hcv->GetNbinsX();

      if (biny < 1) biny = 1;
      if (biny > hcv->GetNbinsY()) biny = hcv->GetNbinsY();

      if (binz < 1) binz = 1;
      if (binz > hcv->GetNbinsZ()) binz = hcv->GetNbinsZ();

      double cv = hcv->GetBinContent(binx, biny, binz);
      double rw = hrw->GetBinContent(binx, biny, binz);

      double weight = 1.0;

      if (cv > 0.0) {
	weight = rw / cv;
	if (!std::isfinite(weight) || weight <= 0.0) {
	  weight = 1.0;
	}
      }

      if (debug_counter < 20) {
	std::cout << "\n[FluxHist 3D DEBUG]"
		  << " parent=" << parent_pdg
		  << " nu=" << nu_pdg
		  << "\n  px=" << px << " binx=" << binx
		  << "\n  py=" << py << " biny=" << biny
		  << "\n  pz=" << pz << " binz=" << binz
		  << "\n  CV=" << cv
		  << "\n  RW=" << rw
		  << "\n  weight=" << weight
		  << std::endl;
      }

      debug_counter++;

      return std::isfinite(weight) ? weight : 1.0;
    }
    /*
    double FluxWeightCalc::FluxHistWeightCalc(int parent_pdg, int nu_pdg,
                                              double px, double py, double pz)
    {
      static int debug_counter = 0;

      int idx = GetFluxHistIndex(parent_pdg, nu_pdg);
      if (idx < 0) return 1.0;

      double comps[3] = {px, py, pz};
      const char* compNames[3] = {"px", "py", "pz"};

      double log_sum = 0.0;
      int n_used = 0;

      if (debug_counter < 20) {
	std::cout << "\n[FluxHist DEBUG] Event " << debug_counter
                  << " parent=" << parent_pdg
                  << " nu=" << nu_pdg << std::endl;
      }

      for (int icomp = 0; icomp < 3; ++icomp) {
        TH1* hcv = fCVHist[idx][icomp];
        TH1* hrw = fRWHist[idx][icomp];

        if (!hcv || !hrw) continue;

        double val = comps[icomp];
        int bin = hcv->GetXaxis()->FindBin(val);

        if (bin < 1) bin = 1;
        if (bin > hcv->GetNbinsX()) bin = hcv->GetNbinsX();

        double cv = hcv->GetBinContent(bin);
        double rw = hrw->GetBinContent(bin);

        double ratio = 1.0;

        if (cv > 0.0) {
          ratio = rw / cv;
          if (!std::isfinite(ratio) || ratio <= 0.0) ratio = 1.0;
        }


        log_sum += std::log(ratio);
        n_used++;

        if (debug_counter < 20) {
	  std::cout << "  " << compNames[icomp]
                    << " = " << val
                    << "  bin=" << bin
                    << "  CV=" << cv
                    << "  RW=" << rw
                    << "  ratio=" << ratio
                    << std::endl;
        }
      }

      if (n_used == 0) return 1.0;

      double weight = std::exp(log_sum / n_used);

      if (debug_counter < 20) {
	std::cout << "  --> FINAL WEIGHT = " << weight << std::endl;
      }

      debug_counter++;

      return std::isfinite(weight) ? weight : 1.0;
    }
    */

    REGISTER_WEIGHTCALC(FluxWeightCalc)
  }  // namespace evwgh
}  // namespace sbn
