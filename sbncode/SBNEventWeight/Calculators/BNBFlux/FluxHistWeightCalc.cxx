#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"

#include <iostream>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TFile.h"
#include "TH1F.h"

namespace evwgh {
  class FluxHistWeightCalc : public WeightCalc
  {
  public:
    FluxHistWeightCalc() = default;
    void Configure(fhicl::ParameterSet const& pset,
                   CLHEP::HepRandomEngine& engine);
    std::vector<std::vector<double> > GetWeight(art::Event & e);
    
  private:    
    std::vector<double> fWeightArray{};
    int fNmultisims{};
    std::string fMode{};
    std::string fGenieModuleLabel{};

    //         pi+-,k+-,k0,mu+- 
    //         |  numu, numubar, nue, nuebar 
    //         |  |   50MeV bins
    //         |  |   |
    double fCV[7][4][200];
    double fRW[7][4][200];
    
    DECLARE_WEIGHTCALC(FluxHistWeightCalc)
  };

  void FluxHistWeightCalc::Configure(fhicl::ParameterSet const& p,
                                     CLHEP::HepRandomEngine& engine)
  {    
    //global config
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");

    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    //calc config
    fNmultisims = pset.get<int>("number_of_multisims");
    fMode       = pset.get<std::string>("mode");
    std::string dataInput1 = pset.get< std::string >("cv_hist_file");
    std::string dataInput2 = pset.get< std::string >("rw_hist_file");

    cet::search_path sp("FW_SEARCH_PATH");
    std::string cvfile = sp.find_file(dataInput1);
    std::string rwfile = sp.find_file(dataInput2);
    
    std::string ptype[] = {"pi+", "pi-", "k+", "k-", "k0", "mu-", "mu+"};
    std::string ntype[] = {"numu", "numubar", "nue", "nuebar"};

    TFile fcv(Form("%s",cvfile.c_str()));
    TFile frw(Form("%s",rwfile.c_str()));
    for (int iptyp=0;iptyp<7;iptyp++) {
      for (int intyp=0;intyp<4;intyp++) {
	for (int ibin=0;ibin<200;ibin++) {
	  fCV[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (fcv.Get(Form("h_%s_%s",ptype[iptyp].c_str(),ntype[intyp].c_str()))))->GetBinContent(ibin+1);
	  fRW[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (frw.Get(Form("h_%s_%s",ptype[iptyp].c_str(),ntype[intyp].c_str()))))->GetBinContent(ibin+1);
	}
      }
    }
    fcv.Close();
    frw.Close();

    fWeightArray.resize(fNmultisims);

    if (fMode.find("multisim") != std::string::npos )
      for (double& weight : fWeightArray) weight = CLHEP::RandGaussQ::shoot(&engine, 0, 1.);
    else
      for (double& weight : fWeightArray) weight = 1.;
  }

  std::vector<std::vector<double> > FluxHistWeightCalc::GetWeight(art::Event & e)
  {
    //calculate weight(s) here 
    std::vector<std::vector<double> > weight;

    // * MC flux information
    auto mcfluxListHandle = e.getHandle< std::vector<simb::MCFlux> >(fGenieModuleLabel);
    if (!mcfluxListHandle) {
      return weight;
    }

    // * MC truth information
    auto mctruthListHandle = e.getHandle< std::vector<simb::MCTruth> >(fGenieModuleLabel);
    if (!mctruthListHandle) {
      return weight;
    }

    std::vector<simb::MCFlux> const& fluxlist = *mcfluxListHandle;
    std::vector<simb::MCTruth> const& mclist = *mctruthListHandle;

    weight.resize(mclist.size());
    for (unsigned int inu=0;inu<mclist.size();inu++) {
      weight[inu].resize(fNmultisims);
     
      int ptype=-9999;
      int ntype=-9999;
      int bin=-9999;
      
      if ( fluxlist[inu].fptype==211 ) ptype = 0;
      else if ( fluxlist[inu].fptype==-211 ) ptype = 1;
      else if ( fluxlist[inu].fptype==321 ) ptype = 2;
      else if ( fluxlist[inu].fptype==-321 ) ptype = 3;
      else if ( fluxlist[inu].fptype==130 ) ptype = 4;
      else if ( fluxlist[inu].fptype==13 ) ptype = 5;
      else if ( fluxlist[inu].fptype==-13 ) ptype = 6;
      else {
        throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<fluxlist[0].fptype<< std::endl;
      }
      
      if ( fluxlist[inu].fntype==14 ) ntype=0;
      else if ( fluxlist[inu].fntype==-14 ) ntype=1;
      else if ( fluxlist[inu].fntype==12 ) ntype=2;
      else if ( fluxlist[inu].fntype==-12 ) ntype=3;
      else {
        throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<fluxlist[0].fptype<< std::endl;
      }
      
      double enu=mclist[inu].GetNeutrino().Nu().E();
      bin=enu/0.05;     
      for (int i=0;i<fNmultisims;i++) {
	double test = 1-(1-fRW[ptype][ntype][bin]/fCV[ptype][ntype][bin])*fWeightArray[i];
	
	// Guards against inifinite weights
	if(std::isfinite(test)){ weight[inu][i] = test;}
	else{weight[inu][i] = 1;}

      }
    }
    return weight;
  }
  REGISTER_WEIGHTCALC(evwgh::FluxHistWeightCalc)
}
//REGISTER_WEIGHTCALC(FluxHistWeightCalc)
