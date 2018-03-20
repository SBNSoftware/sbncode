#include "larsim/EventWeight/Base/WeightCalc.h"

#include <iostream>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
//#include "artextensions/SeedService/SeedService.hh"
#include "nutools/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "nusimdata/SimulationBase/MCFlux.h"
//#include "SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
//#include "SimulationBase/MCTruth.h"

#include "TFile.h"
#include "TH1F.h"

namespace sbncode {
namespace evwgh {
  class FluxHistWeightCalc : public ::evwgh::WeightCalc
  {
  public:
    FluxHistWeightCalc();
    void Configure(fhicl::ParameterSet const& pset);
    std::vector<std::vector<double> > GetWeight(art::Event & e);
    
  private:    
    CLHEP::RandGaussQ *fGaussRandom;
    std::vector<double> fWeightArray;
    int fNmultisims;
    std::string fMode;
    std::string fGenieModuleLabel;

    //         pi+-,k+-,k0,mu+- 
    //         |  numu, numubar, nue, nuebar 
    //         |  |   50MeV bins
    //         |  |   |
    double fCV[4][4][200];
    double fRW[4][4][200];
    
    DECLARE_WEIGHTCALC(FluxHistWeightCalc)
  };
  FluxHistWeightCalc::FluxHistWeightCalc()
  {
  }

  void FluxHistWeightCalc::Configure(fhicl::ParameterSet const& p)
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
    
    std::string ptype[] = {"pi", "k", "k0", "mu"};
    std::string ntype[] = {"numu", "numubar", "nue", "nuebar"};

    TFile fcv(Form("%s",cvfile.c_str()));
    TFile frw(Form("%s",rwfile.c_str()));
    for (int iptyp=0;iptyp<4;iptyp++) {
      for (int intyp=0;intyp<4;intyp++) {
	for (int ibin=0;ibin<200;ibin++) {
	  fCV[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (fcv.Get(Form("h_%s_%s",ptype[iptyp].c_str(),ntype[intyp].c_str()))))->GetBinContent(ibin+1);
	  fRW[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (frw.Get(Form("h_%s_%s",ptype[iptyp].c_str(),ntype[intyp].c_str()))))->GetBinContent(ibin+1);
	}
      }
    }
    fcv.Close();
    frw.Close();

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));
    fWeightArray.resize(fNmultisims);

    if (fMode.find("multisim") != std::string::npos )
      for (int i=0;i<fNmultisims;i++) fWeightArray[i]=fGaussRandom->shoot(&rng->getEngine(GetName()),0,1.);
    else
      for (int i=0;i<fNmultisims;i++) fWeightArray[i]=1.;
  }

  std::vector<std::vector<double> > FluxHistWeightCalc::GetWeight(art::Event & e)
  {
    //calculate weight(s) here 
    // art::ServiceHandle<art::RandomNumberGenerator> rng;
    //    CLHEP::HepRandomEngine &engine = rng->getEngine(GetName());///***avisar Zarko q lo he quitado
    std::vector<std::vector<double> > weight;

    // * MC flux information
    art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
    std::vector<art::Ptr<simb::MCFlux> > fluxlist;
    if (e.getByLabel(fGenieModuleLabel,mcfluxListHandle))
      art::fill_ptr_vector(fluxlist, mcfluxListHandle);
    else{ return weight;}

    // * MC truth information
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (e.getByLabel(fGenieModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);
    else{return weight;}


    weight.resize(mclist.size());
    for (unsigned int inu=0;inu<mclist.size();inu++) {
      weight[inu].resize(fNmultisims);
      
      int ptype=-9999;
      int ntype=-9999;
      int bin=-9999;
      
      if ( fluxlist[inu]->fptype==211 || fluxlist[inu]->fptype==-211 ) ptype = 0;
      else if ( fluxlist[inu]->fptype==321 || fluxlist[inu]->fptype==-321 ) ptype = 1;
      else if ( fluxlist[inu]->fptype==130 ) ptype = 2;
      else if ( fluxlist[inu]->fptype==13 || fluxlist[inu]->fptype==-13 ) ptype = 3;
      else {
	throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<fluxlist[0]->fptype<< std::endl;
      }
      
      if ( fluxlist[inu]->fntype==14 ) ntype=0;
      else if ( fluxlist[inu]->fntype==-14 ) ntype=1;
      else if ( fluxlist[inu]->fntype==12 ) ntype=2;
      else if ( fluxlist[inu]->fntype==-12 ) ntype=3;
      else {
	throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<fluxlist[0]->fptype<< std::endl;
      }
      
      double enu=mclist[inu]->GetNeutrino().Nu().E();
      bin=enu/0.05;     
      for (int i=0;i<fNmultisims;i++) {
	double test = 1-(1-fRW[ptype][ntype][bin]/fCV[ptype][ntype][bin])*fWeightArray[i];
	
	if(test != test){ test = 1; } // Guards against inifinite weights

	weight[inu][i] = test;
      }
    }
    return weight;
  }

}
}

REGISTER_WEIGHTCALC(sbncode::evwgh::FluxHistWeightCalc)

