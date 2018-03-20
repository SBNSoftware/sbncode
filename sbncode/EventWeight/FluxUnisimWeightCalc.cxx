#include "larsim/EventWeight/Base/WeightCalc.h"

#include <iostream>

#include "TH1D.h"
#include "TH1F.h"
#include "TFile.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nutools/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandGaussQ.h"

namespace sbncode {
namespace evwgh {
  class FluxUnisimWeightCalc : public ::evwgh::WeightCalc
  {

  public:
    FluxUnisimWeightCalc();
    double MiniBooNEWeightCalc(double enu, int ptype, int ntype, int uni, bool noNeg);
    double MicroBooNEWeightCalc(double enu, int ptype, int ntype, int uni);
    void Configure(fhicl::ParameterSet const& p);
    std::vector<std::vector<double> > GetWeight(art::Event & e);

  private:
    CLHEP::RandGaussQ *fGaussRandom;
    std::vector<double> fWeightArray;
    int fNuni;   // Number of universes you want to generate
    double fScalePos; // Scale factor to enhance or degrade a given positive systematic uncertanity
    double fScaleNeg; // Scale factor to enhance or degrade a given negative systematic uncertanity
    std::string fMode;   // if you want multisim or +/- 1sigma
    std::string fGenieModuleLabel;
    std::string fWeightCalc; // if you want MiniBooNE or MicroBooNE calculator

    //Weight Arrays
    // These will contain the systematic variations
    // with these in place we can start to calculate weights
    // these will depend on the parent particle, neutrino flavor, and neutrino energy
    //         ptype: pi, k, k0, mu
    //         |  ntype: numu, numubar, nue, nuebar 
    //         |  |   binnig: 50MeV
    //         |  |   |
    double fCV[4][4][200];
    double fRWpos[4][4][200];
    double fRWneg[4][4][200];

    // This is for when there is only one systematic variation (i.e. skin depth)
    bool PosOnly = false;
       
     DECLARE_WEIGHTCALC(FluxUnisimWeightCalc)
  };

  FluxUnisimWeightCalc::FluxUnisimWeightCalc()
  {

  }

  void FluxUnisimWeightCalc::Configure(fhicl::ParameterSet const& p)
  {
    //Collect the fcl parameters from the fhicl file
    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    fNuni  = pset.get<int>("number_of_multisims");
    fScalePos = pset.get<double>("scale_factor_pos");
    fScaleNeg = pset.get<double>("scale_factor_neg");
    fMode  = pset.get<std::string>("mode");	
    fWeightCalc = pset.get<std::string>("weight_calculator");

    // This sets whether we use the MiniBooNE or a modern calculation for the weights  
    if(fWeightCalc != "MicroBooNE" && fWeightCalc != "MiniBooNE"){
      throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown weight calculator "<< fWeightCalc << std::endl;     
    }

    fGenieModuleLabel= p.get< std::string > ("genie_module_label");
    
    /// Grab the histogram related to the CV
    std::string dataInput1       =   pset.get< std::string >("CentralValue_hist_file");

    /// Grab the histogram related to the variation 
    std::string dataInput2pos  =   pset.get< std::string >("PositiveSystematicVariation_hist_file");
    std::string dataInput2neg  =   pset.get< std::string >("NegativeSystematicVariation_hist_file");

    //
    //  If there is only one file supplied use 
    //     it for both positive and negative
    //
    if(dataInput2pos == dataInput2neg){
      PosOnly = true;
    }


    cet::search_path sp("FW_SEARCH_PATH");
    std::string rwfilepos      = sp.find_file(dataInput2pos);
    std::string rwfileneg      = sp.find_file(dataInput2neg);
    std::string cvfile         = sp.find_file(dataInput1);
    //Set up the naming convention of the 
    //   histograms that contain all the information we have 
    int ptype[4] = {1,2,3,4}; //mu, pi, k0, k
    int ntype[4] = {1,2,3,4}; //nue, anue, numu, anumu

        // Define the files that store the histograms
    TFile fcv(Form("%s",cvfile.c_str()));
    TFile frwpos(Form("%s", rwfilepos.c_str()));
    TFile frwneg(Form("%s", rwfileneg.c_str()));

    // > Iterate through all the files and extract 
    //       the relevant histograms, bin by bin
    // > The histograms are labled as
    //       'h5'+ptype[]+ntype[] 

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


    //Setup the random number generator
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));
    
    //
    //   This part is important!!! You want to be sure to use the same random number throughout all the 
    //   events, otherwise you will be smearing over all the correlations. 
    //
    //   While the random number seed it set in the fcl file, you only want to use the same throw PER UNIVERSE
    //
    fWeightArray.resize(fNuni);
    
    // This sets up if we want to reweight events or just assess 1sigma shifts 
    if (fMode.find("multisim") != std::string::npos )
      for (int i=0;i<fNuni;i++) fWeightArray[i]=fGaussRandom->shoot(&rng->getEngine(GetName()),0,1.);
    else
      for (int i=0;i<fNuni;i++) fWeightArray[i]=1.;

  }

  std::vector<std::vector<double> > FluxUnisimWeightCalc::GetWeight(art::Event & e)
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

    // No neutrinos in this event
    if(mclist.size() == 0) return weight;

    //Iterate through each neutrino in the event
    for(unsigned int inu = 0; inu < mclist.size(); inu++){

      //Resize vector to the number of universes you want to generate
      weight[inu].resize(fNuni);
      
      //containers for the parent and neutrino type information
      int ptype = std::numeric_limits<int>::max(); 
      int ntype = std::numeric_limits<int>::max();
 
     // Discover the neutrino parent type
      //     This contains the neutrino's parentage information
      if (      fluxlist[inu].fptype==13  || fluxlist[inu].fptype==-13  ) ptype = 0;
      else if ( fluxlist[inu].fptype==211 || fluxlist[inu].fptype==-211 ) ptype = 1;
      else if ( fluxlist[inu].fptype==130                               ) ptype = 2;
      else if ( fluxlist[inu].fptype==321 || fluxlist[inu].fptype==-321 ) ptype = 3;                                    
      else {
	throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<fluxlist[0].fptype<< std::endl;
      }

      // Discover the neutrino type
      //     This contains the neutrino's flavor information
      if (      fluxlist[inu].fntype==12  ) ntype=0;
      else if ( fluxlist[inu].fntype==-12 ) ntype=1;
      else if ( fluxlist[inu].fntype==14  ) ntype=2;
      else if ( fluxlist[inu].fntype==-14 ) ntype=3;
      else {
	throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<fluxlist[0].fptype<< std::endl;
      }

      // Collect neutrino energy
      double enu=mclist[inu].GetNeutrino().Nu().E();      

      //Let's make a weights based on the calculator you have requested 
      if(fMode.find("multisim") != std::string::npos){
	for (int i=0;i<fNuni;i++) {

	  if(fWeightCalc.find("MicroBooNE") != std::string::npos){
	    weight[inu][i]=MicroBooNEWeightCalc(enu, ptype, ntype, i);
	  }
	  if(fWeightCalc.find("MiniBooNE") != std::string::npos){
	    weight[inu][i]=MiniBooNEWeightCalc(enu, ptype, ntype, i, PosOnly);
	  }

	}//Iterate through the number of universes      
      }
    }
     
    return weight;
  }

  double FluxUnisimWeightCalc::MiniBooNEWeightCalc(double enu, int ptype, int ntype, int uni, bool noNeg)
  {
    // 
    //   This is directly based on the MiniBooNE framework to 
    //   allow us to reproduce their weights before
    //         All this is a reproduction by J. Zennamo in 2017 
    //                      of code written by S. Brice in 2008
    //                                         
    double weight = 1;
    
    int bin = int(enu/0.05); //convert energy (in GeV) into 50 MeV bin


    //  This is based on:
    //    http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/~checkout~/...
    //          miniboone/AnalysisFramework/MultisimMatrix/src/MultisimMatrix_initialise.F?rev=1.18;content-type=text%2Fplain
    //
    //  pseudocode:
    //   Scaled Reweighting = ScaleFactor * Reweighting + ( 1 - ScaleFactor) * Central Value
    //
    double scaled_pos = fScalePos*fRWpos[ptype][ntype][bin] + 
      (1-fScalePos)*fCV[ptype][ntype][bin];
    
    double scaled_neg = fScaleNeg*fRWneg[ptype][ntype][bin] + 
      (1-fScaleNeg)*fCV[ptype][ntype][bin];
    
    // This is based on:
    //    http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/~checkout~/...
    //           miniboone/AnalysisFramework/MultisimMatrix/src/MultisimMatrix_getWeight.F?rev=1.41;content-type=text%2Fplain
    //
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

    if(fCV[ptype][ntype][bin] == 0) weight = 1;
    if(fCV[ptype][ntype][bin] < 1.e-12) weight = 1;

    if(weight < 0) weight = 1; 
    if(weight > 30) weight = 30;
    if(weight != weight) weight = 30;

    if( (ntype == 0 || ntype == 1) && ptype == 1) weight = 1;

    if( (ntype == 1 || ntype == 3) && ptype == 3) weight = 1;

     
    return weight;
  }
  
  double FluxUnisimWeightCalc::MicroBooNEWeightCalc(double enu, int ptype, int ntype, int uni)
  {
   
    double weight = 1;
    /*    int bin = int(enu/0.05); //convert energy (in GeV) into 50 MeV bin
    double Scaled_syst_shift = 0;

    if(fWeightArray[uni] >= 0 ){ 
      Scaled_syst_shift = (fScalePos)*(fRWpos[ptype][ntype][bin]) + (1-fScalePos)*(fCV[ptype][ntype][bin]);
      weight = 1+(fWeightArray[uni]*((Scaled_syst_shift/fCV[ptype][ntype][bin])-1));
    }
    else{  
      Scaled_syst_shift = (fScaleNeg)*(fRWneg[ptype][ntype][bin]) + (1-fScaleNeg)*(fCV[ptype][ntype][bin]);      
      weight = 1-(fWeightArray[uni]*((Scaled_syst_shift/fCV[ptype][ntype][bin])-1));
    }    
    */
    return weight;
  }

}
}

REGISTER_WEIGHTCALC(sbncode::evwgh::FluxUnisimWeightCalc)

