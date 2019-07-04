// Standard script for DUNE MH sensitivity
// Now uses state file created by spec.C

#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/CalcsNuFit.h"

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Progress.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Systs/Systs.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Analysis/Exposures.h"

using namespace ana;

#include "Utilities/rootlogon.C"

#include "OscLib/func/IOscCalculator.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"

const Var kRecoE_nue = SIMPLEVAR(dune.Ev_reco_nue);
const Var kRecoE_numu = SIMPLEVAR(dune.Ev_reco_numu);

// POT/yr * 3.5yrs * mass correction
const double potFD = 3.5 * POT120 * 40/1.13;

const char* stateFname = "spec_state.root";
const char* outputFname = "mh_sens.root";

//Set systematics style by hand for now
bool nosyst = false;
bool normsyst = true;
bool fullsyst = false;

std::vector<const ISyst*> systlist;
std::vector<const ISyst*> normlist_sig = {&kNueFHCSyst, &kNumuFHCSyst, &kNueRHCSyst, &kNumuRHCSyst};
std::vector<const ISyst*> normlist_bg = {&kNueBeamFHCSyst, &kNueBeamRHCSyst, &kNCAppSyst, &kNCDisSyst, &kNutauSyst};

void mh()
{
  rootlogon(); // style
  
  if(TFile(stateFname).IsZombie()){
    std::cout << "No state available, please run spec.C true" << std::endl;
    }
  else{    
    if (normsyst) {
      systlist.insert(systlist.end(), normlist_sig.begin(), normlist_sig.end()); 
      systlist.insert(systlist.end(), normlist_bg.begin(), normlist_bg.end()); 
    }

    // List all of the systematics we'll be using
    std::cout << "Systematics in this fit: " << std::endl;
    for(const ISyst* s: systlist) std::cout << s->ShortName() << "\t\t" << s->LatexName() << std::endl;
    if (systlist.size()==0) {std::cout << "None" << std::endl;}

    std::cout << "Loading state from " << stateFname << std::endl; 
  }

  TFile fin(stateFname);
  PredictionInterp& predInt_FDNumuFHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu_fhc")).release();
  PredictionInterp& predInt_FDNueFHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_nue_fhc")).release();
  PredictionInterp& predInt_FDNumuRHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu_rhc")).release();
  PredictionInterp& predInt_FDNueRHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_nue_rhc")).release();

  fin.Close();
  std::cout << "Done loading state" << std::endl;

  TFile* fout = new TFile(outputFname, "RECREATE");

  for(int hie = -1; hie <= +1; hie += 2){

    bool oscvar = true;
    double dcpstep = 2*TMath::Pi()/36;
    TGraph* gMH_NH = new TGraph();
    TGraph* gMH_IH = new TGraph();
    double thisdcp;

    for(double idcp = 0; idcp < 37; ++idcp) {
	
      thisdcp = -TMath::Pi() + idcp*dcpstep;
	
      osc::IOscCalculatorAdjustable* trueOsc = NuFitOscCalc(hie);
      trueOsc->SetdCP(thisdcp);

      const Spectrum data_nue_fhc_syst = predInt_FDNueFHC.Predict(trueOsc).FakeData(potFD);
      SingleSampleExperiment app_expt_fhc_syst(&predInt_FDNueFHC, data_nue_fhc_syst);
      app_expt_fhc_syst.SetMaskHist(0.5,8.0);

      const Spectrum data_nue_rhc_syst = predInt_FDNueRHC.Predict(trueOsc).FakeData(potFD);
      SingleSampleExperiment app_expt_rhc_syst(&predInt_FDNueRHC, data_nue_rhc_syst);
      app_expt_rhc_syst.SetMaskHist(0.5,8.0);

      const Spectrum data_numu_fhc_syst = predInt_FDNumuFHC.Predict(trueOsc).FakeData(potFD);
      SingleSampleExperiment dis_expt_fhc_syst(&predInt_FDNumuFHC, data_numu_fhc_syst);
      dis_expt_fhc_syst.SetMaskHist(0.5,8.0);

      const Spectrum data_numu_rhc_syst = predInt_FDNumuRHC.Predict(trueOsc).FakeData(potFD);
      SingleSampleExperiment dis_expt_rhc_syst(&predInt_FDNumuRHC, data_numu_rhc_syst);
      dis_expt_rhc_syst.SetMaskHist(0.5,8.0);

      std::vector<const IFitVar*> oscVars =
	{&kFitDmSq32Scaled, &kFitSinSqTheta23,
	 &kFitTheta13, &kFitDeltaInPiUnits, &kFitRho};

      double chisqmin = 99999;
      double thischisq;

      osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(hie);	
      testOsc->SetdCP(thisdcp);
      testOsc->SetDmsq32(-1*testOsc->GetDmsq32());
      for(int ioct = -1; ioct <= 1; ioct += 2) {
	if (ioct < 0) {
	  testOsc->SetTh23(TMath::PiOver2() - testOsc->GetTh23());
	}

	osc::IOscCalculatorAdjustable* cvcalc = testOsc->Copy();	  
	Penalizer_GlbLike penalty(cvcalc,hie);

	MultiExperiment full_expt_syst({&app_expt_fhc_syst, &app_expt_rhc_syst, &dis_expt_fhc_syst, &dis_expt_rhc_syst, &penalty});

	Fitter fit_syst(&full_expt_syst, oscVars, systlist);

	if (nosyst) {
	  continue; //etw figure out how to do this
	}
	else {
	  thischisq = fit_syst.Fit(testOsc);
	}

	chisqmin = TMath::Min(thischisq,chisqmin);
      }

      chisqmin = TMath::Max(chisqmin,1e-6);
      std::cout << thisdcp << " " << TMath::Sqrt(chisqmin) << std::endl;
      if (hie > 0) {
	gMH_NH->SetPoint(gMH_NH->GetN(),thisdcp/TMath::Pi(),TMath::Sqrt(chisqmin));
      }
      else {
	gMH_IH->SetPoint(gMH_IH->GetN(),thisdcp/TMath::Pi(),TMath::Sqrt(chisqmin));
      }
    }

    if (hie > 0) {
      gMH_NH->Draw("ALP");
      gMH_NH->Write("sens_mh_nh");
    }
    else {
      gMH_IH->Draw("ALP");
      gMH_IH->Write("sens_mh_ih");
    }
  }
  fout->Close();
}
