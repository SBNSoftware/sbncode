// ETW May 2018
// Standard script for DUNE CPV sensitivity
// Input files use TensorFlow CVN training from May 2018 

#include "common_fit_definitions.C"

std::string stateFname  = "common_state_ndfd_xsecflux.root";
std::string outputFname = "cpv_sens_fd_xsecflux.root";

//Set systematics style by hand for now
bool normsyst = false;
bool fluxsyst = true;
bool xsecsyst = true;
bool use_nd   = false;

// Need to accept filename, ND/FD, systs and reload as arguments
void cpv_joint(bool reload = false){
  
  gROOT->SetBatch(1);

  // Get the systematics to use
  std::vector<const ISyst*> systlist = GetListOfSysts(fluxsyst, xsecsyst, normsyst);

  // Get the prediction interpolators
  std::vector<unique_ptr<PredictionInterp> > return_list = GetPredictionInterps(stateFname, systlist);
  PredictionInterp& predInterpFDNumuFHC = *return_list[0].release();
  PredictionInterp& predInterpFDNueFHC  = *return_list[1].release();
  PredictionInterp& predInterpFDNumuRHC = *return_list[2].release();
  PredictionInterp& predInterpFDNueRHC  = *return_list[3].release();
  PredictionInterp& predInterpNDNumuFHC = *return_list[4].release();
  PredictionInterp& predInterpNDNumuRHC = *return_list[5].release();

  TFile* fout = new TFile(outputFname.c_str(), "RECREATE");

  for(int hie = -1; hie <= +1; hie += 2){

    const std::string hieStr = (hie > 0) ? "nh" : "ih";
    bool oscvar = true;
    double dcpstep = 2*TMath::Pi()/36;
    TGraph* gCPV_NH = new TGraph();
    TGraph* gCPV_IH = new TGraph();
    double thisdcp;

    for(double idcp = 0; idcp < 37; ++idcp) {
	
      thisdcp = -TMath::Pi() + idcp*dcpstep;

      osc::IOscCalculatorAdjustable* trueOsc = NuFitOscCalc(hie);
      trueOsc->SetdCP(thisdcp);

      const Spectrum data_nue_fhc_syst = predInterpFDNueFHC.Predict(trueOsc).FakeData(pot_fd);
      SingleSampleExperiment app_expt_fhc_syst(&predInterpFDNueFHC, data_nue_fhc_syst);
      app_expt_fhc_syst.SetMaskHist(0.5, 8);
      
      const Spectrum data_nue_rhc_syst = predInterpFDNueRHC.Predict(trueOsc).FakeData(pot_fd);
      SingleSampleExperiment app_expt_rhc_syst(&predInterpFDNueRHC, data_nue_rhc_syst);
      app_expt_rhc_syst.SetMaskHist(0.5, 8);
      
      const Spectrum data_numu_fhc_syst = predInterpFDNumuFHC.Predict(trueOsc).FakeData(pot_fd);
      SingleSampleExperiment dis_expt_fhc_syst(&predInterpFDNumuFHC, data_numu_fhc_syst);
      dis_expt_fhc_syst.SetMaskHist(0.5, 8);
      
      const Spectrum data_numu_rhc_syst = predInterpFDNumuRHC.Predict(trueOsc).FakeData(pot_fd);
      SingleSampleExperiment dis_expt_rhc_syst(&predInterpFDNumuRHC, data_numu_rhc_syst);
      dis_expt_rhc_syst.SetMaskHist(0.5, 8);
      
      const Spectrum nd_data_numu_fhc_syst = predInterpNDNumuFHC.PredictUnoscillated().FakeData(pot_nd);
      SingleSampleExperiment nd_expt_fhc_syst(&predInterpNDNumuFHC, nd_data_numu_fhc_syst);
      nd_expt_fhc_syst.SetMaskHist(0, -1, 0.1, 1);
      
      const Spectrum nd_data_numu_rhc_syst = predInterpNDNumuRHC.PredictUnoscillated().FakeData(pot_nd);
      SingleSampleExperiment nd_expt_rhc_syst(&predInterpNDNumuRHC, nd_data_numu_rhc_syst);
      nd_expt_rhc_syst.SetMaskHist(0, -1, 0.1, 1);
      
      std::vector<const IFitVar*> oscVars =
	{&kFitDmSq32Scaled, &kFitSinSqTheta23,
	 &kFitTheta13, &kFitRho};

      double chisqmin = 99999;
      double thischisq;

      for(int ihie = -1; ihie <= +1; ihie += 2) {
	for (int idcp = 0; idcp < 2; ++idcp) {
	  for (int ioct = -1; ioct <= 1; ioct +=2) {
	    osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(hie);	
	    double dcptest = idcp*TMath::Pi();
	    testOsc->SetdCP(dcptest);

	    if (ihie < 0) testOsc->SetDmsq32(-1*testOsc->GetDmsq32());
	    if (ioct < 0) testOsc->SetTh23(TMath::Pi()/2 - testOsc->GetTh23());

	    osc::IOscCalculatorAdjustable* cvcalc = testOsc->Copy();	  
	    Penalizer_GlbLike penalty(cvcalc,hie);

	    MultiExperiment expt_fd({&app_expt_fhc_syst, &app_expt_rhc_syst, &dis_expt_fhc_syst, &dis_expt_rhc_syst, &penalty});	    
	    MultiExperiment expt_nd_fd({&app_expt_fhc_syst, &app_expt_rhc_syst, &dis_expt_fhc_syst, &dis_expt_rhc_syst,
	       	  &nd_expt_fhc_syst, &nd_expt_rhc_syst, &penalty});
	    
	    Fitter fit_fd(&expt_fd, oscVars, systlist);
	    Fitter fit_nd_fd(&expt_nd_fd, oscVars, systlist);

	    if(use_nd) thischisq = fit_nd_fd.Fit(testOsc);
	    else thischisq = fit_fd.Fit(testOsc);
	    chisqmin = TMath::Min(thischisq,chisqmin);
	  }
	}
      }
      
      chisqmin = TMath::Max(chisqmin,1e-6);
      if (hie > 0) gCPV_NH->SetPoint(gCPV_NH->GetN(),thisdcp/TMath::Pi(),TMath::Sqrt(chisqmin));
      else gCPV_IH->SetPoint(gCPV_IH->GetN(),thisdcp/TMath::Pi(),TMath::Sqrt(chisqmin));
    }

    if (hie > 0) {
      gCPV_NH->Draw("ALP");
      gCPV_NH->Write("sens_cpv_nh");
    } else {
      gCPV_IH->Draw("ALP");
      gCPV_IH->Write("sens_cpv_ih");
    }
  }
  fout->Close();
}
