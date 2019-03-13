// ETW May 2018
// Standard script for DUNE spectra
// Input files use TensorFlow CVN training from May 2018 

#include "common_fit_definitions.C"

std::string stateFname  = "common_state_ndfd_nosyst.root";
std::string outputFname = "spec_hist.root";

//Set systematics style by hand for now
bool normsyst = false;
bool fluxsyst = false;
bool use_nd   = true;

void spec_joint(bool reload = false){
  
  gROOT->SetBatch(1);
  
  // Get the systematics to use
  std::vector<const ISyst*> systlist = {};

  // Get the prediction interpolators
  std::vector<unique_ptr<PredictionInterp> > return_list = GetPredictionInterps(stateFname, systlist);
  PredictionInterp& predInterpFDNumuFHC = *return_list[0].release();
  PredictionInterp& predInterpFDNueFHC  = *return_list[1].release();
  PredictionInterp& predInterpFDNumuRHC = *return_list[2].release();
  PredictionInterp& predInterpFDNueRHC  = *return_list[3].release();
  PredictionInterp& predInterpNDNumuFHC = *return_list[4].release();
  PredictionInterp& predInterpNDNumuRHC = *return_list[5].release();

  // Open 
  TFile* fout = new TFile(outputFname.c_str(), "RECREATE");

  osc::NoOscillations noOsc;

  // Unoscillated FD histograms
  std::vector<TH1*> FD_FHCNumu_uohists = GetMCComponents(&FD_predFHCNumu, &noOsc, "FD_FHC_Numu_unosc", pot_fd);
  std::vector<TH1*> FD_FHCNue_uohists  = GetMCComponents(&FD_predFHCNue, &noOsc, "FD_FHC_Nue_unosc", pot_fd);
  std::vector<TH1*> FD_RHCNumu_uohists = GetMCComponents(&FD_predRHCNumu, &noOsc, "FD_RHC_Numu_unosc", pot_fd);
  std::vector<TH1*> FD_RHCNue_uohists  = GetMCComponents(&FD_predRHCNue, &noOsc, "FD_RHC_Nue_unosc", pot_fd);
  for (auto & hist : FD_FHCNumu_uohists) hist->Write();
  for (auto & hist : FD_FHCNue_uohists)  hist->Write();
  for (auto & hist : FD_RHCNumu_uohists) hist->Write();
  for (auto & hist : FD_RHCNue_uohists)  hist->Write();
  
  // Sort out the ND histograms
  std::vector<TH1*> ND_FHC_hists     = GetMCComponents(&ND_predFHC, &noOsc, "ND_FHC", pot_nd);
  std::vector<TH1*> ND_RHC_hists     = GetMCComponents(&ND_predRHC, &noOsc, "ND_RHC", pot_nd);    

  std::vector<TH1*> ND_FHC_1Dhists   = GetMCComponents(&ND_predFHC, &noOsc, "ND_FHC_1D", pot_nd, true);
  std::vector<TH1*> ND_RHC_1Dhists   = GetMCComponents(&ND_predRHC, &noOsc, "ND_RHC_1D", pot_nd, true);  
  
  for (auto & hist : ND_FHC_hists)   hist->Write();
  for (auto & hist : ND_RHC_hists)   hist->Write();
  for (auto & hist : ND_FHC_1Dhists) hist->Write();
  for (auto & hist : ND_RHC_1Dhists) hist->Write();  
  
  std::string dcpnames[4] = {"0pi","piover2","pi","3piover2"};

  for(int hie = -1; hie <= +1; hie += 2){

    osc::IOscCalculatorAdjustable* inputOsc = NuFitOscCalc(hie);

    const std::string hieStr = (hie > 0) ? "nh" : "ih";

    for(int deltaIdx = 0; deltaIdx < 4; ++deltaIdx){
      inputOsc->SetdCP(deltaIdx/2.*TMath::Pi());
      const std::string dcpStr = dcpnames[deltaIdx];

      // FD components for this set of parameters
      std::vector<TH1*> FD_FHCNumu_hists = GetMCComponents(&FD_predFHCNumu, inputOsc, "FD_FHC_Numu_"+hieStr+"_"+dcpStr, pot_fd);
      std::vector<TH1*> FD_FHCNue_hists  = GetMCComponents(&FD_predFHCNue, inputOsc, "FD_FHC_Nue_"+hieStr+"_"+dcpStr, pot_fd);
      std::vector<TH1*> FD_RHCNumu_hists = GetMCComponents(&FD_predRHCNumu, inputOsc, "FD_RHC_Numu_"+hieStr+"_"+dcpStr, pot_fd);
      std::vector<TH1*> FD_RHCNue_hists  = GetMCComponents(&FD_predRHCNue, inputOsc, "FD_RHC_Nue_"+hieStr+"_"+dcpStr, pot_fd);
            
      for (auto & hist : FD_FHCNumu_hists) hist->Write();
      for (auto & hist : FD_FHCNue_hists)  hist->Write();
      for (auto & hist : FD_RHCNumu_hists) hist->Write();
      for (auto & hist : FD_RHCNue_hists)  hist->Write();
      
    }
  }
  fout->Close();
}

