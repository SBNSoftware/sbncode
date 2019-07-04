#include "common_fit_definitions.C"

std::string stateFname  = "common_state_ndfd_xsecflux.root";
std::string outputFname = "validation_hists.root";

//Set systematics style by hand for now
bool normsyst = false;
bool fluxsyst = true;
bool xsecsyst = true;
bool use_nd   = true;

void validation_joint(bool reload = false){
  
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

  // Open 
  TFile* fout = new TFile(outputFname.c_str(), "RECREATE");

  osc::NoOscillations noOsc;
  // Use normal hierarchy for now
  osc::IOscCalculatorAdjustable* inputOsc = NuFitOscCalc(+1);

  // Loop over the systematics and make dial validations for each one
  for (auto & syst : systlist){

    std::cout << "Making variation histograms for " << syst->ShortName() << std::endl;
    
    std::vector<TH1*> FD_FHCNumu_uohists = GetMCTotalForSystShifts(&predInterpFDNumuFHC, &noOsc, syst, "FD_FHC_Numu_unosc", pot_fd);
    std::vector<TH1*> FD_FHCNue_uohists  = GetMCTotalForSystShifts(&predInterpFDNueFHC,  &noOsc, syst, "FD_FHC_Nue_unosc",  pot_fd);
    std::vector<TH1*> FD_RHCNumu_uohists = GetMCTotalForSystShifts(&predInterpFDNumuRHC, &noOsc, syst, "FD_RHC_Numu_unosc", pot_fd);
    std::vector<TH1*> FD_RHCNue_uohists  = GetMCTotalForSystShifts(&predInterpFDNueRHC,  &noOsc, syst, "FD_RHC_Nue_unosc",  pot_fd);
    for (auto & hist : FD_FHCNumu_uohists) hist->Write();
    for (auto & hist : FD_FHCNue_uohists)  hist->Write();
    for (auto & hist : FD_RHCNumu_uohists) hist->Write();
    for (auto & hist : FD_RHCNue_uohists)  hist->Write();

    std::vector<TH1*> ND_FHC_hists     = GetMCTotalForSystShifts(&predInterpNDNumuFHC, &noOsc, syst, "ND_FHC", pot_nd);
    std::vector<TH1*> ND_RHC_hists     = GetMCTotalForSystShifts(&predInterpNDNumuRHC, &noOsc, syst, "ND_RHC", pot_nd);    
    std::vector<TH1*> ND_FHC_1Dhists   = GetMCTotalForSystShifts(&predInterpNDNumuFHC, &noOsc, syst, "ND_FHC_1D", pot_nd, true);
    std::vector<TH1*> ND_RHC_1Dhists   = GetMCTotalForSystShifts(&predInterpNDNumuRHC, &noOsc, syst, "ND_RHC_1D", pot_nd, true);  
    for (auto & hist : ND_FHC_hists)   hist->Write();
    for (auto & hist : ND_RHC_hists)   hist->Write();
    for (auto & hist : ND_FHC_1Dhists) hist->Write();
    for (auto & hist : ND_RHC_1Dhists) hist->Write();

    std::vector<TH1*> FD_FHCNumu_hists = GetMCTotalForSystShifts(&predInterpFDNumuFHC, inputOsc, syst, "FD_FHC_Numu", pot_fd);
    std::vector<TH1*> FD_FHCNue_hists  = GetMCTotalForSystShifts(&predInterpFDNueFHC,  inputOsc, syst, "FD_FHC_Nue",  pot_fd);
    std::vector<TH1*> FD_RHCNumu_hists = GetMCTotalForSystShifts(&predInterpFDNumuRHC, inputOsc, syst, "FD_RHC_Numu", pot_fd);
    std::vector<TH1*> FD_RHCNue_hists  = GetMCTotalForSystShifts(&predInterpFDNueRHC,  inputOsc, syst, "FD_RHC_Nue",  pot_fd);
    for (auto & hist : FD_FHCNumu_hists) hist->Write();
    for (auto & hist : FD_FHCNue_hists)  hist->Write();
    for (auto & hist : FD_RHCNumu_hists) hist->Write();
    for (auto & hist : FD_RHCNue_hists)  hist->Write();
  }
  
  fout->Close();
}

