#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Vars/FitVarsSterile.h"
#include "CAFAna/Analysis/FitAxis.h"

// New includes
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Analysis/MedianSurface.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/GaussianConstraint.h"
#include "CAFAna/Analysis/ExpInfo.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

#include "toysysts.h"

using namespace ana;

// State file
const char* basicFname = "cafe_state_smear.root";

const double sbndPOT = POTnominal;
const double icarusPOT = POTnominal;

void nus(const char* stateFname = basicFname, int nmock = 0, bool useSysts = true)
{

  if (TFile(stateFname).IsZombie()){
    std:: cout << "Run make_state.C first!" << std::endl;
    return;
  }
  else{

    std::cout << "Loading state from " << stateFname << std::endl; 
    TFile fin(stateFname);
    PredictionInterp& pred_nd_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd_numu")).release();
    PredictionInterp& pred_fd_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu")).release();

    // Calculator
    osc::OscCalculatorSterile* calc = DefaultSterileCalc(4);
    calc->SetL(BaselineSBND); 

    // To make a fit we need to have a "data" spectrum to compare to our MC
    // Prediction object
    const Spectrum data = pred_nd_numu.Predict(calc).FakeData(sbndPOT);
    SingleSampleExperiment expt(&pred_nd_numu, data);

    TFile* fOutput = new TFile("Surfaces_nus.root","RECREATE");

    //Define fit axes
    const FitAxis kAxSinSq2Theta24(&kFitSinSq2Theta24Sterile, 40, 1e-3, 1, true);
    const FitAxis kAxDmSq41(&kFitDmSq41Sterile, 40, 0.1, 100, true);

    if(!useSysts) allSysts.clear();

    // A Surface evaluates the experiment's chisq across a grid
    Surface surf(&expt, calc,
		 kAxSinSq2Theta24,
		 kAxDmSq41,
                 {},
		 allSysts);

    surf.SaveTo(fOutput->mkdir("surf"));

    TCanvas* c1 = new TCanvas("c1");
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.15);

    // In a full Feldman-Cousins analysis you need to provide a critical value
    // surface to be able to draw a contour. But we provide these helper
    // functions to use the gaussian up-values.
    TH2* crit2sig = Gaussian3Sigma1D1Sided(surf);

    surf.DrawContour(crit2sig, kSolid, kBlue);

    c1->SaveAs(useSysts ? "nus_SBND.pdf" : "nus_SBND_statsOnly.pdf");

    // Let's now try adding a second experiment, for instance Icarus
    calc->SetL(BaselineIcarus); // Icarus
    const Spectrum data2 = pred_fd_numu.Predict(calc).FakeData(icarusPOT);
    SingleSampleExperiment expt2(&pred_fd_numu, data2);

    MultiExperimentSBN multiExpt({&expt, &expt2}, {kSBND,kICARUS});

    Surface surf2(&expt2, calc,
                  kAxSinSq2Theta24,
                  kAxDmSq41,
                  {},
		  allSysts);
		  
    surf2.SaveTo(fOutput->mkdir("surf2"));

    c1->Clear(); // just in case

    TH2* crit2sig2 = Gaussian3Sigma1D1Sided(surf2);

    surf2.DrawContour(crit2sig2, kSolid, kGreen+2);

    c1->SaveAs(useSysts ? "nus_ICARUS.pdf" : "nus_ICARUS_statsOnly.pdf");

    Surface surfMulti(&multiExpt, calc,
		      kAxSinSq2Theta24,
		      kAxDmSq41,
		      {},
		      allSysts);

    std::vector<Surface> mockSurfs;
    for(int i = 0; i < nmock; ++i){
      std::cout << "Mock " << i+1 << " / " << nmock << std::endl;
      const FitAxis kCoarseAxSinSq2Theta24(&kFitSinSq2Theta24Sterile, 20, 1e-3, 1, true);
      const FitAxis kCoarseAxDmSq41(&kFitDmSq41Sterile, 20, 0.1, 100, true);

      osc::OscCalculatorSterile* c = DefaultSterileCalc(4);
      c->SetL(BaselineSBND); 
      SingleSampleExperiment e1(&pred_nd_numu, pred_nd_numu.Predict(c).MockData(sbndPOT));

      c->SetL(BaselineIcarus); // Icarus
      SingleSampleExperiment e2(&pred_fd_numu, pred_fd_numu.Predict(c).MockData(icarusPOT));

      MultiExperimentSBN me({&e1, &e2}, {kSBND, kICARUS});

      Surface ms(&me, c,
                 kCoarseAxSinSq2Theta24,
                 kCoarseAxDmSq41,
                 {},
		 allSysts);

      mockSurfs.push_back(ms);
    }

    surfMulti.SaveTo(fOutput->mkdir("surfMulti"));

    c1->Clear(); // just in case

    TH2* crit2sigMulti = Gaussian3Sigma1D1Sided(surfMulti);

    surfMulti.DrawContour(crit2sigMulti, kSolid, kRed);
    
    c1->SaveAs(useSysts ? "nus_BOTH.pdf" : "nus_BOTH_statsOnly.pdf");

    c1->Clear();

    if(!mockSurfs.empty()){
      MedianSurface ms(mockSurfs);
      TH2* crit2 = Gaussian3Sigma1D1Sided(ms);
      ms.DrawEnsemble(crit2);
      ms.DrawContour(crit2, kSolid, kRed);

      ms.SaveTo(fOutput->mkdir("median_surf"));
    }

    surfMulti.DrawContour(crit2sigMulti, nmock > 0 ? 7 : kSolid, kRed);
    surf.DrawContour(crit2sigMulti, kSolid, kBlue);
    surf2.DrawContour(crit2sigMulti, kSolid, kGreen+2);

    TH1F* hDummy1 = new TH1F("","",1,0,1);
    TH1F* hDummy2 = new TH1F("","",1,0,1);
    TH1F* hDummy3 = new TH1F("","",1,0,1);

    hDummy1->SetLineColor(kRed);
    hDummy2->SetLineColor(kBlue);
    hDummy3->SetLineColor(kGreen+2);

    TLegend* leg = new TLegend(0.6,0.6,0.84,0.84);
    leg->SetHeader("3#sigma C.L.","C");
    leg->AddEntry(hDummy2, "SBND only", "l");
    leg->AddEntry(hDummy3, "Icarus only", "l");
    leg->AddEntry(hDummy1, "Combined fit", "l");
    leg->SetLineColor(10);
    leg->SetFillColor(10);
    leg->Draw();

    c1->SaveAs(useSysts ? "nus_ALL.root" : "nus_ALL_statsOnly.root");
    c1->SaveAs(useSysts ? "nus_ALL.pdf" : "nus_ALL_statsOnly.pdf");

    fOutput->Close();
  }
}
