#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "StandardRecord/StandardRecord.h"
#include "TCanvas.h"
#include "TH1.h"
#include "CAFAna/Vars/FitVarsSterile.h"
#include "CAFAna/Analysis/FitAxis.h"

// New includes
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/GaussianConstraint.h"
#include "CAFAna/Analysis/ExpInfo.h"

// Random numbers to fake an efficiency and resolution
#include "TRandom3.h"

#include "TMarker.h"

using namespace ana;

// State file
const char* stateFname = "cafe_state.root";

const double sbndPOT = POTnominal;
const double icarusPOT = POTnominal;

void nus()
{

  if (TFile(stateFname).IsZombie()){
    std:: cout << "Run make_state.C first!" << std::endl;
    return;
  }
  else{

    std::cout << "Loading state from " << stateFname << std::endl; 
    TFile fin(stateFname);
    PredictionNoExtrap& pred_nd_numu = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("pred_nd_numu")).release();
    PredictionNoExtrap& pred_fd_numu = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("pred_fd_numu")).release();

    // Calculator
    osc::OscCalculatorSterile* calc = DefaultSterileCalc(4);
    calc->SetAngle(2, 4, 0);
    calc->SetDm(4, 0); // Make these explicit
    calc->SetL(BaselineSBND); 
    calc->SetAngle(2, 3, M_PI/4);
    calc->SetDm(3, 2.45e-3); // make these explicit

    // TMarker* trueValues = new TMarker(pow(TMath::Sin(2*calc->GetAngle(2,4)),2), calc->GetDm(4), kFullCircle);
    // trueValues->SetMarkerColor(kRed);
    
    // To make a fit we need to have a "data" spectrum to compare to our MC
    // Prediction object
    const Spectrum data = pred_nd_numu.Predict(calc).FakeData(sbndPOT);
    SingleSampleExperiment expt(&pred_nd_numu, data);

    TFile* fOutput = new TFile("4flavour/Surfaces_nus.root","RECREATE");

    //Define fit axes
    const FitAxis kAxSinSq2Theta24(&kFitSinSq2Theta24Sterile, 100, 0.001, 1);
    const FitAxis kAxDmSq41(&kFitDmSq41Sterile, 100, 0.001, 100);

    // A Surface evaluates the experiment's chisq across a grid
    Surface surf(&expt, calc,
		 kAxSinSq2Theta24,
		 kAxDmSq41);

    surf.SaveTo(fOutput->mkdir("surf"));

    TCanvas* c1 = new TCanvas("c1");
    c1->SetLogx();
    c1->SetLogy();
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.15);
    //surf.Draw();
    // surf.DrawBestFit(kBlue);
    // trueValues->Draw();

    // In a full Feldman-Cousins analysis you need to provide a critical value
    // surface to be able to draw a contour. But we provide these helper
    // functions to use the gaussian up-values.
    TH2* crit1sig = Gaussian68Percent2D(surf);
    TH2* crit2sig = Gaussian2Sigma2D(surf);

    surf.DrawContour(crit1sig, 7, kBlue);
    surf.DrawContour(crit2sig, kSolid, kBlue);

    c1->SaveAs("4flavour/nus_plot1.pdf");

    // Let's now try adding a second experiment, for instance Icarus
    calc->SetL(BaselineIcarus); // Icarus
    const Spectrum data2 = pred_fd_numu.Predict(calc).FakeData(icarusPOT);
    SingleSampleExperiment expt2(&pred_fd_numu, data2);

    MultiExperimentSBN multiExpt({&expt, &expt2}, {kSBND,kICARUS});

    Surface surf2(&expt2, calc,
		 kAxSinSq2Theta24,
		 kAxDmSq41);
		  
    surf2.SaveTo(fOutput->mkdir("surf2"));

    c1->Clear(); // just in case
    // surf2.DrawBestFit(kGreen+2);
    // trueValues->Draw();
    TH2* crit1sig2 = Gaussian68Percent2D(surf2);
    TH2* crit2sig2 = Gaussian2Sigma2D(surf2);

    surf2.DrawContour(crit1sig2, 7, kGreen+2);
    surf2.DrawContour(crit2sig2, kSolid, kGreen+2);

    c1->SaveAs("4flavour/nus_plot2.pdf");

    Surface surfMulti(&multiExpt, calc,
		      kAxSinSq2Theta24,
		      kAxDmSq41);


    surfMulti.SaveTo(fOutput->mkdir("surfMulti"));

    c1->Clear(); // just in case
    // surfMulti.DrawBestFit(kBlue);
    // trueValues->Draw();
    TH2* crit1sigMulti = Gaussian68Percent2D(surfMulti);
    TH2* crit2sigMulti = Gaussian2Sigma2D(surfMulti);

    surfMulti.DrawContour(crit1sigMulti, 7, kBlue);
    surfMulti.DrawContour(crit2sigMulti, kSolid, kBlue);
    
    c1->SaveAs("4flavour/nus_plot3.pdf");

    fOutput->Close();
  }
}
