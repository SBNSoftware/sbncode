#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Vars/FitVarsSterileApprox.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/RatioExperiment.h"
#include "CAFAna/Experiment/CountingExperiment.h"
#include "CAFAna/Analysis/ExpInfo.h"
#include "CAFAna/Analysis/Surface.h"
using namespace ana;

#include "OscLib/IOscCalculator.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"

#include <vector>

const double sbndPOT = kPOTnominal/100;
const double icarusPOT = kPOTnominal/100;

void WithBand(const std::vector<TH1*>& hs)
{
  TH1* median = (TH1*)hs[0]->Clone();
  median->SetDirectory(0);

  TGraphAsymmErrors* g = new TGraphAsymmErrors;

  for(int i = 0; i < hs[0]->GetNbinsX()+2; ++i){
    std::vector<double> ys;
    for(unsigned int j = 0; j < hs.size(); ++j)
      ys.push_back(hs[j]->GetBinContent(i));
    std::sort(ys.begin(), ys.end());
    const double y   = ys[.50*ys.size()];
    const double ylo = ys[.16*ys.size()];
    const double yhi = ys[.84*ys.size()];

    median->SetBinContent(i, y);

    const double w = hs[0]->GetXaxis()->GetBinWidth(i);

    g->SetPoint(i, hs[0]->GetXaxis()->GetBinCenter(i), y);
    g->SetPointError(i, w/2, w/2, y-ylo, yhi-y);
  }

  median->SetLineColor(kRed);
  median->Draw("hist ][");
  g->SetFillColor(kRed-10);
  g->Draw("e2 same");
  median->Draw("hist ][ same");

  median->GetYaxis()->SetRangeUser(0, 1.3*median->GetMaximum());
}

void plot_multi()
{
  TFile fin("cafe_state_smear_numu.root");

  const int Nuniv = 1000;

  std::vector<PredictionInterp*> preds_nd(Nuniv), preds_fd(Nuniv);

  for(int i = 0; i < Nuniv; ++i){
    preds_nd[i] = LoadFrom<PredictionInterp>(fin.GetDirectory(TString::Format("pred_nd_numu_%d", i).Data())).release();
    preds_fd[i] = LoadFrom<PredictionInterp>(fin.GetDirectory(TString::Format("pred_fd_numu_%d", i).Data())).release();
  }

  PredictionInterp* nom_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd_numu_nom")).release();
  PredictionInterp* nom_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu_nom")).release();

  //  osc::NoOscillations calc;
  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  std::vector<TH1*> hs_nd;
  for(const IPrediction* p: preds_nd) hs_nd.push_back(p->Predict(calc).ToTH1(sbndPOT));
  WithBand(hs_nd);
  gPad->Print("spect_sbnd_multi.pdf");

  new TCanvas;

  std::vector<TH1*> hs_fd;
  for(const IPrediction* p: preds_fd) hs_fd.push_back(p->Predict(calc).ToTH1(icarusPOT));
  WithBand(hs_fd);
  gPad->Print("spect_icarus_multi.pdf");

  new TCanvas;

  std::vector<TH1*> hs_ratio;
  for(int i = 0; i < Nuniv; ++i){
    hs_ratio.push_back(Ratio(preds_fd[i]->Predict(calc),
                             preds_nd[i]->Predict(calc)).ToTH1());
  }
  WithBand(hs_ratio);
  gPad->Print("ratio_multi.pdf");

  // TODO use of true nominal makes contour disappear. Do I need to zero
  // subtract? Make stats-only contour to show, for sure.

  //  new TCanvas;

  //Define fit axes
  const FitAxis kAxSinSq2ThetaMuMu(&kFitSinSq2ThetaMuMu, 15/*40*/, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 15/*40*/, 2e-2, 1e2, true);

  // We'll call zero nominal
  calc->SetL(kBaselineSBND);
  const Spectrum data_nd = nom_nd->Predict(calc).FakeData(sbndPOT);
  calc->SetL(kBaselineIcarus);
  const Spectrum data_fd = nom_fd->Predict(calc).FakeData(icarusPOT);

  new TCanvas;

  // Nominal stats-only case
  Surface* surf_nom = 0;
  Surface* rsurf_nom = 0;
  Surface* csurf_nom = 0;
  {
    SingleSampleExperiment expt_nd(nom_nd, data_nd);
    SingleSampleExperiment expt_fd(nom_fd, data_fd);

    MultiExperimentSBN multiExpt({&expt_nd, &expt_fd}, {kSBND, kICARUS});
    //    MultiExperimentSBN multiExpt({&expt_fd}, {kICARUS});

    surf_nom = new Surface(&multiExpt, calc,
                           kAxSinSq2ThetaMuMu,
                           kAxDmSq);

    TH2* crit90 = Gaussian90Percent1D1Sided(*surf_nom);
    TH2* crit3sig = Gaussian3Sigma1D1Sided(*surf_nom);
    TH2* crit5sig = Gaussian5Sigma1D1Sided(*surf_nom);
    //    surf_nom->Draw();
    surf_nom->DrawContour(crit90, 2, kBlack);
    surf_nom->DrawContour(crit3sig, kSolid, kBlack);
    surf_nom->DrawContour(crit5sig, 7, kBlack);

    TLegend* leg = new TLegend(.15, .15, .45, .4);
    leg->SetFillStyle(0);
    TH1* dummy = new TH1F("", "", 1, 0, 1);
    dummy->SetLineStyle(2);
    leg->AddEntry(dummy->Clone(), "90% C.L.", "l");
    dummy->SetLineStyle(kSolid);
    leg->AddEntry(dummy->Clone(), "3#sigma C.L.", "l");
    dummy->SetLineStyle(7);
    leg->AddEntry(dummy->Clone(), "5#sigma C.L.", "l");
    leg->Draw();

    RatioExperiment rexp(nom_nd, nom_fd, data_nd, data_fd);
    rsurf_nom = new Surface(&rexp, calc, kAxSinSq2ThetaMuMu, kAxDmSq);

    CountingExperiment cexp(nom_nd, data_nd);
    calc->SetL(kBaselineSBND);
    csurf_nom = new Surface(&cexp, calc, kAxSinSq2ThetaMuMu, kAxDmSq);

    //    rsurf_nom->DrawContour(crit90, 2, kRed);
    rsurf_nom->DrawContour(crit3sig, kSolid, kRed);
    //    rsurf_nom->DrawContour(crit5sig, 7, kRed);

    //    csurf_nom->DrawContour(crit90, 2, kBlue);
    csurf_nom->DrawContour(crit3sig, kSolid, kBlue);
    //    csurf_nom->DrawContour(crit5sig, 7, kBlue);

    TLegend* leg2 = new TLegend(.55, .6, .85, .85);
    leg2->SetFillStyle(0);
    dummy->SetLineStyle(kSolid);
    leg2->AddEntry(dummy->Clone(), "ND+FD joint fit", "l");
    dummy->SetLineColor(kRed);
    leg2->AddEntry(dummy->Clone(), "FD / ND ratio", "l");
    dummy->SetLineColor(kBlue);
    leg2->AddEntry(dummy->Clone(), "ND norm", "l");
    leg2->Draw();

    gPad->Print("conts_nom.pdf");

    new TCanvas;
    rsurf_nom->Draw();
    gPad->SetLogz();
    gPad->Print("cont_nom_ratio_colz.pdf");
    new TCanvas;
    csurf_nom->Draw();
    gPad->SetLogz();
    gPad->Print("cont_nom_count_colz.pdf");
  }


  TH2F* hmin = 0;
  TH2F* hminmap = 0;

  Surface* surf0 = 0;

  for(int univ = 0; univ < Nuniv; ++univ){
    std::cout << "universe " << univ << std::endl;
    // SingleSampleExperiment expt_nd(preds_nd[univ], data_nd);
    // SingleSampleExperiment expt_fd(preds_fd[univ], data_fd);

    // MultiExperimentSBN multiExpt({&expt_nd, &expt_fd}, {kSBND, kICARUS});
    // //    MultiExperimentSBN multiExpt({&expt_fd}, {kICARUS});

    // Surface surf(&multiExpt, calc,
    //              kAxSinSq2ThetaMuMu,
    //              kAxDmSq);

    RatioExperiment rexp(preds_nd[univ], preds_fd[univ], data_nd, data_fd);
    Surface surf(&rexp, calc, kAxSinSq2ThetaMuMu, kAxDmSq);

    //    CountingExperiment cexp(preds_nd[univ], data_nd);
    //    calc->SetL(kBaselineSBND);
    //    Surface surf(&cexp, calc, kAxSinSq2ThetaMuMu, kAxDmSq);

    if(univ == 0){
      //      surf.Draw();
      surf0 = new Surface(surf);
      //      gPad->SetLogz();
    }
    if(hmin == 0){
      hmin = surf.ToTH2(0);
      hminmap = surf.ToTH2(0);
      hminmap->Reset(); // universe zero
    }

    TH2F* h = surf.ToTH2(0); // careful to not subtract this universe's min
    for(int i = 0; i < h->GetNbinsX()+2; ++i){
      for(int j = 0; j < h->GetNbinsY()+2; ++j){
        const double znow = hmin->GetBinContent(i, j);
        const double z = h->GetBinContent(i, j);
        if(z < znow){
          hmin->SetBinContent(i, j, z);
          hminmap->SetBinContent(i, j, univ);
        }
      }
    }
  }

  // new TCanvas;
  // hmin->Draw("colz");
  // gPad->SetLogx();
  // gPad->SetLogy();
  // gPad->SetLogz();

  Surface surfMulti = *surf0;
  surfMulti.fHist = (TH2F*)hmin->Clone();
  surfMulti.fMinChi = 0;
  surfMulti.fLogX = true;
  surfMulti.fLogY = true;

  new TCanvas;
  TH2* crit90 = Gaussian90Percent1D1Sided(*surf0);
  TH2* crit3sig = Gaussian3Sigma1D1Sided(*surf0);
  TH2* crit5sig = Gaussian5Sigma1D1Sided(*surf0);

  surfMulti.Draw();
  //  gPad->SetLogz();
  surfMulti.DrawContour(crit90, 2, kRed, 0);
  surfMulti.DrawContour(crit5sig, kSolid, kRed, 0);
  surfMulti.DrawContour(crit3sig, 7, kRed, 0);

  rsurf_nom->DrawContour(crit90, 2, kBlack, 0);
  rsurf_nom->DrawContour(crit5sig, kSolid, kBlack, 0);
  rsurf_nom->DrawContour(crit3sig, 7, kBlack, 0);

  gPad->Print("conts_multi.pdf");

  new TCanvas;
  hminmap->Draw("colz");
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->Print("minmap.pdf");
}
