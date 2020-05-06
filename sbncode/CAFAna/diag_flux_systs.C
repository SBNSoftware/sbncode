#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionGenerator.h"

#include "CAFAna/Systs/SBNWeightSysts.h"

#include "OscLib/IOscCalculator.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLatex.h"

#include <string>

using namespace ana;

const std::string dir = "/sbnd/data/users/jlarkin/workshop_samples/";

const char* state_name = "flux_diag_state.root";

const int NEnergyBins = 70;
const double MaxEnergy = 7; // GeV

const int NSpectra = 3; // SBND/UB/Icarus

const int NUnivs = 1000;

const int kNumComponents = 100; // arbitrary

void TextHelper(double x, double y, const char* txt, int col, bool horiz)
{
  TLatex* ltx = new TLatex(horiz ? y : x, horiz ? x : y, txt);
  ltx->SetTextAlign(22);
  if(horiz) ltx->SetTextAngle(90);
  ltx->SetNDC();
  ltx->SetTextColor(col);
  ltx->Draw();
}

void draw_edges(bool horiz = false)
{
  for(int i = 1; i < NSpectra; ++i){
    const double x = i*NEnergyBins;

    TGraph* g = new TGraph;
    g->SetLineWidth(2);
    g->SetPoint(0, x, -1e10);
    g->SetPoint(1, x, +1e10);
    g->Draw("l same");

    if(horiz){
      g = new TGraph;
      g->SetLineWidth(2);
      g->SetPoint(0, -1e10, x);
      g->SetPoint(1, +1e10, x);
      g->Draw("l same");
    }
  }

  for(int yaxis = false; yaxis <= horiz; ++yaxis){
    TextHelper(.1+(1/6.)*.8, .05, "SBND", kRed, yaxis);
    TextHelper(.1+(3/6.)*.8, .05, "MicroBooNE", kGreen+2, yaxis);
    TextHelper(.1+(5/6.)*.8, .05, "Icarus", kBlue, yaxis);
  }
}

void SaveEVec(TH1* h, int i)
{
  TDirectory* tmp = gDirectory;
  tmp->mkdir(TString::Format("syst%d", i).Data())->cd();

  int j = 1;

  // These loops must match the actual order in the histogram
  for(std::string detStr: {"nd", "ub", "fd"}){
    for(std::string pdgStr: {"numu"}){//, "nue"}){
      TH1* hout = new TH1F(TString::Format("%s_%s",
                                           detStr.c_str(),
                                           pdgStr.c_str()).Data(),
                           ";True neutrino energy",
                           NEnergyBins, 0, MaxEnergy);

      for(int i = 0; i < hout->GetNbinsX(); ++i){
        hout->SetBinContent(i+1, h->GetBinContent(j++));
      }

      hout->Write();
    }
  }

  tmp->cd();
}

void diag_flux_systs(bool force_rebuild = false)
{
  // TODO does anything in here impact the flux shape?
  //  const Var kWeight = SIMPLEVAR(reco.weight);

  const Cut kOneTrue([](const caf::SRProxy* sr)
         {
           return (sr->truth.size() == 1);
         });

  const Var kTrueE = SIMPLEVAR(truth[0].neutrino.energy);

  const HistAxis axis("True neutrino energy (GeV)",
                      Binning::Simple(NEnergyBins, 0, MaxEnergy), kTrueE);

  // NB beam nues are only about 0.08% of the flux. Anti-numus are as much as 0.8%. Anti-nues basically don't exist
  const Cut kIsNumu = kOneTrue && SIMPLEVAR(truth[0].neutrino.initpdg) == 14;
  const Cut kIsNue  = kOneTrue && SIMPLEVAR(truth[0].neutrino.initpdg) == 12;

  if(force_rebuild || TFile(state_name).IsZombie()){
    TFile fout(state_name, "RECREATE");

    for(std::string det: {"nd", "fd", "ub"}){
      std::string fname = dir+"output_SBNOsc_NumuSelection_Modern_";
      if(det == "nd") fname += "SBND";
      if(det == "fd") fname += "Icarus";
      if(det == "ub") fname += "Uboone";
      fname += ".flat.root";

      SpectrumLoader loader(fname);

      Spectrum numu_nom(loader, axis, kIsNumu);
      Spectrum nue_nom(loader, axis, kIsNue);
      std::vector<Spectrum*> numus, nues;

      for(int univIdx = 0; univIdx < NUnivs; ++univIdx){
        const Var wei = GetUniverseWeight({
            "kminus_PrimaryHadronNormalization",
             "kplus_PrimaryHadronFeynmanScaling",
             "kzero_PrimaryHadronSanfordWang",
            "piplus_PrimaryHadronSWCentralSplineVariation",
           "piminus_PrimaryHadronSWCentralSplineVariation",
              }, univIdx);

        numus.push_back(new Spectrum(loader, axis, kIsNumu, kNoShift, wei));
        nues .push_back(new Spectrum(loader, axis, kIsNue,  kNoShift, wei));
      } // end for univIdx

      loader.Go();

      TDirectory* ddet = fout.mkdir(det.c_str());
      TDirectory* dnumu = ddet->mkdir("numu");
      numu_nom.SaveTo(dnumu->mkdir("nom"));
      for(int i = 0; i < NUnivs; ++i) numus[i]->SaveTo(dnumu->mkdir(TString::Format("univ_%d", i).Data()));
      TDirectory* dnue = ddet->mkdir("nue");
      nue_nom.SaveTo(dnue->mkdir("nom"));
      for(int i = 0; i < NUnivs; ++i) nues[i]->SaveTo(dnue->mkdir(TString::Format("univ_%d", i).Data()));
    } // end for det
  } // end if rebuild


  TFile fin(state_name);

  DontAddDirectory guard;

  Spectrum snom_numu_nd(*LoadFrom<Spectrum>(fin.GetDirectory("nd/numu/nom")));
  Spectrum snom_numu_fd(*LoadFrom<Spectrum>(fin.GetDirectory("fd/numu/nom")));
  Spectrum snom_numu_ub(*LoadFrom<Spectrum>(fin.GetDirectory("ub/numu/nom")));

  Spectrum snom_nue_nd(*LoadFrom<Spectrum>(fin.GetDirectory("nd/nue/nom")));
  Spectrum snom_nue_fd(*LoadFrom<Spectrum>(fin.GetDirectory("fd/nue/nom")));
  Spectrum snom_nue_ub(*LoadFrom<Spectrum>(fin.GetDirectory("ub/nue/nom")));

  TMatrixD mcov(NSpectra*NEnergyBins, NSpectra*NEnergyBins);

  TH2* hcov = new TH2D("", "",
                       NSpectra*NEnergyBins, 0, NSpectra*NEnergyBins,
                       NSpectra*NEnergyBins, 0, NSpectra*NEnergyBins);

  for(int univIdx = 0; univIdx < NUnivs; ++univIdx){
    TH1* hmnd = Ratio(*LoadFrom<Spectrum>(fin.GetDirectory(TString::Format("nd/numu/univ_%d", univIdx).Data())), snom_numu_nd).ToTH1();
    TH1* hmub = Ratio(*LoadFrom<Spectrum>(fin.GetDirectory(TString::Format("ub/numu/univ_%d", univIdx).Data())), snom_numu_ub).ToTH1();
    TH1* hmfd = Ratio(*LoadFrom<Spectrum>(fin.GetDirectory(TString::Format("fd/numu/univ_%d", univIdx).Data())), snom_numu_fd).ToTH1();

    // We don't use these for anything. They're basically insignificant, and
    // make the matrix decomposition fail.
    /*
    TH1* hend = Ratio(*LoadFrom<Spectrum>(fin.GetDirectory(TString::Format("nd/nue/univ_%d", univIdx).Data())), snom_nue_nd).ToTH1();
    TH1* heub = Ratio(*LoadFrom<Spectrum>(fin.GetDirectory(TString::Format("ub/nue/univ_%d", univIdx).Data())), snom_nue_ub).ToTH1();
    TH1* hefd = Ratio(*LoadFrom<Spectrum>(fin.GetDirectory(TString::Format("fd/nue/univ_%d", univIdx).Data())), snom_nue_fd).ToTH1();
    */

    TH1* htot = new TH1D("", "", NSpectra*NEnergyBins, 0, NSpectra*NEnergyBins);
    for(int i = 0; i < NEnergyBins; ++i){
      htot->SetBinContent(              i+1, hmnd->GetBinContent(i+1));
      htot->SetBinContent(  NEnergyBins+i+1, hmub->GetBinContent(i+1));
      htot->SetBinContent(2*NEnergyBins+i+1, hmfd->GetBinContent(i+1));
    }
    for(int i = 1; i <= NSpectra*NEnergyBins; ++i){
      if(htot->GetBinContent(i) == 0) htot->SetBinContent(i, 1);
    }
    htot->Smooth(); // this makes a big difference to the quality of the evecs
    htot->Draw(univIdx ? "hist same" : "hist");
    if(univIdx == 0){
      htot->GetXaxis()->SetLabelSize(0);
      htot->GetYaxis()->SetRangeUser(.5, 2);
    }

    for(int i = 0; i < htot->GetNbinsX(); ++i){
      for(int j = 0; j < htot->GetNbinsX(); ++j){
        // We'll do the actual mathematical work in logarithmic space, since
        // these systs presumably operate as scale factors. Note that since
        // log(1) = 0 we're baking in 1 as the true mean here.
        const double z = log(htot->GetBinContent(i+1)) * log(htot->GetBinContent(j+1));
        hcov->Fill(i+.5, j+.5, z);
        mcov(i, j) += z;
      }
    }
  }

  for(int i = 0; i < mcov.GetNrows(); ++i) mcov(i, i) += 1e-6;

  draw_edges();
  gPad->Print("plots/univs.pdf");

  new TCanvas("cov", "cov", 800, 800);

  hcov->Scale(1./NUnivs);
  mcov *= 1./NUnivs;

  hcov->GetXaxis()->SetLabelSize(0);
  hcov->GetYaxis()->SetLabelSize(0);
  hcov->Draw("col");
  draw_edges(true);
  hcov->GetZaxis()->SetRangeUser(-.01, +.01);
  gPad->Print("plots/covmat.pdf");


  TVectorD evals;
  TMatrixD evecs = mcov.EigenVectors(evals);

  new TCanvas;

  TH1* heval = new TH1F("", ";Eigen number;Eigenvalue",
                        NSpectra*NEnergyBins, 0, NSpectra*NEnergyBins);
  for(int i = 0; i < NSpectra*NEnergyBins; ++i) heval->Fill(i+.5, evals[i]);
  heval->Draw("hist");
  heval->GetYaxis()->SetRangeUser(1e-9, 1);
  gPad->SetLogy();
  gPad->Print("plots/evals.pdf");

  std::cout << "Converting eigenvectors to histograms..." << std::endl;

  std::vector<TH1*> hevecs_vec;
  for(int i = 0; i < evals.GetNrows(); ++i){
    hevecs_vec.push_back(new TH1F("", "", NSpectra*NEnergyBins, 0, NSpectra*NEnergyBins));
    for(int j = 0; j < evals.GetNrows(); ++j){
      // Convert the logarithmic eigenvalue back to a linear fraction away from 1
      const double z = exp(sqrt(evals[i])*evecs(j, i))-1;
      hevecs_vec.back()->SetBinContent(j+1, z);
    }
  }

  new TCanvas;

  TFile* fout = new TFile("flux_shifts.root", "RECREATE");

  for(int i = 0; i < kNumComponents; ++i){
    hevecs_vec[i]->GetYaxis()->SetRangeUser(-.2, +.2);
    hevecs_vec[i]->SetTitle(TString::Format("Eigenvector %d - #sqrt{eigenvalue}\
 = %g;;Fractional error", i, sqrt(evals[i])).Data());
    hevecs_vec[i]->GetXaxis()->SetLabelSize(0);
    hevecs_vec[i]->Draw();
    draw_edges();
    gPad->Update();
    gPad->Print(TString::Format("plots/eig%d.pdf", i).Data());

    SaveEVec(hevecs_vec[i], i);
  }

}



