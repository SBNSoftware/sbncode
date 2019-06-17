// This is the macro that produces Systs/flux_shifts.root"
//
// The flux files are found at
// http://home.fnal.gov/~ljf26/DUNEFluxes/
// and the covariance matrix at
// http://docs.dunescience.org:8080/cgi-bin/ShowDocument?docid=1517

#include "TRandom3.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLatex.h"

#include <iostream>

const int Nbins_numu = 19;
const double edges_numu[Nbins_numu+1] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 12, 16, 20, 40, 100};

const int Nbins_nue = 7;
const double edges_nue[Nbins_nue+1] = {0, 2, 4, 6, 8, 10, 20, 100};

// Binning for each individual spectrum
const double max_E = 8;
//const int N = 80; // bin according to finest flux bins
const int N = 16; // bin according to finest covariance bins
const double bin_width = max_E/N;

const int Nspectra = 16; // Near/Far * FHC/RHC * numu/nue * nu/nubar

void TextHelper(double x, double y, const char* txt, int col, bool horiz)
{
  TLatex* ltx = new TLatex(horiz ? y : x, horiz ? x : y, txt);
  ltx->SetTextAlign(horiz ? 23 : 21);
  if(horiz) ltx->SetTextAngle(90);
  ltx->SetNDC();
  ltx->SetTextColor(col);
  ltx->Draw();
}

void draw_edges(bool horiz = false)
{
  for(int i = 1; i < Nspectra; ++i){
    const double x = i*N;

    TGraph* g = new TGraph;
    g->SetLineColor(kRed);
    g->SetLineWidth(2);
    if(i%4 != 0) g->SetLineStyle(2);
    g->SetPoint(0, x, -1e10);
    g->SetPoint(1, x, +1e10);
    g->Draw("l same");

    if(horiz){
      g = new TGraph;
      g->SetLineColor(kRed);
      g->SetLineWidth(2);
      g->SetPoint(0, -1e10, x);
      g->SetPoint(1, +1e10, x);
      g->Draw("l same");
    }
  }

  std::vector<std::string> labels = {"#nu_{#mu}", "#bar#nu_{#mu}",
                                     "#nu_{e}", "#bar#nu_{e}"};

  for(int yaxis = false; yaxis <= horiz; ++yaxis){
    for(int i = 0; i < Nspectra; ++i){
      TLatex* ltx = 0;
      if(yaxis)
        ltx = new TLatex(0.06, 0.1 + 0.8 * (i+.5)/Nspectra, labels[i%4].c_str());
      else
        ltx = new TLatex(0.1 + 0.8 * (i+.5)/Nspectra, 0.06, labels[i%4].c_str());
      ltx->SetNDC();
      ltx->SetTextAlign(yaxis ? 12 : 21);

      if(i < 4) ltx->SetTextColor(kRed);
      if(i >= 4 && i < 8) ltx->SetTextColor(kGreen+2);
      if(i >= 8 && i < 12) ltx->SetTextColor(kBlue);
      if(i >= 12) ltx->SetTextColor(kMagenta);

      ltx->Draw();
    }

    TextHelper(.1+.125*.8, .01, "ND FHC", kRed, yaxis);
    TextHelper(.1+.375*.8, .01, "ND RHC", kGreen+2, yaxis);
    TextHelper(.1+.625*.8, .01, "FD FHC", kBlue, yaxis);
    TextHelper(.1+.875*.8, .01, "FD RHC", kMagenta, yaxis);
  }
}

double chisq(const TMatrixD& m, const TVectorD& v)
{
  return m.Similarity(v);

  //  return ROOT::Math::SimilarityT(m, v);

  //  TVectorD v_trans = v;
  //  v.Transpose;
  //  return (v_trans * m * v)(0, 0);
}

void FillInto(TH1* out, TH1* in, int i0)
{
  // Normalize here
  in->Scale(1/in->Integral(0, -1));

  for(int i = 0; i < N; ++i){
    out->Fill(i0+i+.5, in->GetBinContent(in->FindBin((i+.5)*bin_width)));
  }
}

TH1* flux_hist()
{
  TH1* ret = new TH1F("", "", Nspectra*N, 0, Nspectra*N);

  int i0 = 0;

  for(std::string hcStr: {"neutrino", "antineutrino"}){
    for(std::string detStr: {"ND", "FD"}){
      TFile* f = new TFile(("histos_g4lbne_v3r4p2_QGSP_BERT_CP_run15_12388_80GeV_"+hcStr+"_LBNE"+detStr+"_fastmc.root").c_str());

      for(std::string flavStr: {"numu", "numubar", "nue", "nuebar"}){
        TH1* h = (TH1*)f->Get((flavStr+"_flux").c_str());

        FillInto(ret, h, i0);
        i0 += N;
      }
    }
  }

  return ret;
}

int find_bin(const double* begin, const double* end, double x)
{
  if(x < *begin) return -1;
  if(x > *(end-1)) return -1;

  return std::lower_bound(begin, end, x)-begin-1;
}

int find_bin(double x, bool nue)
{
  if(nue) return find_bin(edges_nue, edges_nue +Nbins_nue, x);
  return find_bin(edges_numu, edges_numu+Nbins_numu, x);
}

TH2* rebin_cov(TH2* h)
{
  TH2* ret = new TH2F("", "", Nspectra*N, 0, Nspectra*N, Nspectra*N, 0, Nspectra*N);

  int xbin0 = 0;
  for(int i = 0; i < Nspectra; ++i){
    const bool inue = (i%4 >= 2);
    for(int ix = 0; ix < N; ++ix){
      const double x = (ix+.5)*bin_width;

      const int xbin = xbin0+find_bin(x, inue);

      if(xbin == -1) continue;

      int ybin0 = 0;
      for(int j = 0; j < Nspectra; ++j){
        const bool jnue = (j%4 >= 2);
        for(int iy = 0; iy < N; ++iy){
          const double y = (iy+.5)*bin_width;

          const int ybin = ybin0+find_bin(y, jnue);

          if(ybin == -1) continue;

          ret->Fill(i*N+ix+.5, j*N+iy+.5,
                    h->GetBinContent(xbin+1, ybin+1));
        } // end for iy
        if(jnue) ybin0 += Nbins_nue; else ybin0 += Nbins_numu;
      } // end for j
    } // end for ix
    if(inue) xbin0 += Nbins_nue; else xbin0 += Nbins_numu;
  } // end for i
  return ret;
}

TMatrixD mat(TH2* h)
{
  TMatrixD m(N*Nspectra, N*Nspectra);
  for(int i = 1; i <= N*Nspectra; ++i){
    for(int j = 1; j <= N*Nspectra; ++j){
      m(i-1, j-1) = h->GetBinContent(i, j);
    }
  }
  return m;
}

TH2* th2(TMatrixD& m)
{
  TH2* h = new TH2F("", "", N*Nspectra, 0, N*Nspectra, N*Nspectra, 0, N*Nspectra);
  for(int i = 0; i < N*Nspectra; ++i){
    for(int j = 0; j < N*Nspectra; ++j){
      h->SetBinContent(i+1, j+1, m(i, j));
    }
  }
  return h;
}

TH1* SaveEVec(TH1* h, int i)
{
  TDirectory* tmp = gDirectory;
  tmp->mkdir(TString::Format("syst%d", i).Data())->cd();

  int j = 1;

  for(int det = 0; det <= 1; ++det){
    const std::string detStr = (det == 0) ? "ND" : "FD";
    for(int hc = 0; hc <= 1; ++hc){
      const std::string hcStr = (hc == 0) ? "FHC" : "RHC";
      for(int pdg = 14; pdg >= 12; pdg -= 2){
        std::string pdgStr = (pdg == 14) ? "numu" : "nue";
        for(bool anti: {false, true}){
          if(anti) pdgStr += "bar";

          TH1* hout = new TH1F(TString::Format("%s_%s_%s",
                                               detStr.c_str(),
                                               pdgStr.c_str(),
                                               hcStr.c_str()).Data(),
                               ";True neutrino energy",
                               N, 0, max_E);

          for(int i = 0; i < hout->GetNbinsX(); ++i){
            hout->SetBinContent(i+1, h->GetBinContent(j++));
          }

          hout->Write();
        }
      }
    }
  }

  tmp->cd();
}

void diag_flux_systs()
{
  TCanvas* c = new TCanvas("a", "b", 1600, 400);
  TH1* flux = flux_hist();
  flux->GetXaxis()->SetLabelSize(0);
  flux->Draw();
  draw_edges();
  gPad->Print("plots/flux.pdf");


  c = new TCanvas("c", "d", 800, 800);
  TFile* f = new TFile("total_covariance_DUNE_opt.root");
  TH2* hcov = rebin_cov((TH2*)f->Get("total_covariance"));

  hcov->GetXaxis()->SetLabelSize(0);
  hcov->GetYaxis()->SetLabelSize(0);
  hcov->Draw("colz");
  draw_edges(true);
  gPad->SetLogz();
  gPad->Print("plots/cov_rel.pdf");

  c = new TCanvas("h", "i", 800, 800);

  // Normalize to absolute flux uncertainty
  for(int i = 1; i <= N*Nspectra; ++i){
    for(int j = 1; j <= N*Nspectra; ++j){
      const double z = hcov->GetBinContent(i, j);
      hcov->SetBinContent(i, j, z * flux->GetBinContent(i) * flux->GetBinContent(j));
    }
  }


  hcov->GetXaxis()->SetLabelSize(0);
  hcov->GetYaxis()->SetLabelSize(0);
  hcov->Draw("colz");
  draw_edges(true);
  gPad->SetLogz();
  gPad->Print("plots/cov_unfix.pdf");


  TMatrixD mcov = mat(hcov);

  std::cout << "Diagonalizing matrix to fix eigenvalues. This could take some time..." << std::endl;

  TVectorD evals;
  TMatrixD evecs = mcov.EigenVectors(evals);
  TMatrixD evalmat(evals.GetNrows(), evals.GetNrows());
  for(int i = 0; i < evals.GetNrows(); ++i){
    evalmat(i, i) = std::max(1e-14, evals[i]);
  }

  new TCanvas;

  TH1* heval = new TH1F("", ";Eigen number;Eigenvalue",
                        N*Nspectra, 0, N*Nspectra);
  for(int i = 0; i < N*Nspectra; ++i) heval->Fill(i+.5, evals[i]);
  heval->Draw();
  heval->GetYaxis()->SetRangeUser(1e-16, 1e-2);
  gPad->SetLogy();
  gPad->Print("plots/evals_unfix.pdf");

  c = new TCanvas("g", "h", 800, 800);

  //  TMatrixD evecs_inv = evecs;
  //  evecs_inv.Invert();
  // Cov was real-symmetric
  TMatrixD evecs_inv(TMatrixD::kTransposed, evecs);
  mcov = evecs*evalmat*evecs_inv;

  hcov = th2(mcov);

  hcov->GetXaxis()->SetLabelSize(0);
  hcov->GetYaxis()->SetLabelSize(0);
  hcov->Draw("colz");
  draw_edges(true);
  gPad->SetLogz();

  gPad->Print("plots/cov.pdf");

  std::cout << "Finding eigenvectors. This could take some time..." << std::endl;
  evecs = mcov.EigenVectors(evals);
  for(int i = 0; i < evals.GetNrows(); ++i){
    if(evals[i] < 0 || std::isinf(evals[i]) || std::isnan(evals[i])) std::cout << i << " " << evals[i] << std::endl;
   }

  new TCanvas;

  heval = new TH1F("", ";Eigen number;Eigenvalue",
                   N*Nspectra, 0, N*Nspectra);
  for(int i = 0; i < N*Nspectra; ++i) heval->Fill(i+.5, evals[i]);
  heval->Draw();
  heval->GetYaxis()->SetRangeUser(1e-16, 1e-2);
  gPad->SetLogy();
  gPad->Print("plots/evals.pdf");

  std::cout << "Converting eigenvectors to histograms..." << std::endl;

  std::vector<TH1*> hevecs_vec;
  std::vector<TVectorD> evecs_vec;
  for(int i = 0; i < evals.GetNrows(); ++i){
    hevecs_vec.push_back(new TH1F("", "", Nspectra*N, 0, Nspectra*N));
    evecs_vec.emplace_back(evals.GetNrows());
    for(int j = 0; j < evals.GetNrows(); ++j){
      hevecs_vec.back()->SetBinContent(j+1, sqrt(evals[i])*evecs(j, i) / flux->GetBinContent(j+1));
      evecs_vec.back()[j] = sqrt(evals[i])*evecs(j, i);
    }
  }

  c = new TCanvas("e", "f", 1600, 400);

  TFile* fout = new TFile("flux_shifts.root", "RECREATE");

  // Arbitrary number of shifts to save - in practice they go to near zero
  // around 30.
  for(int i = 0; i < N*Nspectra; ++i){
    hevecs_vec[i]->GetYaxis()->SetRangeUser(-.2, +.2);
    hevecs_vec[i]->SetTitle(TString::Format("Eigenvector %d - #sqrt{eigenvalue} = %g;;Fractional error", i, sqrt(evals[i])).Data());
    hevecs_vec[i]->GetXaxis()->SetLabelSize(0);
    hevecs_vec[i]->Draw();
    draw_edges();
    gPad->Update();
    gPad->Print(TString::Format("plots/eig%d.pdf", i).Data());

    SaveEVec(hevecs_vec[i], i);
  }

  std::cout << "Checking normalization and orthogonality..." << std::endl;

  TMatrixD mcov_inv = mcov;
  mcov_inv.Invert();

  std::cout << chisq(mcov_inv, evecs_vec[0]) << std::endl;
  std::cout << chisq(mcov_inv, evecs_vec[1]) << std::endl;
  std::cout << chisq(mcov_inv, evecs_vec[2]) << std::endl;
  std::cout << chisq(mcov_inv, evecs_vec[0]+evecs_vec[1]) << std::endl;
  std::cout << chisq(mcov_inv, evecs_vec[0]+evecs_vec[2]) << std::endl;
  std::cout << chisq(mcov_inv, evecs_vec[1]+evecs_vec[2]) << std::endl;
  std::cout << chisq(mcov_inv, evecs_vec[0]+evecs_vec[1]+evecs_vec[2]) << std::endl;
  for(int i = 0; i < 30; ++i){
    std::cout << i << " " << chisq(mcov_inv, evecs_vec[i]) << std::endl;
  }
}
