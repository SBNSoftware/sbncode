// This makes xsec_shifts.root, which nothing uses. For now it's just for fun.
//
// Input covariance matrix comes from
// /dune/data/users/marshalc/total_covariance_XS.root (I know...)

#include "TRandom3.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLatex.h"
#include "TROOT.h"

#include <iostream>

void label_bins(TH1* h, bool horiz = false)
{
  const std::string valor_categories[33] = { "#nu CCQE 1",       // 0
                                             "#nu CCQE 2",       // 1
                                             "#nu CCQE 3",       // 2
                                             "#bar{#nu} CCQE 1",    // 3
                                             "#bar{#nu} CCQE 2",    // 4
                                             "#bar{#nu} CCQE 3",    // 5
                                             "#nu MEC",    // 6
                                             "#bar{#nu} MEC", // 7
                                             "#nu CC 1#pi^{0} 1",     // 8
                                             "#nu CC 1#pi^{0} 2",     // 9
                                             "#nu CC 1#pi^{0} 3",     // 10
                                             "#nu CC 1#pi^{#pm} 1",     // 11
                                             "#nu CC 1#pi^{#pm} 2",     // 12
                                             "#nu CC 1#pi^{#pm} 3",     // 13
                                             "#bar{#nu} CC 1#pi^{0} 1",  // 14
                                             "#bar{#nu} CC 1#pi^{0} 2",  // 15
                                             "#bar{#nu} CC 1#pi^{0} 3",  // 16
                                             "#bar{#nu} CC 1#pi^{#pm} 1",  // 17
                                             "#bar{#nu} CC 1#pi^{#pm} 2",  // 18
                                             "#bar{#nu} CC 1#pi^{#pm} 3",  // 19
                                             "#nu 2#pi",          // 20
                                             "#bar{#nu} 2#pi",       // 21
                                             "#nu DIS 1",        // 22
                                             "#nu DIS 2",        // 23
                                             "#nu DIS 3",        // 24
                                             "#bar{#nu} DIS 1",     // 25
                                             "#bar{#nu} DIS 2",     // 26
                                             "#bar{#nu} DIS 3",     // 27
                                             "#nu COH",          // 28
                                             "#bar{#nu} COH",       // 29
                                             "#nu NC",           // 30
                                             "#bar{#nu} NC",        // 31
                                             "#nu_{e}/#nu_{#mu}" }; // 32

  for(int yaxis = false; yaxis <= horiz; ++yaxis){
    TAxis* ax = yaxis ? h->GetYaxis() : h->GetXaxis();

    for(int i = 0; i <= 32; ++i){
      ax->SetBinLabel(i+1, valor_categories[i].c_str());
    }
    // Otherwise everything disappears in TH2
    for(int i = 33; i <= 42; ++i){
      ax->SetBinLabel(i+1, " ");
    }

    ax->LabelsOption(yaxis ? "h" : "v");
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

TMatrixD mat(TH2* h)
{
  const int N = h->GetNbinsX();
  TMatrixD m(N, N);
  for(int i = 1; i <= N; ++i){
    for(int j = 1; j <= N; ++j){
      m(i-1, j-1) = h->GetBinContent(i, j);
    }
  }
  return m;
}

TH2* th2(TMatrixD& m)
{
  const int N = m.GetNrows();
  TH2* h = new TH2F("", "", N, 0, N, N, 0, N);
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < N; ++j){
      h->SetBinContent(i+1, j+1, m(i, j));
    }
  }
  return h;
}

void diag_xsec_systs()
{
  gROOT->ForceStyle();

  TCanvas* c = new TCanvas("c", "d", 800, 800);
  TFile* f = new TFile("total_covariance_XS.root");
  TH2* hcov = (TH2*)f->Get("xs_covariance");

  const int N = hcov->GetNbinsX();

  hcov->Draw("colz");

  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.15);
  label_bins(hcov, true);
  gPad->SetLogz();
  hcov->SetMinimum(1e-4);
  gPad->Print("plots/xsec_cov_rel.pdf");


  TMatrixD mcov = mat(hcov);

  std::cout << "Diagonalizing matrix to fix eigenvalues. This could take some time..." << std::endl;

  TVectorD evals;
  TMatrixD evecs = mcov.EigenVectors(evals);
  TMatrixD evalmat(N, N);
  for(int i = 0; i < N; ++i){
    evalmat(i, i) = std::max(1e-14, evals[i]);
    std::cout << i << " " << evals[i] << std::endl;
  }

  new TCanvas;

  TH1* heval = new TH1F("", ";Eigen number;Eigenvalue", N, 0, N);
  for(int i = 0; i < N; ++i) heval->Fill(i+.5, evals[i]);
  heval->Draw();
  heval->GetYaxis()->SetRangeUser(1e-7, 10);
  gPad->SetLogy();
  gPad->Print("plots/xsec_evals.pdf");

  //  TMatrixD evecs_inv = evecs;
  //  evecs_inv.Invert();
  // Cov was real-symmetric
  TMatrixD evecs_inv(TMatrixD::kTransposed, evecs);
  mcov = evecs*evalmat*evecs_inv;

  hcov = th2(mcov);

  std::cout << "Finding eigenvectors. This could take some time..." << std::endl;
  evecs = mcov.EigenVectors(evals);
  for(int i = 0; i < N; ++i){
    if(evals[i] < 0 || std::isinf(evals[i]) || std::isnan(evals[i])) std::cout << i << " " << evals[i] << std::endl;
   }

  std::cout << "Converting eigenvectors to histograms..." << std::endl;

  std::vector<TH1*> hevecs_vec;
  std::vector<TVectorD> evecs_vec;
  for(int i = 0; i < N; ++i){
    hevecs_vec.push_back(new TH1F("", "", N, 0, N));
    evecs_vec.emplace_back(N);
    for(int j = 0; j < N; ++j){
      hevecs_vec.back()->SetBinContent(j+1, sqrt(evals[i])*evecs(j, i));
      evecs_vec.back()[j] = sqrt(evals[i])*evecs(j, i);
    }
  }

  new TCanvas;

  TFile* fout = new TFile("xsec_shifts.root", "RECREATE");

  for(int i = 0; i < N; ++i){
    hevecs_vec[i]->GetYaxis()->SetRangeUser(-1, +1);
    hevecs_vec[i]->SetTitle(TString::Format("Eigenvector %d - #sqrt{eigenvalue} = %g;;Fractional error", i, sqrt(evals[i])).Data());

    hevecs_vec[i]->Draw();
    gPad->SetBottomMargin(0.15);
    label_bins(hevecs_vec[i]);

    gPad->Update();
    gPad->Print(TString::Format("plots/xsec_eig%d.pdf", i).Data());

    hevecs_vec[i]->Write(TString::Format("syst%d", i).Data());
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
