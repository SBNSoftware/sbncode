////////////////////////////////////////////////
//  bars.C  - Makes cut-flow bar plots 
//
////////////////////////////////////////////////
#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

void bars(int opt = 1) {

  string filename[2] = { "/sbnd/data/users/gputnam/NuMuReco/combined-1c/output_sbnd_combined_old.root",
			 "/sbnd/data/users/gputnam/NuMuReco/combined-1c/output_icarus_combined_old.root" };

  const unsigned int ncuts = 8;
  const unsigned int nmods = 6;

  string cutnames[ncuts] = { "R_length",
  			     "R_contained",
  			     "R_crtactive",
  			     "R_crthit",
  			     "R_crttrack",		 
  			     "R_fid",
  			     "R_goodmcs",
 			     "Reco" };
  string modnames[nmods] = { "CC",
			     "NC",
			     "NC-Other", 
			     "CC-Other",
			     "InTime-Cosmic",
			     "Cosmic" };
  string varnames[1]     = { "reco_momentum" };
  //  unsigned int col[nmods]= { kPink-5,kPink+2,kAzure+4,kAzure+5,kGray+2,kGray+1 };
  unsigned int col[nmods]= { kPink+2,kAzure+3,kAzure+4,kAzure+5,kGray+2,kGray+1 };

  TFile *f = new TFile(filename[opt].c_str());

  // Count Events
  int evtcounts[ncuts][nmods];
  for ( unsigned int iCut = 0; iCut < ncuts; iCut++ ){
    for ( unsigned int iMod = 0; iMod < nmods; iMod++ ){

      string pre = "_Primary_";
      string post = "_all_all__";
      string hname= varnames[0] + pre + modnames[iMod] + post+ cutnames[iCut];
      TH1F * hist = (TH1F*)f->Get(hname.c_str());

      evtcounts[iCut][iMod] = hist->Integral();
      std::cout << hname << std::endl;
    }
  }

  // Draw Bars 
  TCanvas *c = new TCanvas("Cut efficiency","c",1000,1200);
  c->SetLogx();
  c->SetGridx();
  c->SetLeftMargin(0.2);
  c->SetRightMargin(0.05);

  TH1F *h1 = new TH1F("signal","",ncuts+2,0,ncuts+2);
  h1->SetFillColor(col[0]);
  h1->SetLineColor(col[0]);
  h1->SetBarWidth(0.14);
  //  h1->SetBarOffset(0.654);
  h1->SetBarOffset(0.85);
  h1->GetXaxis()->SetLabelSize(0.045);
  h1->GetXaxis()->SetTickLength(0);
  h1->GetYaxis()->SetTitle("Events");
  h1->GetYaxis()->CenterTitle();
  h1->SetStats(0);
  h1->SetMinimum(1);
  h1->SetMaximum(10000000);
  for(unsigned int i=1; i<=ncuts; i++) {
    h1->SetBinContent(i+1, evtcounts[i-1][0]);
    h1->GetXaxis()->SetBinLabel(i+1,cutnames[i-1].c_str());
  }
  h1->Draw("hbar0");
  h1->Draw("hbar0 same");

  TH1F *h[nmods];
  for ( unsigned int iMod = 1; iMod < nmods; iMod++ ){

    h[iMod] = new TH1F( modnames[iMod].c_str(),
			modnames[iMod].c_str(), ncuts+2, 0, ncuts+2);
    
    h[iMod]->SetFillColor(col[iMod]);
    h[iMod]->SetLineColor(col[iMod]);
    h[iMod]->SetBarWidth(0.14);
    h[iMod]->SetBarOffset(0.85-(iMod*0.14));
    h[iMod]->SetStats(0);
    for (unsigned int i=1;i<=ncuts;i++) 
      h[iMod]->SetBinContent(i+1, evtcounts[i-1][iMod]);
    h[iMod]->Draw("hbar0 same");

  }

  // Add Legends and stuff 
  TLegend *l = new TLegend(0.807114,0.165633,0.964429,0.280779);
  l->AddEntry(h1, modnames[0].c_str(),"f");
  for ( unsigned int iMod = 1; iMod < nmods; iMod++ ){
    l->AddEntry(h[iMod], modnames[iMod].c_str(),"f");
  }
  l->SetLineWidth(0);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetTextSize(0.0159433);
  l->Draw();

  TLatex *pot = new TLatex(0.925,0.88,"1x10^{20} POT normalized");
  pot->SetNDC();
  pot->SetTextSize(1/50.);
  pot->SetTextAlign(32);
  pot->Draw();

  TLatex *time = new TLatex(0.925,0.86,"XXX seconds livetime");
  time->SetNDC();
  time->SetTextSize(1/50.);
  time->SetTextAlign(32);
  //  time->Draw();

  TLatex* prelim = new TLatex(.95, .92, "SBN Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/55.);
  prelim->SetTextAlign(32);
  prelim->Draw();

  if(opt == 0) c->Print("cutFlow_SBND.pdf");
  if(opt == 1) c->Print("cutFlow_ICARUS.pdf");

  TLatex *sample;
  if(opt==0 ) sample = new TLatex(0.22,0.91,"SBND Sample");
  if(opt==1 ) sample = new TLatex(0.22,0.91,"ICARUS Sample");
  sample->SetNDC();
  sample->SetTextSize(1/40.);
  sample->SetTextAlign(12);
  sample->Draw();

  if(opt == 0) c->Print("cutFlow_SBND_label.pdf");
  if(opt == 1) c->Print("cutFlow_ICARUS_label.pdf");

}
