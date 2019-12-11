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


void bars(int opt = 0) {

  string filename = "/sbnd/data/users/gputnam/NuMuReco/combined-1c/output_sbnd_combined.root";

  const unsigned int ncuts = 7;
  const unsigned int nmods = 6;

  // string cutnames[ncuts] = { "R_goodmcs",
  // 			     "R_crttrack",
  // 			     "R_crthit",
  // 			     "R_crtactive",
  // 			     "R_contained",
  // 			     "R_length" 
  // };
  			     //			     "R_flashmatch",

  string cutnames[ncuts] = { "R_length",
  			     "R_contained",
  			     "R_crtactive",
  			     "R_crthit",
  			     "R_crttrack",		 
  			     "R_goodmcs",
 			     "Reco"		
  };

  string varnames[1] = { "reco_momentum" };

  string modnames[nmods] = { "CC",
			     "CC-Other",
			     "NC-Other", 
			     "NC",
			     "InTime-Cosmic",
			     "Cosmic"};

  unsigned int col[nmods] ={kPink-5,kPink+2,kAzure+4,kAzure+5,kGray+2,kGray+1};

  TFile *f = new TFile(filename.c_str());

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

  // if(opt == 0){ // Core
  //   const Int_t nx = 8;
  //   string name[nx]   = { "Nearest Slice","CVN","Preselection",
  // 			  "p_{T}/p","Backward Photon",
  // 			  "Containment","Event Quality","Veto"};
  //   float cosm[nx] = {2.02,2.34,15734.7,71734.5,
  // 		      90340.2,94172,2880000,3260000};
  //   float beam[nx] = {12.6,12.66,161.43,329.95,347.02,355.08,545.92,569.6};
  //   float sig[nx] = {30.98,31.09,38.28,40.02,40.59,40.92,53.5,54.3};
  // }


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

  // TH1F *h2 = new TH1F("beam","beam",ncuts+2,0,ncuts+2);
  // h2->SetFillColor(kViolet+6);
  // h2->SetLineColor(kViolet+6);
  // h2->SetBarWidth(0.25);
  // h2->SetBarOffset(0.375);
  // h2->SetStats(0);
  // for (unsigned int i=1;i<=ncuts;i++) 
  //   h2->SetBinContent(i+1, evtcounts[i-1][1]);
  // h2->Draw("hbar0 same");

  // TH1F *h3 = new TH1F("cosmic","cosmic",ncuts+2,0,ncuts+2);
  // h3->SetFillColor(kAzure+2);
  // h3->SetLineColor(kAzure+2);
  // h3->SetBarWidth(0.25);
  // h3->SetBarOffset(0.125);
  // h3->SetStats(0);
  // for (unsigned int i=1;i<=ncuts;i++) 
  //   h3->SetBinContent(i+1, evtcounts[i-1][2]);
  // h3->Draw("hbar0 same");

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


  TLegend *l = new TLegend(0.807114,0.165633,0.964429,0.280779);
  l->AddEntry(h1, modnames[0].c_str(),"f");
  for ( unsigned int iMod = 1; iMod < nmods; iMod++ ){
    l->AddEntry(h[iMod], modnames[iMod].c_str(),"f");
  }
  // l->AddEntry(h2,"Beam Background","f");
  // l->AddEntry(h3,"Cosmic Background","f");
  // l->AddEntry(h[3], modnames[3].c_str(),"f");
  // l->AddEntry(h[4], modnames[4].c_str(),"f");
  // l->AddEntry(h[5], modnames[5].c_str(),"f");
  l->SetLineWidth(0);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetTextSize(0.0159433);
  l->Draw();

  TLatex *pot = new TLatex(0.925,0.88,"XXX POT normalized");
  pot->SetNDC();
  pot->SetTextSize(1/50.);
  pot->SetTextAlign(32);
  pot->Draw();

  TLatex *time = new TLatex(0.925,0.86,"XXX seconds livetime");
  time->SetNDC();
  time->SetTextSize(1/50.);
  time->SetTextAlign(32);
  time->Draw();

  TLatex* prelim = new TLatex(.95, .925, "SBN Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/55.);
  prelim->SetTextAlign(32);
  prelim->Draw();

  if(opt == 0) c->Print("cutFlow_SBND.pdf");
  if(opt == 1) c->Print("cutFlow_ICARUS.pdf");

  TLatex *sample;
  if(opt==0 ) sample = new TLatex(0.25,0.86,"SBND Sample");
  if(opt==1 ) sample = new TLatex(0.25,0.86,"ICARUS Sample");
  sample->SetNDC();
  sample->SetTextSize(1/40.);
  sample->SetTextAlign(12);
  sample->Draw();

  if(opt == 0) c->Print("cutFlow_SBND_label.pdf");
  if(opt == 1) c->Print("cutFlow_ICARUS_label.pdf");

}
