#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/LoadFromFile.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"

//----------------------------------------------------------------------
// split canvas in 2
void SplitCanvas2(TCanvas *& c1, TPad *& pad1, TPad *& pad2){

  c1 = new TCanvas("c1","",500,700);
  c1->cd();

  pad1 = new TPad("pad1","pad1",0,0,1,1);
  pad1->Draw();
  pad1->SetTopMargin(0.1);
  pad1->SetBottomMargin(0.4);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.03);
  pad1->SetFillStyle(0);
  c1->cd();

  pad2 = new TPad("pad2","pad2",0,0,1,1);// x1 y1 x2 y2
  pad2->Draw();
  pad2->SetTopMargin(0.6);
  pad2->SetBottomMargin(0.1);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.03);
  pad2->SetFillStyle(0);
  c1->cd();

}

float GetHistMax(std::vector<TH1*> histos){

  float hmax = 0.;
  for(unsigned int hId=0; hId<histos.size(); hId++){
    float thismax = histos[hId]->GetMaximum();
    if(thismax>hmax) hmax=thismax;
  }
  return hmax;
}


void PimpHist(TH1* histo, Color_t color, Style_t linestyle, int linewidth, Style_t markerstyle=8, double markersize=8){

  histo->SetLineColor(color);
  histo->SetLineStyle(linestyle);
  histo->SetLineWidth(linewidth);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(markerstyle);
  histo->SetMarkerSize(markersize);

}


// Legends
//--------------------------------------------------
void DrawSigBkgLegend(TH1* h1, char *name1, TH1* h2, char *name2){

  TLegend *leg = new TLegend(.60,.60,.8,.8);
  leg->AddEntry(h1, name1,"l");
  leg->AddEntry(h2, name2,"l");
  leg->SetBorderSize(0); //no border for legend
  leg->SetFillColor(0);  //fill colour is white
  leg->SetFillStyle(0);  //fill colour is white
  leg->SetTextSize(0.04);
  leg->Draw();

}


void DrawSigBkgIntLegend(TH1* h1, char *name1, double iSig, TH1* h2, char *name2, double iBac){

  //double iRatio = iBac/iSig;

  TLegend *leg = new TLegend(.60,.60,.8,.8);
  leg->AddEntry(h1, name1,"l");
  leg->AddEntry(h1, TString::Format("%.2f",iSig),"");
  leg->AddEntry(h2, name2,"l");
  leg->AddEntry(h2, TString::Format("%.2f",iBac),"");
  //leg->AddEntry(h2, TString::Format("Ratio=%.2f",iRatio),"");
  leg->SetBorderSize(0); //no border for legend                                                                                                                                                                                                                             
  leg->SetFillColor(0);  //fill colour is white                                                                                                                                                                                                                             
  leg->SetFillStyle(0);  //fill colour is white                                                                                                                                                                                                                             
  leg->SetTextSize(0.04);
  leg->Draw();

}

//--------------------------------------------------
void DrawEffPurLegend(TGraph* g1, char *name1, TGraph* g2, char *name2){

  TLegend *leg = new TLegend(.6,.62,.8,.74);
  leg->AddEntry(g1, name1,"l");
  leg->AddEntry(g2, name2,"l");
  leg->SetBorderSize(0); //no border for legend
  leg->SetFillColor(0);  //fill colour is white
  leg->SetFillStyle(0);  //fill colour is white
  leg->SetTextSize(0.04);
  leg->Draw();

}

//--------------------------------------------------
void DrawEffPurSMLegend(TGraph* g1, char *name1, TGraph* g2, char *name2){

  TLegend *leg = new TLegend(.69,.61,.84,.81);
  leg->AddEntry(g1, name1,"l");
  leg->AddEntry(g2, name2,"l");
  leg->SetBorderSize(0); //no border for legend
  leg->SetFillColor(0);  //fill colour is white
  leg->SetFillStyle(0);  //fill colour is white
  leg->SetTextSize(0.06);
  leg->Draw();

}

// Efficiency and purity graphs 
//-----------------------------------------------------------------

TGraph* SelEFForPURvsX(TH1* hSelSignal, TH1* hSelBack, TH1* hSignal, bool effpur) {

  //
  // Make a ROC TGraph for the given signal and background histos
  //

  const int NBins = hSignal->GetNbinsX();

  TString xTitle = hSignal->GetXaxis()->GetTitle();
  TString yTitle = hSignal->GetYaxis()->GetTitle();

  double eff[NBins], pur[NBins], val[NBins];
  double sb[NBins], ssb[NBins];

  // loop over bins to calculate eff and pur as functions of the bin
  for(unsigned int i = 1; i <= (unsigned int)NBins; ++i) {
    double allsig = hSignal   ->GetBinContent(i);
    double selsig = hSelSignal->GetBinContent(i);
    double selbac = hSelBack  ->GetBinContent(i);
    /* double S = hSignal->Integral(i,NBins); */
    /* double B = hBack  ->Integral(i,NBins); */

    val[i-1] = hSignal->GetBinCenter(i);
    double EFF = selsig;
    double PUR = selsig;

    if ( (selsig + selbac) > 0 ){
      if ( allsig > 0 ) EFF = EFF/allsig;
      PUR = PUR / (selsig + selbac);
    }

    eff[i-1] = EFF;
    pur[i-1] = PUR;
    if ( PUR > 1 || EFF > 1 ){
      std::cout<< " >>> \n >>> \n>>> \n>>> \n>>> \n>>> \n>>> \n" << " >>> EFF "<< EFF << "\t PUR "<< PUR << std::endl;
    }
    
  }

  TString n = hSignal->GetName();
  TGraph *graphPur= new TGraph(NBins,val,pur);
  graphPur->SetName("SelPur_"+n);
  graphPur->GetXaxis()->SetTitle(xTitle);
  graphPur->GetYaxis()->SetTitle("Pur.");
  TGraph *graphEff= new TGraph(NBins,val,eff);
  graphEff->SetName("SelEff_"+n);
  graphEff->GetXaxis()->SetTitle(xTitle);
  graphEff->GetYaxis()->SetTitle("Eff.");

  if ( effpur == 0 )
    return graphEff;
  else
    return graphPur;

}
