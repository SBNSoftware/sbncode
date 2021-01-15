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