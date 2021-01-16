#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/LoadFromFile.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"

ofstream output("nuesel_icarus.txt");

namespace ana{

  // ----------------------------------------------------------------------
  // Tables

  // Highlight cell depending on range
  TString thisCellColor(double weird){
    TString mystring = "";
    double thisval = abs(weird);
    if( thisval>=0.05 && thisval<0.10) mystring = "\\cellcolor{green!25}";
    else if( thisval>=0.10 && thisval<0.15) mystring = "\\cellcolor{yellow!25}";
    else if( thisval>=0.15) mystring = "\\cellcolor{red!25}";
    return mystring;
  }

  TString fixLatexName(TString mystring){
    // latex doesnt like underscores
    std::vector<TString> in = {"#"," ",".","_"};
    std::vector<TString> out = {"","","","\\_"};

    for(unsigned int i=0;i<in.size();i++)
      mystring.ReplaceAll(in[i],out[i]);
    return mystring;
  }

  void printTableHeader(int quantId=0)
  {
    std::setprecision(3);
    output << "\\begin{table}[H]\n";
    output << "\\centering\n";
    // output << "\\resizebox{\\textwidth}{!}{\n";
    output << "\\begin{tabular}{|l||c|c|c|c|c||c|c|}\n";
    output << "\\hline \n";
    output << "\\multicolumn{1}{|c||}{} & \\multicolumn{5}{c||}{Number of interactions} & \\multicolumn{2}{c|}{Contribution} \\ \\hline \n
    \\multicolumn{1}{|c||}{Cut} & $\nu_{e}$ CC & $\nu_{\\mu}$ CC & NC & Cosmic & Other bkg & Signal & Total background \\ \\hline \n";
  }// printTableHeader

  void printTableFooter(){
    output << "\\end{tabular}}\n";
    output << "\\end{table}";
    output << "\n\n\n";
    //output.close();
  }

  void printEventsLine(std::string cutname, float nue, float numu, float nc, float cos, float other){
    float total = nue+numu+nc+cos+other;
    float percnue = nue/total;
    float perctotbkg = (numu+nc+cos+other)/total;
    output << std::fixed << std::setw(6) << std::setprecision(3) << cutname << "&" << nue << "&" << numu << "&" << nc << "&" << cos << "&" << other "&" << percnue << "&" << perctotbkg << "\\ \\hline \n";
  }


  // ----------------------------------------------------------------------
  // Canvases
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

  void FillWithDimColor(TH1* h, bool usealpha=false, float dim=0.8)
  {
    if ( usealpha ){
      h->SetFillColorAlpha(h->GetLineColor(),dim);
      return;
    }
    TColor *color = gROOT->GetColor(h->GetLineColor());
    float R,G,B,hR,hG,hB,hHue,hSat,hVal;
    color->GetRGB(hR,hG,hB);
    color->RGB2HSV(hR,hG,hB,hHue,hSat,hVal);
    color->HSV2RGB(hHue,dim*hSat,hVal,R,G,B);
    h->SetFillColor(color->GetColor(R,G,B));
  }

  // Legends and Texts
  //--------------------------------------------------
  void DrawComponentsLegend(TH1* hnue, TH1* hnumu, TH1* hnc, TH1* hcos, TH1* hother){
  // TLegend *l = new TLegend(0.60, 0.65, 0.85, 0.85, NULL,"brNDC");
    TLegend *l = new TLegend(0.70, 0.70, 0.85, 0.85, NULL,"brNDC");
    l->SetFillStyle(0);
    l->SetTextSize(0.035);
    // l->SetHeader(pot_tag.c_str());
    // l->SetHeader("Integral");
    // l->AddEntry(hnue,   Form("#nu_{e} CC:   %.2f", hnue->Integral()),   "l");
    // l->AddEntry(hnumu,  Form("#nu_{#mu} CC: %.2f", hnumu->Integral()),  "l");
    // l->AddEntry(hnc,    Form("NC:           %.2f", hnc->Integral()),    "l");
    // l->AddEntry(hcos,   Form("Cosmics:      %.2f", hcos->Integral()),   "l");
    // l->AddEntry(hother, Form("Other bkg:    %.2f", hother->Integral()), "l");
    l->AddEntry(hnue,   "#nu_{e} CC",   "l");
    l->AddEntry(hnumu,  "#nu_{#mu} CC", "l");
    l->AddEntry(hnc,    "NC",           "l");
    l->AddEntry(cos,    "Cosmics",      "l");
    l->AddEntry(hother, "Other bkg",    "l");
    l->Draw("");
  }

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


  void DrawSigBkgIntLegend(TH1* h1, char *name1, double iSig, TH1* h2, char *name2, double iBkg){

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

  void DrawSigBkgIntText(TH1* hsig, TH1* hbkg){

    float isig = hsig->Integral();
    float ibkg = hbkg->Integral();
    float psig = 100 * insig / (isig + ibkg);
    float pbkg = 100 * ibkg / (isig + ibkg);

    TPaveText *pText1 = new TPaveText(0.15, 0.78, 0.30, 0.85, "brNDC");
    TText *text1 = (pText1->AddText("6.6 #times 10^{20} POT");
    text1->SetTextSize(0.04);
    pText1->SetBorderSize(0);
    pText1->SetFillStyle(0);
    pText1->Draw();
    TPaveText *pText2 = new TPaveText(0.15, 0.65, 0.30, 0.75, "brNDC");
    TText *text2 = pText2->AddText(Form("Sig: %2.f = %2.f %%", isig, psig));
    text2->SetTextAlign(11);
    text2->SetTextSize(0.04);
    TText *text3 = pText2->AddText(Form("Bkg: %2.f = %2.f %%", ibkg, pbkg));
    text3->SetTextAlign(11);
    text3->SetTextSize(0.04);
    pText2->SetBorderSize(0);
    pText2->SetFillStyle(0);
    pText2->Draw();
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

}
