// This is a straight ROOT macro (ie no need to run under cafe)

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPad.h"

#include <string>

void plot_fitter_validation()
{
  TFile* fin = new TFile("fitter_validation_cafana.root");

  const int dcpCols[] = {kBlack, kRed, kGreen+2, kBlue};
  const std::string dcpStrs[] = {"0", "#pi/2", "#pi", "3#pi/2"};

  for(std::string hcStr: {"fhc", "rhc"}){
    const std::string HCStr = (hcStr == "fhc") ? "FHC" : "RHC";

    for(std::string chanStr: {"numu", "nue"}){
      const std::string CHANStr = (chanStr == "numu") ? "#nu_{#mu}" : "#nu_{e}";

      for(std::string hieStr: {"nh", "ih"}){
        const std::string HIEStr = (hieStr == "nh") ? "NH" : "IH";

        new TCanvas;

        for(int deltaIdx2 = 0; deltaIdx2 < 4; ++deltaIdx2){
          // For neutrinos 3pi/2 is the tallest histogram, draw it first, for
          // antineutrinos we need pi/2 first.
          const int deltaIdx = (hcStr == "fhc") ? 3-deltaIdx2 : (deltaIdx2+1)%4;

          const std::string dcpStr = TString::Format("%gpi", deltaIdx/2.).Data();

          TH1* h = (TH1*)fin->Get((chanStr+"_"+hcStr+"_"+hieStr+"_"+dcpStr).c_str());
          h->SetLineColor(dcpCols[deltaIdx]);
          h->Draw("same");

          h->SetTitle(("5 yrs "+HCStr+" "+CHANStr+" "+HIEStr).c_str());
        } // end for deltaIdx

        TLegend* leg = new TLegend(.6, .6, .85, .85);
        leg->SetFillStyle(0);
        for(int deltaIdx = 0; deltaIdx < 4; ++deltaIdx){
          TH1* dummy = new TH1F("", "", 1, 0, 1);
          dummy->SetLineColor(dcpCols[deltaIdx]);
          leg->AddEntry(dummy, ("#delta_{CP}="+dcpStrs[deltaIdx]).c_str(), "l");
        }
        leg->Draw("same");

        gPad->Print((hcStr+"_"+chanStr+"_"+hieStr+".pdf").c_str());
      } // end for hieStr
    } // end for chanStr
  } // end for hcStr


  TGraph* gNH = (TGraph*)fin->Get("sens_nh");
  TGraph* gIH = (TGraph*)fin->Get("sens_ih");
  TGraph* gNHOscErr = (TGraph*)fin->Get("sens_nh_oscerr");
  TGraph* gIHOscErr = (TGraph*)fin->Get("sens_ih_oscerr");

  TGraph* gNHFlux[10];
  TGraph* gIHFlux[10];
  for(int i = 0; i < 10; ++i){
    gNHFlux[i] = (TGraph*)fin->Get(TString::Format("sens_nh_flux%d", i).Data());
    gIHFlux[i] = (TGraph*)fin->Get(TString::Format("sens_ih_flux%d", i).Data());
  }

  TGraph* gNHXSec[10];
  TGraph* gIHXSec[10];
  for(int i = 0; i < 10; ++i){
    gNHXSec[i] = (TGraph*)fin->Get(TString::Format("sens_nh_xsec%d", i).Data());
    gIHXSec[i] = (TGraph*)fin->Get(TString::Format("sens_ih_xsec%d", i).Data());
  }

  TH2* axes = new TH2F("", ";#delta_{CP} / #pi;#sigma = #sqrt{#Delta#chi^{2}}", 100, 0, 2, 100, 0, 8);
  axes->GetXaxis()->CenterTitle();
  axes->GetYaxis()->SetTitleOffset(.75);
  axes->GetYaxis()->CenterTitle();
  axes->Draw();

  gNH->Draw("l same");
  gIH->Draw("l same");
  gNHOscErr->Draw("l same");
  gIHOscErr->Draw("l same");

  TLegend* leg = new TLegend(.4, .65, .6, .875);
  leg->SetFillStyle(0);
  leg->AddEntry(gNH, "NH", "l");
  leg->AddEntry(gIH, "IH", "l");
  leg->AddEntry(gNHOscErr, "NH osc err", "l");
  leg->AddEntry(gIHOscErr, "IH osc err", "l");
  leg->Draw();

  gPad->Print("mcd.pdf");


  new TCanvas;
  axes->Draw();
  for(int i = 0; i < 10; ++i){
    gNHFlux[i]->Draw("l same");
    gIHFlux[i]->Draw("l same");
  }
  gNH->Draw("l same");
  gIH->Draw("l same");

  leg = new TLegend(.4, .65, .6, .875);
  leg->SetFillStyle(0);
  leg->AddEntry(gNH, "NH", "l");
  leg->AddEntry(gIH, "IH", "l");
  leg->AddEntry(gNHFlux[0], "NH flux err", "l");
  leg->AddEntry(gIHFlux[0], "IH flux err", "l");
  leg->Draw();

  gPad->Print("mcd_flux.pdf");


  new TCanvas;
  axes->Draw();
  for(int i = 0; i < 10; ++i){
    gNHXSec[i]->Draw("l same");
    gIHXSec[i]->Draw("l same");
  }
  gNH->Draw("l same");
  gIH->Draw("l same");

  leg = new TLegend(.4, .65, .6, .875);
  leg->SetFillStyle(0);
  leg->AddEntry(gNH, "NH", "l");
  leg->AddEntry(gIH, "IH", "l");
  leg->AddEntry(gNHXSec[0], "NH xsec err", "l");
  leg->AddEntry(gIHXSec[0], "IH xsec err", "l");
  leg->Draw();

  gPad->Print("mcd_xsec.pdf");
}
