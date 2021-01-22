#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Spectrum.h"

#include "helper_pur_slc_cumulative.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>

using namespace ana;

void plot_spectra_pur_slc_cumulative(std::string name = "nue")
{
  // Arbitrary POT to scale all plots to
  // Spill tree currently not being filled, so there will be a warning
  std::string inFile = "spectra_pur_slc_cumulative_" + name + ".root";
  double POT = 6.6E20;

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  const unsigned int kNType = types.size();

  std::ofstream texEff, texEffRel, texNum, texHdr;

  texEff.open("tex/" + name + "_pur_slc_cumulative_Eff.tex");
  texEffRel.open("tex/" + name + "_pur_slc_cumulative_EffRel.tex");
  texNum.open("tex/" + name + "_pur_slc_cumulative_Num.tex");
  texHdr.open("tex/" + name + "_pur_slc_cumulative_Hdr.tex");

  texEff << "True Type: ";
  texEffRel << "True Type: ";
  texNum << "True Type: ";
  texHdr << "True Type: ";

  texNum.imbue(std::locale(""));

  std::vector<float> baseEvents;
  for (unsigned int lType = 0; lType < kNType; lType++) {
    std::string thissuffix = plots[0].suffix + "_" + sels[0].suffix + "_" + types[lType].suffix;
    Spectrum* spec = LoadFromFile<Spectrum>(inFile, thissuffix).release();

    const double thisPOT(spec->POT());
    if (thisPOT < std::numeric_limits<double>::epsilon())
      spec->OverridePOT(POT);

    baseEvents.push_back(spec->Integral(POT));
    std::cout << "TEST: " << types[lType].label.c_str() << ": " << spec->Integral(POT) << std::endl;

    texEff << " & " << types[lType].label.c_str();
    texEffRel << " & " << types[lType].label.c_str();
    texNum << " & " << types[lType].label.c_str();
    texHdr << " & " << types[lType].label.c_str();
  }
  std::vector<float> prevEvents = baseEvents;

  texEff << " \\\\\n\\hline \\hline" << std::endl;
  texEffRel << " \\\\\n\\hline \\hline" << std::endl;
  texNum << " \\\\\n\\hline \\hline" << std::endl;
  texHdr << " \\\\\n\\hline \\hline" << std::endl;

  TFile* newFile = new TFile(("plots/plots_pur_slc_cumulative_" + name + ".root").c_str(), "RECREATE");

  // I want to make a plot for each var
  for (unsigned int iVar = 0; iVar < kNVar; ++iVar) {

    for (unsigned int jSel = 0; jSel < kNSel; ++jSel) {

      std::string mysuffix = plots[iVar].suffix + "_" + sels[jSel].suffix;
      TCanvas* c = new TCanvas(("c" + mysuffix).c_str(), mysuffix.c_str());

      const std::string myTitle(plots[iVar].label + ": " + sels[jSel].label + ";" + plots[iVar].label);

      THStack* stack = new THStack(mysuffix.c_str(), myTitle.c_str());
      std::vector<TH1*> hists;
      TLegend* l = new TLegend(0.15, 0.6, 0.4, 0.8);
      l->SetFillStyle(0);

      std::ofstream thisTexEff, thisTexEffRel, thisTexNum;
      thisTexNum.imbue(std::locale(""));

      if (iVar == 0) {
        std::cout << "Cut: " << sels[jSel].suffix << std::endl;

        texEff << sels[jSel].label;
        texEffRel << sels[jSel].label;
        texNum << sels[jSel].label;

        thisTexEff.open("tex/cuts_" + name + "_" + sels[jSel].suffix + "_pur_slc_cumulative_Eff.tex");
        thisTexEffRel.open("tex/cuts_" + name + "_" + sels[jSel].suffix + "_pur_slc_cumulative_EffRel.tex");
        thisTexNum.open("tex/cuts_" + name + "_" + sels[jSel].suffix + "_pur_slc_cumulative_Num.tex");

        thisTexEff << sels[jSel].label;
        thisTexEffRel << sels[jSel].label;
        thisTexNum << sels[jSel].label;
      }

      for (unsigned int lType = 0; lType < kNType; lType++) {
        std::string thissuffix = mysuffix + "_" + types[lType].suffix;
        Spectrum* spec = LoadFromFile<Spectrum>(inFile, thissuffix).release();

        const double thisPOT(spec->POT());
        if (thisPOT < std::numeric_limits<double>::epsilon())
          spec->OverridePOT(POT);

        TH1* h = spec->ToTH1(POT);
        h->SetLineColor(types[lType].color);

        hists.push_back(h);
        l->AddEntry(h, types[lType].label.c_str(), "l");
        stack->Add(hists.at(lType));

        newFile->cd();
        h->SetName(("h" + thissuffix).c_str());
        h->Write();

        if (iVar == 0) {
          std::cout << "TEST: " << types[lType].label.c_str() << ": " << spec->Integral(POT) << " (" << spec->Integral(POT) / baseEvents.at(lType) << ")" << std::endl;

          texEff << " & " << std::setprecision(2) << std::fixed << 100 * spec->Integral(POT) / baseEvents.at(lType);
          texEffRel << " & " << std::setprecision(2) << std::fixed << 100 * spec->Integral(POT) / prevEvents.at(lType);
          texNum << " & " << std::setprecision(0) << std::fixed << spec->Integral(POT);

          thisTexEff << " & " << std::setprecision(2) << std::fixed << 100 * spec->Integral(POT) / baseEvents.at(lType);
          thisTexEffRel << " & " << std::setprecision(2) << std::fixed << 100 * spec->Integral(POT) / prevEvents.at(lType);
          thisTexNum << " & " << std::setprecision(0) << std::fixed << spec->Integral(POT);

          prevEvents.at(lType) = spec->Integral(POT);
        }
      }

      if (iVar == 0) {
        texEff << " \\\\\n\\hline" << std::endl;
        texEffRel << " \\\\\n\\hline" << std::endl;
        texNum << " \\\\\n\\hline" << std::endl;

        thisTexEff << "\\\\\n\\hline" << std::endl;
        thisTexEffRel << "\\\\\n\\hline" << std::endl;
        thisTexNum << "\\\\\n\\hline" << std::endl;
      }

      stack->Draw("hist nostack e");
      l->Draw();

      c->Print(("plots/pur_slc_cumulative_" + plots[iVar].suffix + "_" + sels[jSel].suffix + "_" + name + ".png").c_str());

      newFile->cd();
      stack->Write();
      c->Write();
      delete c;
    } // iSel
  }   // iVar
  newFile->Close();
}
