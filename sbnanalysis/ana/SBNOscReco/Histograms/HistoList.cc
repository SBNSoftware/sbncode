#include "HistoList.h"
#include "TH1.h"
#include "TDirectory.h"

namespace ana {

 namespace SBNOsc {

void HistoList::Write() {
  for (unsigned i = 0; i < fAllHistos.size(); i++) {
    fLocations[i]->cd();
    fAllHistos[i]->Write();
  }
}

void HistoList::Scale(double scale) {
  for (TH1 *hist: fAllHistos) hist->Scale(scale);
}

void HistoList::StoreHisto(TH1 *hist) {
  fAllHistos.push_back(hist);
  fLocations.push_back(gDirectory);
}

void HistoList::Add(const HistoList &other) {
  for (unsigned i = 0; i < fAllHistos.size(); i++) {
    fAllHistos[i]->Add(other.fAllHistos[i]);
  }
}

void HistoList::Merge(const HistoList &merge) {
  for (unsigned i = 0; i < merge.fAllHistos.size(); i++) {
    fAllHistos.push_back(merge.fAllHistos[i]);
    fLocations.push_back(merge.fLocations[i]);
  }
}
  }
}
