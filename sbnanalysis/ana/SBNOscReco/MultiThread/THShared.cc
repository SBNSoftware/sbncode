#include "THShared.h"

TH1Shared::TH1Shared(TH1 *_hist): 
  hist(_hist) {}

void TH1Shared::Fill(double x, double scale) {
  HistoManager::Fill(Get(), x, scale);
}

const TH1 *TH1Shared::Get() const {
  return hist.get();
}

TH1 *TH1Shared::Get() {
  return hist.get();
}

TH2Shared::TH2Shared(TH2 *_hist): 
  hist(_hist) {}

void TH2Shared::Fill(double x, double y, double scale) {
  HistoManager::Fill(Get(), x, y, scale);
}

TH2 *TH2Shared::Get() {
  return hist.get();
}

const TH2 *TH2Shared::Get() const {
  return hist.get();
}

TH3Shared::TH3Shared(TH3 *_hist): 
  hist(_hist) {}

void TH3Shared::Fill(double x, double y, double z, double scale) {
  HistoManager::Fill(Get(), x, y, z, scale);
}

TH3 *TH3Shared::Get() {
  return hist.get();
}

const TH3 *TH3Shared::Get() const {
  return hist.get();
}
