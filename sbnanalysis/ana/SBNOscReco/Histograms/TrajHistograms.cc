#include "TrajHistograms.h"
#include "TH1D.h"
#include "TFile.h"

void ana::SBNOsc::TrajHistograms::Initialize() {
  for (unsigned i = 0; i < TrackHistos::nPDGs; i++) {
    dEdx[i] = new TH1D(("dEdx_" + std::string(TrackHistos::trackHistoPDGs[i])).c_str(), "dEdx", 200, 0., 20.);
  }
}

void ana::SBNOsc::TrajHistograms::Fill(const numu::RecoTrack &track, const anab::Calorimetry &collection_calo) {
  unsigned pdg_index = TrackHistos::PDGIndex(track);
  if (collection_calo.dEdx().size() == 0) return;
  // Ignore the first and last points as in the LArSoft ParticleID module
  for (unsigned i = 1; i < collection_calo.dEdx().size()-1; i++) {
    dEdx[pdg_index]->Fill(collection_calo.dEdx()[i]);
  } 
}

void ana::SBNOsc::TrajHistograms::Write() {
  for (unsigned i = 0; i < TrackHistos::nPDGs; i++) {
    dEdx[i]->Write();
  }
}

void ana::SBNOsc::TrajHistograms::Get(TFile &f) {
  for (unsigned i = 0; i < TrackHistos::nPDGs; i++) {
    dEdx[i] = (TH1D *) f.Get(("dEdx_" + std::string(TrackHistos::trackHistoPDGs[i])).c_str());
  }
}

void ana::SBNOsc::TrajHistograms::Add(const ana::SBNOsc::TrajHistograms &other) {
  for (unsigned i = 0; i < TrackHistos::nPDGs; i++) {
    dEdx[i]->Add(other.dEdx[i]);
  }
}

ana::SBNOsc::TrajHistograms::~TrajHistograms() {}

