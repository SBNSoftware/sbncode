#include "TrajHistograms.h"
#include "TH1D.h"
#include "TFile.h"

void ana::SBNOsc::TrajHistograms::Initialize(const std::vector<std::string> &names) {
  for (unsigned i = 0; i < names.size(); i++) {
    dEdx.push_back( new TH1D(("dEdx_" + names[i]).c_str(), "dEdx", 200, 0., 20.));
  }
}

void ana::SBNOsc::TrajHistograms::Fill(const numu::RecoTrack &track, const anab::Calorimetry &collection_calo, const numu::RecoEvent &event, const std::vector<numu::TrackSelector> &selectors) {
  if (collection_calo.dEdx().size() == 0) return;
  for (unsigned i = 0; i < selectors.size(); i++) {
    if (selectors[i](track, event)) {
      // Ignore the first and last points as in the LArSoft ParticleID module
      for (unsigned j = 1; j < collection_calo.dEdx().size()-1; j++) {
        dEdx[i]->Fill(collection_calo.dEdx()[i]);
      } 
    }
  }
}

void ana::SBNOsc::TrajHistograms::Write() {
  for (unsigned i = 0; i < dEdx.size(); i++) {
    dEdx[i]->Write();
  }
}

void ana::SBNOsc::TrajHistograms::Get(TFile &f, const std::vector<std::string> &names) {
  for (unsigned i = 0; i < names.size(); i++) {
    dEdx.push_back( (TH1D *) f.Get(("dEdx_" + names[i]).c_str()) );
  }
}

void ana::SBNOsc::TrajHistograms::Add(const ana::SBNOsc::TrajHistograms &other) {
  for (unsigned i = 0; i < dEdx.size(); i++) {
    dEdx[i]->Add(other.dEdx[i]);
  }
}

ana::SBNOsc::TrajHistograms::~TrajHistograms() {}

