#include "ROC.h"

#include "TGraph.h"

void ana::SBNOsc::ROC::Initialize() {
  crt_track_angle.Initialize("crt_track_angle", 0, 6, 60);
  fAllPrimitives.push_back(&crt_track_angle);
  crt_hit_distance.Initialize("crt_hit_distance", 0, 100, 100);
  fAllPrimitives.push_back(&crt_hit_distance);
}

void ana::SBNOsc::ROC::Write() const {
  for (const NormalizedPrimitive *prim: fAllPrimitives) {
    prim->Write();
  }
}

void ana::SBNOsc::ROC::Normalize(float neutrino_scale, float cosmic_scale) {
  for (NormalizedPrimitive *prim: fAllPrimitives) {
    prim->Normalize(neutrino_scale, cosmic_scale);
  }
}

void ana::SBNOsc::ROC::BestCuts() const {
  for (const NormalizedPrimitive *prim: fAllPrimitives) {
    std::cout << "Cut (" << prim->name << ") value: " << prim->BestCut() << std::endl;
  }
}

void ana::SBNOsc::ROC::Fill(const ana::SBNOsc::Cuts cuts, const numu::RecoEvent &event, bool file_is_neutrino) {
  for (const numu::RecoInteraction &reco: event.reco) {
    bool is_signal = reco.match.mode == numu::mCC;
    bool is_bkg = reco.match.mode == numu::mCosmic || reco.match.mode == numu::mIntimeCosmic;
    if (!is_signal && !is_bkg) continue; // ignore NC stuff for now
    
    const numu::RecoTrack &track = event.reco_tracks.at(reco.slice.primary_track_index);
    
    // CRT stuff
    if (track.crt_match.track.present) {
      if (file_is_neutrino) crt_track_angle.FillNeutrino(is_signal, track.crt_match.track.angle);
      else crt_track_angle.FillCosmic(is_signal, track.crt_match.track.angle);
    }
    else {
      if (file_is_neutrino) crt_track_angle.FillNeverNeutrino(is_signal);
      else crt_track_angle.FillNeverCosmic(is_signal);
    }
    if (track.crt_match.hit_match.present && !cuts.TimeInSpill(track.crt_match.hit_match.time)) {
      if (file_is_neutrino) crt_hit_distance.FillNeutrino(is_signal, track.crt_match.hit_match.distance);
      else crt_hit_distance.FillCosmic(is_signal, track.crt_match.hit_match.distance);
    }
    else {
      if (file_is_neutrino) crt_hit_distance.FillNeverNeutrino(is_signal);
      else crt_hit_distance.FillNeverCosmic(is_signal);
    }
  }
}

void ana::SBNOsc::ROC::NormalizedPrimitive::Initialize(const std::string &this_name, float cut_low, float cut_high, unsigned n_bin) {
  fNeutrino.Initialize(this_name + "neutrino_", cut_low, cut_high, n_bin);
  fCosmic.Initialize(this_name + "cosmic_", cut_low, cut_high, n_bin);
  name = this_name;
}

void ana::SBNOsc::ROC::Primitive::Initialize(const std::string &this_name, float cut_low, float cut_high, unsigned n_bin) {
  signal = new TH1D((this_name + "signal").c_str(), this_name.c_str(), n_bin, cut_low, cut_high);
  background = new TH1D((this_name + "background").c_str(), this_name.c_str(), n_bin, cut_low, cut_high);
  n_signal = 0;
  n_background = 0;
  name = this_name;
}

void ana::SBNOsc::ROC::Primitive::Fill(bool is_signal, float value) {
  TH1D *hist = is_signal ? signal : background;
  unsigned bin = 1;
  while (bin <= hist->GetNbinsX() && value > hist->GetXaxis()->GetBinCenter(bin)) {
    hist->Fill(hist->GetXaxis()->GetBinCenter(bin));
    bin += 1;
  }
  n_signal += is_signal;
  n_background += !is_signal;
}

void ana::SBNOsc::ROC::Primitive::FillNever(bool is_signal) {
  TH1D *hist = is_signal ? signal : background;
  unsigned bin = 1;
  for (unsigned bin = 1; bin <= hist->GetNbinsX(); bin++) {
    hist->Fill(hist->GetXaxis()->GetBinCenter(bin));
  }
  n_signal += is_signal;
  n_background += !is_signal;
}

void ana::SBNOsc::ROC::Primitive::FillAlways(bool is_signal) {
  n_signal += is_signal;
  n_background += !is_signal;
}

void ana::SBNOsc::ROC::Primitive::Scale(float scale) {
  signal->Scale(scale);
  background->Scale(scale);
  n_signal *= scale;
  n_background *= scale;
}

float ana::SBNOsc::ROC::NormalizedPrimitive::Normalize(float scale_neutrino, float scale_cosmic) {
  fNeutrino.Scale(scale_neutrino);
  fCosmic.Scale(scale_cosmic);
}

float ana::SBNOsc::ROC::NormalizedPrimitive::Signal(unsigned bin) const {
  return fNeutrino.signal->GetBinContent(bin) + fCosmic.signal->GetBinContent(bin);
}

float ana::SBNOsc::ROC::NormalizedPrimitive::Background(unsigned bin) const {
  return fNeutrino.background->GetBinContent(bin) + fCosmic.background->GetBinContent(bin);
}

unsigned ana::SBNOsc::ROC::NormalizedPrimitive::NCutVals() const {
  return fNeutrino.signal->GetNbinsX();
}

float ana::SBNOsc::ROC::NormalizedPrimitive::GetCutVal(unsigned bin) const {
  return fNeutrino.signal->GetBinCenter(bin);
}

float ana::SBNOsc::ROC::NormalizedPrimitive::NSignal() const {
  return fNeutrino.n_signal + fCosmic.n_signal;
}

float ana::SBNOsc::ROC::NormalizedPrimitive::NBackground() const {
  return fNeutrino.n_background + fCosmic.n_background;
}

float ana::SBNOsc::ROC::NormalizedPrimitive::BestCut() const {
  float max_significance = 0.;
  float best_cut = 0.;
  for (unsigned bin = 1; bin <= NCutVals(); bin++) {
    float sig = Signal(bin);
    float bkg = Background(bin);
    float this_significance = 0.;
    if (sig + bkg > 1.e-4) {
      this_significance = sig / sqrt(sig + bkg);
    }
    if (this_significance > max_significance) {
      max_significance = this_significance;
      best_cut = GetCutVal(bin);
    }
  } 
  return best_cut;
}

void ana::SBNOsc::ROC::NormalizedPrimitive::Write() const {
  std::vector<float> eff;
  std::vector<float> rejection;
  for (unsigned bin = 1; bin <= NCutVals(); bin++) {
    float this_eff = Signal(bin) / NSignal();
    float this_rej = 1 - Background(bin) / NBackground();
    eff.push_back(this_eff);
    rejection.push_back(this_rej);
  }
  TGraph *ROC = new TGraph(NCutVals(), &eff[0], &rejection[0]);
  ROC->SetTitle(name.c_str());
  ROC->SetName(name.c_str());
  ROC->Write();

  delete ROC;
}

ana::SBNOsc::ROC::Primitive::~Primitive() {
  delete signal;
  delete background;
}


