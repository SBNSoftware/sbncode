#include "ROC.h"

#include "TGraph.h"

void ana::SBNOsc::ROC::Initialize() {
  crt_track_angle.Initialize("crt_track_angle", 0, 6, 360);
  fAllPrimitives.push_back(&crt_track_angle);
  crt_hit_distance.Initialize("crt_hit_distance", 0, 50, 500);
  fAllPrimitives.push_back(&crt_hit_distance);
}

void ana::SBNOsc::ROC::Write() const {
  for (const Primitive *prim: fAllPrimitives) {
    prim->Write();
  }
}

void ana::SBNOsc::ROC::BestCuts(float scale_signal, float scale_background, float n_background_data) const {
  for (const Primitive *prim: fAllPrimitives) {
    std::cout << "Cut (" << prim->name << ") value: " << prim->BestCut(scale_signal, scale_background, n_background_data);
  }
}

void ana::SBNOsc::ROC::Fill(const numu::RecoEvent &event) {
  for (const numu::RecoInteraction &reco: event.reco) {
    bool is_signal = reco.match.has_match && (reco.match.mode == numu::mCC || reco.match.mode == numu::mNC);
    const numu::RecoTrack &track = event.reco_tracks.at(reco.slice.primary_track_index);
    
    // CRT stuff
    if (track.crt_match.size() > 0) {
      if (track.crt_match[0].track.present) {
        crt_track_angle.Fill(is_signal, track.crt_match[0].track.angle);
      }
      else {
        crt_track_angle.FillNever(is_signal);
      }
      if (track.crt_match[0].hit.present) {
        crt_hit_distance.Fill(is_signal, track.crt_match[0].hit.distance);
      }
      else {
        crt_hit_distance.FillNever(is_signal);
      }
    }
    else {
      crt_track_angle.FillNever(is_signal);
      crt_hit_distance.FillNever(is_signal);
    }
  }
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
  while (bin <= hist->GetNbinsX() && value < hist->GetXaxis()->GetBinCenter(bin)) {
    hist->Fill(hist->GetXaxis()->GetBinCenter(bin));
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


float ana::SBNOsc::ROC::Primitive::BestCut(float scale_signal, float scale_background, float n_background_data) const {
  float max_significance = 0.;
  float best_cut = 0.;
  for (unsigned bin = 1; bin <= signal->GetNbinsX(); bin++) {
    float sig = signal->GetBinContent(bin) * scale_signal;
    float bkg = background->GetBinContent(bin) * (scale_background * n_background / n_background_data);
    float this_significance = sig / sqrt(sig + bkg);
    if (this_significance > max_significance) {
      max_significance = this_significance;
      best_cut = signal->GetBinCenter(bin); 
    }
  } 
  return best_cut;
}

void ana::SBNOsc::ROC::Primitive::Write() const {
  std::vector<float> eff;
  std::vector<float> purity;
  for (unsigned bin = 1; bin <= signal->GetNbinsX(); bin++) {
    float this_eff = signal->GetBinContent(bin) / n_signal;
    float this_purity = signal->GetBinContent(bin) / (background->GetBinContent(bin) + signal->GetBinContent(bin));
    eff.push_back(this_eff);
    purity.push_back(this_purity);
  }
  TGraph *ROC = new TGraph(signal->GetNbinsX(), &purity[0], &eff[0]);
  ROC->SetTitle(name.c_str());
  ROC->Write();

  delete ROC;
}

ana::SBNOsc::ROC::Primitive::~Primitive() {
  delete signal;
  delete background;
}


