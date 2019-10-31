#ifndef __sbnanalysis_TrajHistogram_HH
#define __sbnanalysis_TrajHistogram_HH

#include <string>
#include <vector>

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "../Data/RecoTrack.h"

#include "DynamicSelector.h"

class TH1D;
class TH2D;
class TFile;
namespace ana {
 namespace SBNOsc {

/**
 * Histograms that require filling over the full trajectory of a reco track
 */
struct TrajHistograms {
  std::vector<TH1D *> dEdx;

  void Initialize(const std::vector<std::string> &names);
  void Fill(const numu::RecoTrack &track, const anab::Calorimetry &collection_calo, const numu::RecoEvent &event, const std::vector<numu::TrackSelector> &selectors);
  void Write();
  void Get(TFile &f, const std::vector<std::string> &names);
  void Add(const ana::SBNOsc::TrajHistograms &other);
  ~TrajHistograms();
};


  } // namespace SBNOSc
} // namespace ana

#endif
