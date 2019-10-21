#ifndef __sbnanalysis_TrajHistogram_HH
#define __sbnanalysis_TrajHistogram_HH

#include <string>
#include <vector>

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "../Data/RecoTrack.h"

#include "Histograms.h"

class TH1D;
class TH2D;
namespace ana {
 namespace SBNOsc {

/**
 * Histograms that require filling over the full trajectory of a reco track
 */
struct TrajHistograms {
  TH1D *dEdx[TrackHistos::nPDGs];

  void Initialize();
  void Fill(const numu::RecoTrack &track, const anab::Calorimetry &collection_calo);
  void Write();
  void Get(TFile &f);
  void Add(const ana::SBNOsc::TrajHistograms &other);
  ~TrajHistograms();
};


  } // namespace SBNOSc
} // namespace ana

#endif
