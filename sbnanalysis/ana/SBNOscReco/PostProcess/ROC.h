#ifndef _sbncode_ROC_hh
#define _sbncode_ROC_hh

#include <vector>

#include "TH1D.h"

#include "../Data/RecoEvent.h"

namespace ana {
  namespace SBNOsc {

  struct ROC {
    void Initialize();
    void Fill(const numu::RecoEvent &event);
    void Write() const;
    void BestCuts(float scale_signal, float scale_background, float n_background_data) const;

    struct Primitive {
      void Initialize(const std::string &name, float cut_low, float cut_high, unsigned n_bin);
      void Fill(bool is_signal, float value);
      void FillNever(bool is_signal);
      void Write() const;
      float BestCut(float scale_signal, float scale_background, float n_background_data) const;
      ~Primitive();

      std::string name;
      TH1D *signal;
      unsigned n_signal;
      TH1D *background;
      unsigned n_background;
    };

    Primitive crt_track_angle;
    Primitive crt_hit_distance;
    std::vector<Primitive *> fAllPrimitives;
  };

  } // namespace SBNOsc
} // namespace ana



#endif
