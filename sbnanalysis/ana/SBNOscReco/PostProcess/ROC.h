#ifndef _sbncode_ROC_hh
#define _sbncode_ROC_hh

#include <vector>

#include "TH1D.h"

#include "../Data/RecoEvent.h"

#include "Cuts.h"

namespace ana {
  namespace SBNOsc {

  /**
 * Class to handle making of ROC plots for cut values
 */ 
  struct ROC {
    /** Initialize the class */
    void Initialize();
    /** Fill the class with a reconstructed event 
   * \param event The RecoEvent information
   */
    void Fill(const Cuts cuts, const numu::RecoEvent &event);
    /** Write the ROC plots to disk */
    void Write() const;
    /** Get the optimal cut values by S/sqrt(S+B) 
   */
    void BestCuts() const;
    void Normalize(float neutrino_scale, float cosmic_scale);

    /**
   * Primitive used to make ROC plots 
   */
    struct Primitive {
      void Initialize(const std::string &name, float cut_low, float cut_high, unsigned n_bin);
      void Fill(bool is_signal, float value);
      void FillNever(bool is_signal);
      void FillAlways(bool is_signal);
      void Scale(float scale);
      ~Primitive();

      std::string name;
      TH1D *signal;
      unsigned n_signal;
      TH1D *background;
      unsigned n_background;
    };

    struct NormalizedPrimitive {
      Primitive fNeutrino;
      Primitive fCosmic;
      std::string name;

      void Initialize(const std::string &name, float cut_low, float cut_high, unsigned n_bin);
      void FillNeutrino(bool is_signal, float value) {fNeutrino.Fill(is_signal, value);}
      void FillNeverNeutrino(bool is_signal) {fNeutrino.FillNever(is_signal);}
      void FillAlwaysNeutrino(bool is_signal) {fNeutrino.FillAlways(is_signal);}

      void FillCosmic(bool is_signal, float value) {fCosmic.Fill(is_signal, value);}
      void FillNeverCosmic(bool is_signal) {fCosmic.FillNever(is_signal);}
      void FillAlwaysCosmic(bool is_signal) {fCosmic.FillAlways(is_signal);}

      float Normalize(float scale_neutrino, float scale_cosmic);
      void Write() const;
      // float BestCut(float scale_neutrino, float scale_cosmic, float n_background_data) const;
      float BestCut() const;

      float Signal(unsigned bin) const;
      float Background(unsigned bin) const;
      float NSignal() const;
      float NBackground() const;
      unsigned NCutVals() const;
      float GetCutVal(unsigned bin) const;
    };

    NormalizedPrimitive crt_track_angle;
    NormalizedPrimitive crt_hit_distance;
    std::vector<NormalizedPrimitive *> fAllPrimitives;
  };

  } // namespace SBNOsc
} // namespace ana



#endif
