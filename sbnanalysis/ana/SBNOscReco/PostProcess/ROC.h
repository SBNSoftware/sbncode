#ifndef _sbncode_ROC_hh
#define _sbncode_ROC_hh

#include <vector>

#include "TH1D.h"

#include "../Data/RecoEvent.h"

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
    void Fill(const numu::RecoEvent &event);
    /** Write the ROC plots to disk */
    void Write() const;
    /** Get the optimal cut values by S/sqrt(S+B) 
   *  \param scale_signal Amount to scale signal interactions by
   *  \param scale_background Amount to scale background interactions by
   *  \param n_background_data Number of background data events to calculate statistical uncertainties
   */
    void BestCuts(float scale_signal, float scale_background, float n_background_data) const;

    /**
   * Primitive used to make ROC plots 
   */
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
