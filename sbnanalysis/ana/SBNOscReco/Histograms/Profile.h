#ifndef _sbnnumu_Profile_h__
#define _sbnnumu_Profile_h__

#include "DynamicSelector.h"
#include "HistoList.h"
#include "../Data/RecoEvent.h"

class TH2D;

namespace ana {
 namespace SBNOsc {

struct TrackProfiles : public HistoList {
  TH2D *range_bias;
  TH2D *range_var;
  TH2D *mcs_bias;
  TH2D *mcs_var;

  std::vector<TH1 *> all_histos;

  void Initialize(const std::string &postfix, unsigned nbinsx, double xlo, double xhi);
  void Fill(const numu::ROOTValue &val, const numu::RecoTrack &track, const numu::RecoEvent &event);
};
  }
}

#endif
