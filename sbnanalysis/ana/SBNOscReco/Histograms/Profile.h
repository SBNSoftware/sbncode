#ifndef _sbnnumu_Profile_h__
#define _sbnnumu_Profile_h__

#include "DynamicSelector.h"
#include "HistoList.h"
#include "../Data/RecoEvent.h"
#include "../MultiThread/THShared.h"

namespace ana {
 namespace SBNOsc {

struct TrackProfiles : public HistoList {
  TH3Shared range_minus_true;
  TH3Shared range_v_true_mom;
  TH3Shared mcs_minus_true;
  TH3Shared mcs_v_true_mom;
  TH3Shared pid_confusion_tr;

  void Initialize(const std::string &postfix, unsigned nbinsx, double xlo, double xhi);
  void Fill(const numu::ROOTValue &val, const numu::RecoTrack &track, const numu::RecoEvent &event);
};
  }
}

#endif
