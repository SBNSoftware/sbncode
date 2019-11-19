#ifndef _sbnanalysis_numu_TriggerEmulator_hh
#define _sbnanalysis_numu_TriggerEmulator_hh

#include <vector>

#include "lardataobj/RawData/OpDetWaveform.h"

#include "../Data/DetInfo.h"

namespace numu {

std::vector<FlashTriggerPrimitive> TriggerPrimitives(
  const std::vector<raw::OpDetWaveform> &waveforms, 
  double tick_period, 
  std::pair<double, double> &window, 
  int thresh,
  bool is_sbnd);

bool HasTrigger(const std::vector<FlashTriggerPrimitive> &primitives, int threshold, unsigned n_above_threshold);
}

#endif
