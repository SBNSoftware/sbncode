
#include "sbncode/CAFMaker/FillTrigger.h"

#include <iostream>

namespace caf
{
  void FillTrigger(const sbn::ExtraTriggerInfo& addltrig_info,
		   const std::vector<raw::Trigger>& trig_info,
		   std::vector<caf::SRTrigger>& triggerInfo)
  {
    uint64_t beam_ts = 0;
    uint64_t trigger_ts = 0;
    triggerInfo.emplace_back();
    triggerInfo.back().global_trigger_time = addltrig_info.triggerTimestamp;
    triggerInfo.back().beam_gate_time_abs = addltrig_info.beamGateTimestamp;
    beam_ts = addltrig_info.beamGateTimestamp;
    trigger_ts = addltrig_info.triggerTimestamp;
    std::cout << addltrig_info.triggerTimestamp << " " << addltrig_info.beamGateTimestamp << std::endl;
    int64_t diff_ts = trigger_ts - beam_ts;
    triggerInfo.back().trigger_within_gate = diff_ts;
    for(const raw::Trigger& trig: trig_info)
    {
      triggerInfo.back().beam_gate_det_time = trig.BeamGateTime();
      triggerInfo.back().global_trigger_det_time = trig.TriggerTime();
    }
  }
    
}
