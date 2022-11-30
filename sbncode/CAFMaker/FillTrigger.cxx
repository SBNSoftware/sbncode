
#include "sbncode/CAFMaker/FillTrigger.h"

namespace caf
{
  void FillTrigger(const sbn::ExtraTriggerInfo& addltrig_info,
		   const raw::Trigger& trig,
		   caf::SRTrigger& triggerInfo)
  {
    uint64_t beam_ts = 0;
    uint64_t trigger_ts = 0;
    triggerInfo.global_trigger_time = addltrig_info.triggerTimestamp;
    triggerInfo.beam_gate_time_abs = addltrig_info.beamGateTimestamp;
    beam_ts = addltrig_info.beamGateTimestamp;
    trigger_ts = addltrig_info.triggerTimestamp;
    triggerInfo.beam_gate_det_time = trig.BeamGateTime();
    triggerInfo.global_trigger_det_time = trig.TriggerTime();
    double diff_ts = triggerInfo.global_trigger_det_time - triggerInfo.beam_gate_det_time;
    triggerInfo.trigger_within_gate = diff_ts;
  }
    
}
