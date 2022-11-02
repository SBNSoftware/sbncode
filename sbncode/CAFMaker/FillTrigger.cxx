
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
    double diff_ts = (trigger_ts - beam_ts) / 1000.;
    triggerInfo.trigger_within_gate = diff_ts;
    triggerInfo.beam_gate_det_time = trig.BeamGateTime();
    triggerInfo.global_trigger_det_time = trig.TriggerTime();
  }

  void FillTriggerMC(double absolute_time, caf::SRTrigger& triggerInfo) {
    triggerInfo.global_trigger_time = absolute_time;
    triggerInfo.beam_gate_time_abs = absolute_time;

    // Set this to 0 since the "MC" trigger is (for now) always at the spill time
    triggerInfo.trigger_within_gate = 0.; 

    // TODO: fill others?
  }
    
}
