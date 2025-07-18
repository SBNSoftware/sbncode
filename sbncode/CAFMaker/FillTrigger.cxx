
#include "sbncode/CAFMaker/FillTrigger.h"

namespace caf
{
  void FillTrigger(const sbn::ExtraTriggerInfo& addltrig_info,
		   const raw::Trigger& trig,
		   caf::SRTrigger& triggerInfo,
		   const double time_offset = 0.0)
  {
    triggerInfo.global_trigger_time = addltrig_info.triggerTimestamp;
    triggerInfo.beam_gate_time_abs = addltrig_info.beamGateTimestamp;
    triggerInfo.beam_gate_det_time = trig.BeamGateTime() + time_offset;
    triggerInfo.global_trigger_det_time = trig.TriggerTime() + time_offset;
    double diff_ts = triggerInfo.global_trigger_det_time - triggerInfo.beam_gate_det_time;
    triggerInfo.trigger_within_gate = diff_ts;

    triggerInfo.prev_global_trigger_time = addltrig_info.previousTriggerTimestamp;
    triggerInfo.source_type = sbn::bits::value(addltrig_info.sourceType);
    triggerInfo.trigger_type = sbn::bits::value(addltrig_info.triggerType);
    triggerInfo.trigger_id = addltrig_info.triggerID;
    triggerInfo.gate_id = addltrig_info.gateID;
    triggerInfo.trigger_count = addltrig_info.triggerCount;
    triggerInfo.gate_count = addltrig_info.gateCount;
    triggerInfo.gate_delta = addltrig_info.gateCountFromPreviousTrigger;
  }

  void FillTriggerMC(double absolute_time, caf::SRTrigger& triggerInfo) {
    triggerInfo.global_trigger_time = absolute_time;
    triggerInfo.beam_gate_time_abs = absolute_time;

    // Set this to 0 since the "MC" trigger is (for now) always at the spill time
    triggerInfo.trigger_within_gate = 0.;

    // TODO: fill others?
  }

}
