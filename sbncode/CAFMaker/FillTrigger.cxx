
#include "sbncode/CAFMaker/FillTrigger.h"
#include "sbnobj/Common/Trigger/BeamBits.h"

namespace caf
{
  void FillTrigger(const sbn::ExtraTriggerInfo& addltrig_info,
		   const raw::Trigger& trig,
		   caf::SRTrigger& triggerInfo,
		   caf::TimeRefShifter<> const& shifter /* = {} */
		   )
  {
    triggerInfo.global_trigger_time = addltrig_info.triggerTimestamp;
    triggerInfo.beam_gate_time_abs = addltrig_info.beamGateTimestamp;
    triggerInfo.beam_gate_det_time = shifter.shiftedTime(trig.BeamGateTime());
    triggerInfo.global_trigger_det_time = shifter.shiftedTime(trig.TriggerTime());
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

  void FillTriggerSBND(caf::SRSBNDTimingInfo& timingInfo, caf::SRTrigger& triggerInfo){
    // updates the existing triggerInfo record
    triggerInfo.global_trigger_time = timingInfo.hltEtrig;
    triggerInfo.beam_gate_time_abs = timingInfo.hltBeamGate;
  }

}
