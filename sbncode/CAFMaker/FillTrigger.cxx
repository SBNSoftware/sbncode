#include<iostream>
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
	std::cout << "Gen triggerInfo.global_trigger_det_time: " << triggerInfo.global_trigger_det_time << std::endl;
	std::cout << "Gen triggerInfo.beam_gate_det_time: " << triggerInfo.beam_gate_det_time << std::endl;
    double diff_ts = triggerInfo.global_trigger_det_time - triggerInfo.beam_gate_det_time;
    std::cout << "diff_ts: " << diff_ts << std::endl;
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

  }

  void FillTriggerSBND(caf::SRSBNDTimingInfo& timingInfo, caf::SRTrigger& triggerInfo){
      
    triggerInfo.global_trigger_time = timingInfo.hltEtrig;
    triggerInfo.beam_gate_time_abs = timingInfo.hltBeamGate;
    std::cout << "triggerInfo.global_trigger_det_time: " << triggerInfo.global_trigger_det_time << std::endl;
	std::cout << "triggerInfo.beam_gate_det_time: " << triggerInfo.beam_gate_det_time << std::endl;
    double diff_ts = triggerInfo.global_trigger_det_time - triggerInfo.beam_gate_det_time;
	std::cout << "diff_ts: " << diff_ts << std::endl;
    triggerInfo.trigger_within_gate = diff_ts;
  }

  void FillTriggerEmulation(art::Handle<std::vector<int>> const& monpulsesFlat,
                             art::Handle<std::vector<int>> const& monpulseSizes,
                             art::Handle<int> const& numPairs,
                             art::Handle<bool> const& passedTrig,
                             caf::SRTrigger& triggerInfo) {

    triggerInfo.monpulses_flat = *monpulsesFlat;
    triggerInfo.monpulse_sizes = *monpulseSizes;
    triggerInfo.num_pairs_over_threshold = *numPairs;
    triggerInfo.passed_trigger = *passedTrig;
  }

  void FillSoftwareTrigger(const sbnd::trigger::pmtSoftwareTrigger& softInfo, caf::SRSoftwareTrigger& caf_softInfo){
    caf_softInfo.npmts = softInfo.nAboveThreshold;
    caf_softInfo.flash_peakpe = softInfo.peakPE;
    caf_softInfo.flash_peaktime = softInfo.peaktime + softInfo.trig_ts*1e-3; 
  }
}
