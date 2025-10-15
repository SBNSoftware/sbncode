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

  void FillTriggerMC(double absolute_time, caf::SRTrigger& triggerInfo) 
  {
    triggerInfo.global_trigger_time = absolute_time;
    triggerInfo.beam_gate_time_abs = absolute_time;
    triggerInfo.trigger_within_gate = 0.; // Set this to 0 since the "MC" trigger is (for now) always at the spill time
  }

  void FillTriggerSBND(caf::SRSBNDTimingInfo& timingInfo, 
                       caf::SRTrigger& triggerInfo)
  {
    triggerInfo.global_trigger_time = timingInfo.hltEtrig;
    triggerInfo.beam_gate_time_abs = timingInfo.hltBeamGate;
    double diff_ts = triggerInfo.global_trigger_det_time - triggerInfo.beam_gate_det_time;
    triggerInfo.trigger_within_gate = diff_ts;
  }

  void FillTriggerICARUS(const sbn::ExtraTriggerInfo& addltrig_info,
                         caf::SRTrigger& triggerInfo) 
  {
    // Per-cryostat additional trigger information, straight from the trigger hardware
    const auto& addltrig_info_cryoE = addltrig_info.cryostats[sbn::ExtraTriggerInfo::EastCryostat];
    const auto& addltrig_info_cryoW = addltrig_info.cryostats[sbn::ExtraTriggerInfo::WestCryostat];

    // Choose the cryostat that triggered first (if both are available)
    if (addltrig_info_cryoE.hasTrigger() && (!addltrig_info_cryoW.hasTrigger() || (addltrig_info_cryoE.beamToTrigger <= addltrig_info_cryoW.beamToTrigger))) {
      triggerInfo.trigger_cryo_source  = 0; ///< East
      triggerInfo.trigger_logic_bits   = addltrig_info_cryoE.triggerLogicBits;
      triggerInfo.beam_to_trigger_time = addltrig_info_cryoE.beamToTrigger;
    }
    else if (addltrig_info_cryoW.hasTrigger()) {
      triggerInfo.trigger_cryo_source  = 1; ///< West
      triggerInfo.trigger_logic_bits   = addltrig_info_cryoW.triggerLogicBits;
      triggerInfo.beam_to_trigger_time = addltrig_info_cryoW.beamToTrigger;
    }
  }

  void FillTriggerEmulation(art::Handle<std::vector<int>> const& monpulsesFlat,
                            art::Handle<std::vector<int>> const& monpulseSizes,
                            art::Handle<int> const& numPairs,
                            art::Handle<bool> const& passedTrig,
                            caf::SRTrigger& triggerInfo) 
  {
    triggerInfo.monpulses_flat = *monpulsesFlat;
    triggerInfo.monpulse_sizes = *monpulseSizes;
    triggerInfo.num_pairs_over_threshold = *numPairs;
    triggerInfo.passed_trigger = *passedTrig;
  }

}
