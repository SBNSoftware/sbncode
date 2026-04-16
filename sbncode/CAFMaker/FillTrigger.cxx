#include<iostream>
#include <bitset>
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

  void FillTriggerICARUS(const sbn::ExtraTriggerInfo& addltrig_info,
                         caf::SRTrigger& triggerInfo) 
  {
    // Choose the cryostat that triggered first (if both are available)
    int const cryo = addltrig_info.cryostats[sbn::ExtraTriggerInfo::EastCryostat].beamToTrigger < addltrig_info.cryostats[sbn::ExtraTriggerInfo::WestCryostat].beamToTrigger
      ? sbn::ExtraTriggerInfo::EastCryostat 
      : sbn::ExtraTriggerInfo::WestCryostat;

    sbn::ExtraTriggerInfo::CryostatInfo const& cryoInfo = addltrig_info.cryostats[cryo];
    if (cryoInfo.hasTrigger()) {
      triggerInfo.trigger_cryo_source  = cryo;
      triggerInfo.trigger_logic_bits   = cryoInfo.triggerLogicBits;
      triggerInfo.gate_to_trigger_time = static_cast<float>(cryoInfo.beamToTrigger) / 1000.0f; // [us]
    }  
  }

  void FillTriggerSBND(caf::SRSBNDTimingInfo& timingInfo, caf::SRTrigger& triggerInfo){
      
    triggerInfo.global_trigger_time = timingInfo.hltEtrig;
    triggerInfo.beam_gate_time_abs = timingInfo.hltBeamGate;
    double diff_ts = triggerInfo.global_trigger_det_time - triggerInfo.beam_gate_det_time;
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

  void FillSoftwareTriggerSBND(const sbnd::trigger::pmtSoftwareTrigger& softInfo, caf::SRSoftwareTrigger& caf_softInfo){
    caf_softInfo.npmts = softInfo.nAboveThreshold;
    caf_softInfo.flash_peakpe = softInfo.peakPE;
    caf_softInfo.flash_peaktime = softInfo.peaktime + softInfo.trig_ts*1e-3; 
  }

  void FillPTBTriggersSBND(const std::vector<sbn::pot::PTBInfo_t>& ptb_triggers, caf::SRTrigger& triggerInfo) {
    triggerInfo.ptb_hlt_timestamp.clear();
    triggerInfo.ptb_hlt_bit.clear();
    triggerInfo.ptb_llt_timestamp.clear();
    triggerInfo.ptb_llt_bit.clear();
    
    // Decode trigger words: each set bit becomes a separate entry with the same timestamp
    for(const auto& trig : ptb_triggers) {
      // Choose destination vectors based on trigger type
      auto& ptb_timestamp = trig.isHLT ? triggerInfo.ptb_hlt_timestamp : triggerInfo.ptb_llt_timestamp;
      auto& ptb_bit = trig.isHLT ? triggerInfo.ptb_hlt_bit : triggerInfo.ptb_llt_bit;
      
      std::bitset<64> const triggerWord { trig.triggerWord };
      // currPTBTimeStamp is already in nanoseconds (uint64_t), use directly
      // Loop variable is uint8_t since we know bit values are 0-63 (fits in uint8_t)
      for(std::uint8_t bit = 0; bit < triggerWord.size(); ++bit) {
        if(!triggerWord[bit]) continue;
        ptb_timestamp.push_back(trig.currPTBTimeStamp);
        ptb_bit.push_back(bit);
      }
    }
  }
}
