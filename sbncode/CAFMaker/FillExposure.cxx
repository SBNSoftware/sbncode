
#include "sbncode/CAFMaker/FillExposure.h"

namespace caf
{

  caf::SRBNBInfo makeSRBNBInfo(sbn::BNBSpillInfo const& info)
  {
    caf::SRBNBInfo single_store;
    single_store.spill_time_sec = info.spill_time_s;
    single_store.spill_time_nsec = info.spill_time_ns;
    single_store.event = info.event;
    single_store.TOR860 = info.TOR860;
    single_store.TOR875 = info.TOR875;
    single_store.LM875A = info.LM875A;
    single_store.LM875B = info.LM875B;
    single_store.LM875C = info.LM875C;
    single_store.HP875 = info.HP875;
    single_store.VP875 = info.VP875;
    single_store.HPTG1 = info.HPTG1;
    single_store.VPTG1 = info.VPTG1;
    single_store.HPTG2 = info.HPTG2;
    single_store.VPTG2 = info.VPTG2;
    single_store.BTJT2 = info.BTJT2;
    single_store.THCURR = info.THCURR;
    single_store.M875BB = info.M875BB;
    single_store.M876BB = info.M876BB;
    single_store.MMBTBB = info.MMBTBB;
    single_store.M875BB_spill_time_diff = info.M875BB_spill_time_diff;
    single_store.M876BB_spill_time_diff = info.M876BB_spill_time_diff;
    single_store.MMBTBB_spill_time_diff = info.MMBTBB_spill_time_diff;
    return single_store;
  }

  void FillExposure(const std::vector<sbn::BNBSpillInfo>& bnb_spill_info,
		    std::vector<caf::SRBNBInfo>& BNBInfo,
		    double& subRunPOT)
  {
    for(const sbn::BNBSpillInfo& info: bnb_spill_info)
      {
	subRunPOT += info.POT();
        BNBInfo.push_back(makeSRBNBInfo(info));
      }
  }
  void FillExposureNuMI(const std::vector<sbn::NuMISpillInfo>& numi_spill_info,
		    std::vector<caf::SRNuMIInfo>& NuMIInfo,
		    double& subRunPOT) {
    for (const sbn::NuMISpillInfo &info: numi_spill_info) {
      subRunPOT += info.POT();

      NuMIInfo.emplace_back();
      NuMIInfo.back().HP121 = info.HP121;
      NuMIInfo.back().VP121 = info.VP121;
      NuMIInfo.back().HPTGT = info.HPTGT;
      NuMIInfo.back().VPTGT = info.VPTGT;
      NuMIInfo.back().HITGT = info.HITGT;
      NuMIInfo.back().VITGT = info.VITGT;
      NuMIInfo.back().MTGTDS = info.MTGTDS;
      NuMIInfo.back().HRNDIR = info.HRNDIR;
      NuMIInfo.back().NSLINA = info.NSLINA;
      NuMIInfo.back().NSLINB = info.NSLINB;
      NuMIInfo.back().NSLINC = info.NSLINC;
      NuMIInfo.back().NSLIND = info.NSLIND;
      NuMIInfo.back().TRTGTD = info.TRTGTD;
      NuMIInfo.back().TR101D = info.TR101D;
      NuMIInfo.back().TORTGT = info.TORTGT;
      NuMIInfo.back().TOR101 = info.TOR101;
      NuMIInfo.back().time = info.time;
      NuMIInfo.back().spill_time_s = info.spill_time_s;
      NuMIInfo.back().spill_time_ns = info.spill_time_ns;
      NuMIInfo.back().event = info.event;
      NuMIInfo.back().daq_gates = info.daq_gates;
    }
  }

  caf::SRNuMIInfo makeSRNuMIInfo(sbn::NuMISpillInfo const& info)
  {
    caf::SRNuMIInfo single_store;

    single_store.HP121 = info.HP121;
    single_store.VP121 = info.VP121;
    single_store.HPTGT = info.HPTGT;
    single_store.VPTGT = info.VPTGT;
    single_store.HITGT = info.HITGT;
    single_store.VITGT = info.VITGT;
    single_store.MTGTDS = info.MTGTDS;
    single_store.HRNDIR = info.HRNDIR;
    single_store.NSLINA = info.NSLINA;
    single_store.NSLINB = info.NSLINB;
    single_store.NSLINC = info.NSLINC;
    single_store.NSLIND = info.NSLIND;
    single_store.TRTGTD = info.TRTGTD;
    single_store.TR101D = info.TR101D;
    single_store.TORTGT = info.TORTGT;
    single_store.TOR101 = info.TOR101;
    single_store.time = info.time;
    single_store.spill_time_s = info.spill_time_s;
    single_store.spill_time_ns = info.spill_time_ns;
    single_store.event = info.event;
    single_store.daq_gates = info.daq_gates;

    return single_store;
  }

}
