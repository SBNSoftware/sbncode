
#include "FillExposure.h"

namespace caf
{
  void FillExposure(const std::vector<art::Ptr<sbn::BNBSpillInfo> > bnb_spill_info,
		    std::vector<caf::SRBNBInfo>& BNBInfo,
		    double& subRunPOT)
  {
    for(const art::Ptr<sbn::BNBSpillInfo>& info: bnb_spill_info)
      {
	caf::SRBNBInfo single_store;
        subRunPOT += info->TOR875;
        single_store.spill_time_sec = info->spill_time_s;
        single_store.spill_time_nsec = info->spill_time_ns;
        single_store.event = info->event;
        single_store.TOR860 = info->TOR860;
        single_store.TOR875 = info->TOR875;
        single_store.LM875A = info->LM875A;
        single_store.LM875B = info->LM875B;
        single_store.LM875C = info->LM875C;
        single_store.HP875 = info->HP875;
        single_store.VP875 = info->VP875;
        single_store.HPTG1 = info->HPTG1;
        single_store.VPTG1 = info->VPTG1;
        single_store.HPTG2 = info->HPTG2;
        single_store.VPTG2 = info->VPTG2;
        single_store.BTJT2 = info->BTJT2;
        single_store.THCURR = info->THCURR;
        single_store.M875BB = info->M875BB;
        single_store.M876BB = info->M876BB;
        single_store.MMBTBB = info->MMBTBB;
        single_store.M875BB_spill_time_diff = info->M875BB_spill_time_diff;
	single_store.M876BB_spill_time_diff = info->M876BB_spill_time_diff;
        single_store.MMBTBB_spill_time_diff = info->MMBTBB_spill_time_diff;
        BNBInfo.push_back(single_store);
      }
  }
}
