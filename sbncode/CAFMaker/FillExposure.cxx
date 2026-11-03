
#include "sbncode/CAFMaker/FillExposure.h"

#include <cmath>
#include <limits>

namespace caf
{

  caf::SRBNBInfo makeSRBNBInfo(sbn::BNBSpillInfo const& info)
  {
    // Legacy FOM field: the FOM with the fitted multiwire width, else with the
    // database width, else with the nominal width. Valid values are in [0, 1] (a
    // FOM of exactly 0, beam fully off target, is valid); getBNBqualityFOM()
    // flags failures with -999, -1, 2, 3.
    auto const inRange = [](float v){ return (v >= 0.0f) && (v <= 1.0f); };
    // input written by a release without the beam-quality result (fom_version 0)
    bool const olderInput = (info.fom_version == 0);
    // back-up FOMs present: written by a release that had them (SBN2025A, or
    // later); otherwise they read back as their default -999
    bool const hasBackupFOMs = (info.PreFitFOM != -999.0f) || (info.NoMultiWireFOM != -999.0f);
    float finalFOM = -999.;
    if (olderInput && !hasBackupFOMs) {
      // the stored value, as the production CAFs of that input have it
      finalFOM = info.FOM;
    }
    else if (inRange(info.FOM)) finalFOM = info.FOM;
    else if (inRange(info.PreFitFOM)) finalFOM = info.PreFitFOM;
    else if (inRange(info.NoMultiWireFOM)) {
      // nominal width: used as is (fom_source 5, fom_status bit 64). The releases
      // without the beam-quality result stored it as 100 + FOM, so that the
      // standard 0.98 < FOM <= 1 cut rejects it: kept for their input, so that a
      // re-made CAF has the same FOM field as their production CAFs.
      finalFOM = olderInput? 100.0f + info.NoMultiWireFOM: info.NoMultiWireFOM;
    }

    caf::SRBNBInfo single_store;
    single_store.spill_time_sec = info.spill_time_s;
    single_store.spill_time_nsec = info.spill_time_ns;
    single_store.event = info.event;
    single_store.TOR860 = info.TOR860;
    single_store.TOR875 = info.TOR875;
    single_store.LM875A = info.LM875A;
    single_store.LM875B = info.LM875B;
    single_store.LM875C = info.LM875C;
    single_store.HP873 = info.HP873;
    single_store.VP873 = info.VP873;
    single_store.HP875 = info.HP875;
    single_store.VP875 = info.VP875;
    single_store.HPTG1 = info.HPTG1;
    single_store.VPTG1 = info.VPTG1;
    single_store.HPTG2 = info.HPTG2;
    single_store.VPTG2 = info.VPTG2;
    single_store.FOM = finalFOM;
    single_store.VP873 = info.VP873;
    single_store.VP875Offset = info.VP875Offset;
	single_store.VP873Offset = info.VP873Offset;
    single_store.HP875Offset = info.HP875Offset;
    single_store.VPTG1Offset = info.VPTG1Offset;
    single_store.HPTG1Offset = info.HPTG1Offset;
    single_store.VPTG2Offset = info.VPTG2Offset;
    single_store.HPTG2Offset = info.HPTG2Offset;
    single_store.BTJT2 = info.BTJT2;
    single_store.THCURR = info.THCURR;
    single_store.M875HS = info.M875HS;
    single_store.M875VS = info.M875VS;
    single_store.M875HM = info.M875HM;
    single_store.M875VM = info.M875VM;
    single_store.M876HS = info.M876HS;
    single_store.M876VS = info.M876VS;
    single_store.M876HM = info.M876HM;
    single_store.M876VM = info.M876VM;
    single_store.M875BB = info.M875BB;
    single_store.M876BB = info.M876BB;
    single_store.MMBTBB = info.MMBTBB;
    single_store.M875BB_spill_time_diff = info.M875BB_spill_time_diff;
    single_store.M876BB_spill_time_diff = info.M876BB_spill_time_diff;
    single_store.MMBTBB_spill_time_diff = info.MMBTBB_spill_time_diff;

    // beam-quality result: missing values (-999) become NaN, so that a cut on
    // FOM_best rejects spills without a FOM
    auto const orNaN = [](float v)
      { return (std::isfinite(v) && v != -999.f)? v: std::numeric_limits<float>::quiet_NaN(); };
    single_store.FOM_best = orNaN(info.FOM_best);
    single_store.fom_source = info.fom_source;
    single_store.fom_status = info.fom_status;
    single_store.tor_used = orNaN(info.tor_used);
    single_store.hpos_tgt = orNaN(info.hpos_tgt);
    single_store.hang_tgt = orNaN(info.hang_tgt);
    single_store.vpos_tgt = orNaN(info.vpos_tgt);
    single_store.vang_tgt = orNaN(info.vang_tgt);
    single_store.tgt_sigx = orNaN(info.tgt_sigx);
    single_store.tgt_sigy = orNaN(info.tgt_sigy);
    single_store.fom_version = info.fom_version;
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
