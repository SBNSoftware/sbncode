
/**
 * @file   sbncode/CAFMaker/FillBlip.h
 * @brief  CAFMaker utilities for filling "blip" objects into CAF.
 * @author Jacob McLaughlin (jmclaughlin2@illinoistech.edu)
 */
#ifndef CAF_FILLBLIP_H
#define CAF_FILLBLIP_H
#include "sbnanaobj/StandardRecord/SRBlip.h"
#include "sbnobj/SBND/Blip/BlipDataTypes.h"

#include <vector>

namespace caf
{

  void FillBlip(   const std::vector<blip::Blip>& LAr_Blips,  std::vector<caf::SRBlip>& CAF_Blips);
  void FillMCTruthBlip( blip::TrueBlip const & TrueLAr_Blip, caf::SRTrueBlip &TrueCAF_Blip );
  void FillBlipRealtedHitCluster(blip::HitClust const & LAr_HitClust, caf::SRBlipHitClust &CAF_HitClust);
}

#endif
