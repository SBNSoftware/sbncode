#ifndef CAF_FILLBLIP_H
#define CAF_FILLBLIP_H
#include<iostream>
#include "sbnanaobj/StandardRecord/SRBlip.h"
#include "sbndcode/BlipRecoSBND/Utils/BlipUtils.h"

#include <vector>

namespace caf
{

  void FillBlip(   const std::vector<blip::Blip>& LarBlips,  std::vector<caf::SRBlip>& CAF_Blips);
  void FillMCTruthBlip(blip::Blip& LarBlip,  caf::SRBlip& CAF_Blip );
  void FillBlipRealtedHitCluster(blip::Blip& LarBlip,  caf::SRBlip& CAF_Blip);
}

#endif
