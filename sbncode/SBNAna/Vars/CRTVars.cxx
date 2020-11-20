#include "SBNAna/Vars/CRTVars.h"

#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  const SpillMultiVar kCRTHit_Times ([](const caf::SRSpillProxy *sr)
      {
        std::vector<double> crtHitTimes;
        for (int crtHitIter=0; crtHitIter<sr->ncrt_hits; crtHitIter++){
          const auto& crtHit(sr->crt_hits[crtHitIter]);
          crtHitTimes.push_back(crtHit.time);
        }
        return crtHitTimes;
      });

  const SpillMultiVar kCRTHit_PEs ([](const caf::SRSpillProxy *sr)
      {
        std::vector<double> crtHitPE;
        for (int crtHitIter=0; crtHitIter<sr->ncrt_hits; crtHitIter++){
          const auto& crtHit(sr->crt_hits[crtHitIter]);
          crtHitPE.push_back(crtHit.pe);
        }
        return crtHitPE;
      });

  const SpillMultiVar kCRTHit_X ([](const caf::SRSpillProxy *sr)
      {
        std::vector<double> crtHitX;
        for (int crtHitIter=0; crtHitIter<sr->ncrt_hits; crtHitIter++){
          const auto& crtHit(sr->crt_hits[crtHitIter]);
          crtHitX.push_back(crtHit.position.x);
        }
        return crtHitX;
      });

  const SpillMultiVar kCRTHit_Y ([](const caf::SRSpillProxy *sr)
      {
        std::vector<double> crtHitY;
        for (int crtHitIter=0; crtHitIter<sr->ncrt_hits; crtHitIter++){
          const auto& crtHit(sr->crt_hits[crtHitIter]);
          crtHitY.push_back(crtHit.position.x);
        }
        return crtHitY;
      });

  const SpillMultiVar kCRTHit_Z ([](const caf::SRSpillProxy *sr)
      {
        std::vector<double> crtHitZ;
        for (int crtHitIter=0; crtHitIter<sr->ncrt_hits; crtHitIter++){
          const auto& crtHit(sr->crt_hits[crtHitIter]);
          crtHitZ.push_back(crtHit.position.x);
        }
        return crtHitZ;
      });

}

