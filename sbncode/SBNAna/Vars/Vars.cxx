#include "SBNAna/Vars/Vars.h"

#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <vector>

namespace ana
{

  const SpillVar kRun = SIMPLESPILLVAR(hdr.run);
  const SpillVar kEvt = SIMPLESPILLVAR(hdr.evt);
  //  const Var kSlc = SIMPLEVAR(hdr.subevt);

  const Var kCounting = kUnweighted;
  const SpillVar kSpillCounting = kSpillUnweighted;


  // There can be more than one crt hit per spill so the variables are vectors.
  // Use SpillMultiVar for now.
  const SpillMultiVar kCRTHitX([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> positions;
    for(const auto& hit : sr->crt_hits){
      positions.push_back(hit.position.x);
    }
    return positions;
  });

  const SpillMultiVar kCRTHitY([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> positions;
    for(const auto& hit : sr->crt_hits){
      positions.push_back(hit.position.y);
    }
    return positions;
  });

  const SpillMultiVar kCRTHitZ([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> positions;
    for(const auto& hit : sr->crt_hits){
      positions.push_back(hit.position.z);
    }
    return positions;
  });

  const SpillMultiVar kCRTHitPE([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> pes;
    for(const auto& hit : sr->crt_hits){
      pes.push_back(hit.pe);
    }
    return pes;
  });

  const SpillMultiVar kCRTHitTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> times;
    for(const auto& hit : sr->crt_hits){
      times.push_back(hit.time);
    }
    return times;
  });

  // // For when we have spill beam mode info
  // const SpillCut kIsRHC([](const caf::SRProxy* sr) {return sr->spill.isRHC;});

  const Var kSlcVtxX([](const caf::SRSliceProxy *slc)
       {
         return slc->vertex.x;
       });

  const Var kSlcVtxY([](const caf::SRSliceProxy *slc)
       {
         return slc->vertex.y;
       });

  const Var kSlcVtxZ([](const caf::SRSliceProxy *slc)
       {
         return slc->vertex.z;
       });

  const Var kSlcNuScore([](const caf::SRSliceProxy *slc)
       {
         return slc->nu_score;
       });

  const Var kSlcHasFlash([](const caf::SRSliceProxy *slc)
       {
         return slc->fmatch.present;
       });

  const Var kSlcFlashScore([](const caf::SRSliceProxy *slc)
       {
         return ((bool)kSlcHasFlash(slc) ? (float)slc->fmatch.score : -5.f);
       });

  const Var kSlcIsRecoNu([](const caf::SRSliceProxy *slc)
       {
         return !slc->is_clear_cosmic;
       });
}
