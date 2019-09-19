#ifndef __sbnanalysis_ana_AutoVetoAnalysis_AutoVetoTools__
#define __sbnanalysis_ana_AutoVetoAnalysis_AutoVetoTools__

#include <TVector3.h>
#include "core/ProviderManager.hh"

/**
 * \file AutoVetoTools.h
 *
 * This is some auxiliary code that is not a selection, but does a piece
 * of the analysis. We can define any number of other functions, classes,
 * etc. which we use in the selection.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

namespace ana {
  namespace AutoVetoAnalysis {

      //geo::CryostatGeo& fCryo0;// = fGeometryService->Cryostat(0);
      //geo::CryostatGeo& fCryo1;// = fGeometryService->Cryostat(1);
      //geo::TPCGeo& fTpc00;
      //geo::TPCGeo& fTpc01;
      //geo::TPCGeo& fTpc10;
      //geo::TPCGeo& fTpc11;
      bool IsAV(const core::ProviderManager &manager, const TVector3 &point);
      bool IsFV(const core::ProviderManager &manager, const TVector3 &point, const TVector3 &fid1, const TVector3 &fid2);
      int MacToADReg(int mac);

  }  // namespace ExampleAnalysis
}  // namespace ana

#endif  // __sbnanalysis_ana_ExampleAnalysis_ExampleTools__

