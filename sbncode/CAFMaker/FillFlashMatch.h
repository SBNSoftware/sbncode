
#ifndef CAF_FILLFLASHMATCH_H
#define CAF_FILLFLASHMATCH_H

#include <array>

#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

//Include new flash match class
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRFlashMatch.h"

namespace caf
{

  void FillSliceFlashMatch(const sbn::SimpleFlashMatch* fmatch,
                           caf::SRSlice& srslice,
                           bool allowEmpty = false);

  void FillSliceFlashMatchA(const sbn::SimpleFlashMatch* fmatch,
                            caf::SRSlice& srslice,
                            bool allowEmpty = false);

  void FillSliceFlashMatchB(const sbn::SimpleFlashMatch* fmatch,
                            caf::SRSlice& srslice,
                            bool allowEmpty = false);


  void FillSliceFlashMatchOp(const sbn::SimpleFlashMatch* fmatchop,
                           caf::SRSlice& srslice,
                           bool allowEmpty = false);

  void FillSliceFlashMatchOpA(const sbn::SimpleFlashMatch* fmatchop,
                            caf::SRSlice& srslice,
                            bool allowEmpty = false);

  void FillSliceFlashMatchOpB(const sbn::SimpleFlashMatch* fmatchop,
                            caf::SRSlice& srslice,
                            bool allowEmpty = false);


  void FillSliceFlashMatchARA(const sbn::SimpleFlashMatch* fmatchara,
                           caf::SRSlice& srslice,
                           bool allowEmpty = false);

  void FillSliceFlashMatchARAA(const sbn::SimpleFlashMatch* fmatchara,
                            caf::SRSlice& srslice,
                            bool allowEmpty = false);

  void FillSliceFlashMatchARAB(const sbn::SimpleFlashMatch* fmatchara,
                            caf::SRSlice& srslice,
                            bool allowEmpty = false);


  void FillSliceFlashMatchOpARA(const sbn::SimpleFlashMatch* fmatchopara,
                           caf::SRSlice& srslice,
                           bool allowEmpty = false);

  void FillSliceFlashMatchOpARAA(const sbn::SimpleFlashMatch* fmatchopara,
                            caf::SRSlice& srslice,
                            bool allowEmpty = false);

  void FillSliceFlashMatchOpARAB(const sbn::SimpleFlashMatch* fmatchopara,
                            caf::SRSlice& srslice,
                            bool allowEmpty = false);
}
#endif
