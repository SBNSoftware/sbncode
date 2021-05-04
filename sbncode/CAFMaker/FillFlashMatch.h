
#ifndef CAF_FILLFLASHMATCH_H
#define CAF_FILLFLASHMATCH_H

#include <array>

#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataobj/AnalysisBase/T0.h"
//Include new flash match class
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"

namespace caf
{

  void FillSliceFlashMatch(const sbn::SimpleFlashMatch &fmatch,
                           caf::SRSlice &srslice,
                           bool allowEmpty = false);

  void FillSliceFlashMatchA(const sbn::SimpleFlashMatch &fmatch,
			    caf::SRSlice &srslice,
			    bool allowEmpty = false);

  void FillSliceFlashMatchB(const sbn::SimpleFlashMatch &fmatch,
			    caf::SRSlice &srslice,
			    bool allowEmpty = false);
}

#endif
