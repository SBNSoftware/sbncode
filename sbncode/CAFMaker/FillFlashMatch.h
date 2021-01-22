
#ifndef CAF_FILLFLASHMATCH_H
#define CAF_FILLFLASHMATCH_H

#include <array>

#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "sbncode/StandardRecord/StandardRecord.h"


namespace caf
{

  void FillSliceFlashMatch(const anab::T0 *fmatch,
                           caf::SRSlice &srslice,
                           bool allowEmpty = false);

  void FillSliceFlashMatchA(const anab::T0 *fmatch,
			    caf::SRSlice &srslice,
			    bool allowEmpty = false);

  void FillSliceFlashMatchB(const anab::T0 *fmatch,
			    caf::SRSlice &srslice,
			    bool allowEmpty = false);
}

#endif
