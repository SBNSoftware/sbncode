
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
                           caf::SRFlashMatch& srflash,
                           bool allowEmpty = false);

}
#endif
