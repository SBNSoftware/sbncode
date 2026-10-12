
#ifndef CAF_FILLFLASHMATCH_H
#define CAF_FILLFLASHMATCH_H

//Include new flash match class
#include "sbncode/CAFMaker/TimeRefShifters.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnanaobj/StandardRecord/SRFlashMatch.h"

namespace caf
{

  void FillSliceFlashMatch(const sbn::SimpleFlashMatch* fmatch,
                           caf::SRFlashMatch& srflash,
                           caf::TimeRefShifter<> const& shifter = {},
                           bool allowEmpty = false);

}
#endif
