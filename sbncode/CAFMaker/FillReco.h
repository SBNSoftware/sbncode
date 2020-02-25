
#ifndef CAF_FILLRECO_H
#define CAF_FILLRECO_H

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"

#include "sbncode/StandardRecord/SRSlice.h"
#include "sbncode/StandardRecord/StandardRecord.h"

namespace caf
{
  void FillSliceVars(const recob::Slice& slice,
                     caf::SRSlice& srslice,
                     bool allowEmpty = false);

  void FillTrackVars(const recob::Track& track,
                     caf::SRTrack& srtrk,
                     bool allowEmpty = false);
}

#endif
