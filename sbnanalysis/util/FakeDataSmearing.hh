#ifndef __sbnanalysis_util_FakeDataSmearing__
#define __sbnanalysis_util_FakeDataSmearing__

#include "core/Event.hh"
#include "core/Experiment.hh"

namespace util {

typedef enum {
  kMode1 = 1, kMode2 = 2, kMode3 = 3, kMode4 = 4, kMode5 = 5, kMode6 = 6, kMode7 = 7
} SmearingMode;

/**
 * Apply fake data smearing to the truth to build a new reco interaction
 */
void ApplyFakeDataSmearing(SmearingMode mode,
                           Experiment exp,
                           const event::Interaction* truth,
                           event::RecoInteraction* reco);


}  // namespace util

#endif  // __sbnanalysis_util_FakeDataSmearing__

