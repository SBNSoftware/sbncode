#include "Covariance.h"
#include "Chi2Sensitivity.h"

namespace ana {
  namespace SBNOsc {

Chi2Sensitivity::Chi2Sensitivity(std::vector<EventSample> samples) {
  Covariance cov(samples);
  // ...
}

  }  // namespace SBNOsc
}  // namespace ana

