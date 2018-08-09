#ifndef __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__
#define __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__

/**
 * \file Chi2Sensitivity.h
 */

#include "Covariance.h"

namespace ana {
  namespace SBNOsc {

class Chi2Sensitivity {
public:
  Chi2Sensitivity(std::vector<EventSample> samples);
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Chi2Sensitivity__

