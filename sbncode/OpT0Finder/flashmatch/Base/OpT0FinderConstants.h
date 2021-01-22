#ifndef OPT0FINDER_OPT0FINDERCONSTANTS_H
#define OPT0FINDER_OPT0FINDERCONSTANTS_H

#include <utility>
#include <climits>
#include <limits>
#include <cstdlib>

namespace flashmatch {

  /// Utility: invalid value for double 
  const double kINVALID_DOUBLE = std::numeric_limits<double>::max();

  /// Utility: invalid value for int
  const int    kINVALID_INT = std::numeric_limits<int>::max();

  /// Utility: invalid value for size
  const size_t kINVALID_SIZE = std::numeric_limits<size_t>::max();

} 

#endif
