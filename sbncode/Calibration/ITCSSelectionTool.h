#include "sbnobj/Common/Calibration/TrackCaloSkimmerObj.h"

namespace sbn {

class ITCSSelectionTool {
public:
  virtual ~ITCSSelectionTool() noexcept = default;

  virtual bool Select(const TrackInfo &t) = 0;
};

} // end namespace sbn
