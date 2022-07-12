#include "sbnobj/Common/Calibration/TrackCaloSkimmerObj.h"
#include "fhiclcpp/ParameterSet.h"

#include <stdlib.h>
#include <time.h> 

namespace sbn {

class ITCSSelectionTool {
public:
  ITCSSelectionTool(const fhicl::ParameterSet &p):
    fAllowT0(p.get<std::vector<bool>>("AllowT0", {})),
    fInvert(p.get<bool>("Invert", false)),
    fNPreScale(p.get<unsigned>("PreScale", 1)),
    fISelect(0) 
  {
      // Set ISelect to a random number
      // so that events get prescaled uniformly across runs
      srand (time(NULL)); // random seed
      fISelect = rand();
  }

  virtual ~ITCSSelectionTool() noexcept = default;

  /// For external modules to call: run the actual selection
  bool DoSelect(const TrackInfo &t) {

    // Check if the T0 is allowed
    if (!fAllowT0.empty()) {
      if (t.whicht0 >= 0 && ((unsigned)t.whicht0 >= fAllowT0.size() || !fAllowT0[t.whicht0])) {
        return false;
      }
    }

    bool selected = Select(t) == !fInvert /* implement inversion */;

    // update index
    fISelect += selected;

    return selected && (fISelect % fNPreScale == 0) /* implement prescale */;
  }

  unsigned GetPrescale() const { return fNPreScale; }

protected:
  /// For children to implement: Whether to select a given track
  virtual bool Select(const TrackInfo &t) = 0;

  std::vector<bool> fAllowT0;
  bool fInvert;
  unsigned fNPreScale;
  unsigned fISelect;
};

} // end namespace sbn
