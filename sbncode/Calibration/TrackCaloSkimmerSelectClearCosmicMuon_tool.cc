// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "ITCSSelectionTool.h"

namespace sbn {

class TrackCaloSkimmerSelectClearCosmicMuon: public ITCSSelectionTool {
public:

  TrackCaloSkimmerSelectClearCosmicMuon(const fhicl::ParameterSet &p);
  ~TrackCaloSkimmerSelectClearCosmicMuon() {}

  bool Select(const TrackInfo &t) override;

private:
};

TrackCaloSkimmerSelectClearCosmicMuon::TrackCaloSkimmerSelectClearCosmicMuon(const fhicl::ParameterSet &p):
  ITCSSelectionTool(p)
{}

bool TrackCaloSkimmerSelectClearCosmicMuon::Select(const TrackInfo &t) {
  return t.clear_cosmic_muon;
}

DEFINE_ART_CLASS_TOOL(TrackCaloSkimmerSelectClearCosmicMuon)

} // end namespace sbn
