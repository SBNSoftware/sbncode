// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "ITCSSelectionTool.h"

namespace sbn {

class TrackCaloSkimmerSelectAllTrack: public ITCSSelectionTool {
public:

  TrackCaloSkimmerSelectAllTrack(const fhicl::ParameterSet &p);
  ~TrackCaloSkimmerSelectAllTrack() {}

  bool Select(const TrackInfo &t) override;

private:
};

TrackCaloSkimmerSelectAllTrack::TrackCaloSkimmerSelectAllTrack(const fhicl::ParameterSet &p):
  ITCSSelectionTool(p)
{}

bool TrackCaloSkimmerSelectAllTrack::Select(const TrackInfo &t) {
  return true;
}

DEFINE_ART_CLASS_TOOL(TrackCaloSkimmerSelectAllTrack)

} // end namespace sbn
