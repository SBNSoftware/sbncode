// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "ITCSSelectionTool.h"

namespace sbn {

class TrackCaloSkimmerSelectAnode2CathodeTrack: public ITCSSelectionTool {
public:

  TrackCaloSkimmerSelectAnode2CathodeTrack(const fhicl::ParameterSet &p);
  ~TrackCaloSkimmerSelectAnode2CathodeTrack() {}

  bool Select(const TrackInfo &t) override;

private:
  // config
  double fTickCut;
};

TrackCaloSkimmerSelectAnode2CathodeTrack::TrackCaloSkimmerSelectAnode2CathodeTrack(const fhicl::ParameterSet &p):
  ITCSSelectionTool(p),
  fTickCut(p.get<double>("TickCut"))
{}

bool TrackCaloSkimmerSelectAnode2CathodeTrack::Select(const TrackInfo &t) {
  // use the collection plane
  return std::max(abs(t.hit_max_time_p2_tpcE - t.hit_min_time_p2_tpcE), abs(t.hit_max_time_p2_tpcW - t.hit_min_time_p2_tpcW)) > fTickCut;

}

DEFINE_ART_CLASS_TOOL(TrackCaloSkimmerSelectAnode2CathodeTrack)

} // end namespace sbn
