// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "ITCSSelectionTool.h"

namespace sbn {

class TrackCaloSkimmerSelectAnocde2CathodeTrack: public ITCSSelectionTool {
public:

  TrackCaloSkimmerSelectAnocde2CathodeTrack(const fhicl::ParameterSet &p);
  ~TrackCaloSkimmerSelectAnocde2CathodeTrack() {}

  bool Select(const TrackInfo &t) override;

private:
  // config
  double fTickCut;
};

TrackCaloSkimmerSelectAnocde2CathodeTrack::TrackCaloSkimmerSelectAnocde2CathodeTrack(const fhicl::ParameterSet &p):
  fTickCut(p.get<double>("TickCut"))
{}

bool TrackCaloSkimmerSelectAnocde2CathodeTrack::Select(const TrackInfo &t) {
  // use the collection plane
  return std::max(abs(t.hit_max_time_p2_tpcE - t.hit_min_time_p2_tpcE), abs(t.hit_max_time_p2_tpcW - t.hit_min_time_p2_tpcW)) > fTickCut;

}

DEFINE_ART_CLASS_TOOL(TrackCaloSkimmerSelectAnocde2CathodeTrack)

} // end namespace sbn
