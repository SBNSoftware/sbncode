#include "Cuts.h"
#include "NumuRecoSelection.h"

namespace ana {
 namespace SBNOsc {

void Cuts::Initialize(double trackMatchCompletionCut) {
  fConfig.trackMatchCompletionCut = trackMatchCompletionCut;
}

std::array<bool, Cuts::nCuts> Cuts::ProcessRecoCuts(const NumuRecoSelection::RecoEvent &event, unsigned reco_vertex_index) {
  bool is_reco = true;

  // require close to truth
  bool v_quality = event.reco[reco_vertex_index].match.event_vertex_id >= 0 && event.reco[reco_vertex_index].match.truth_vertex_distance < 10.;

  bool t_quality = false;
  int t_mcparticle_id = event.reco[reco_vertex_index].slice.tracks.at(event.reco[reco_vertex_index].slice.primary_track_index).match.mcparticle_id;
  for (const NumuRecoSelection::RecoInteraction &truth: event.truth) {
    if (truth.slice.primary_track_index >= 0 && 
      truth.slice.tracks.at(truth.slice.primary_track_index).match.mcparticle_id == t_mcparticle_id) {
      t_quality = true;
      break;
    } 
  }
  // require completion
  t_quality = t_quality && (fConfig.trackMatchCompletionCut < 0 ||
    event.reco[reco_vertex_index].slice.tracks.at(event.reco[reco_vertex_index].slice.primary_track_index).match.completion > fConfig.trackMatchCompletionCut);

  bool is_contained = event.reco[reco_vertex_index].slice.tracks.at(event.reco[reco_vertex_index].slice.primary_track_index).is_contained;

  return {
    is_reco,
    v_quality,
    t_quality && v_quality,
    is_contained
  };
}

/*
bool Cuts::SelectReco(std::array<bool, Cuts::nCuts> &cuts) {
  return 
    cuts[0] && 
    (cuts[1] || !fConfig.requireTrack) &&
    (cuts[3] || !fConfig.requireMatched) &&
    (cuts[4] || !fConfig.requireContained);
}
*/
  }
}
