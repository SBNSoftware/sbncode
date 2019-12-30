#ifndef _sbnanalysis_CosmicHisto_hh
#define _sbnanalysis_CosmicHisto_hh

#include <string>
#include <vector>
#include <map>

#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "../Data/RecoTrack.h"
#include "../Data/Mode.h"
#include "../MultiThread/THShared.h"

#include "HistoList.h"

namespace ana {
 namespace SBNOsc {

/**
 *  Histograms associated with neutrino interactions. Filled for the list of all true
 *  and reco vertices. These histograms are constructed per interaction mode per cut.
 */
struct CosmicHistos : public HistoList {

  TH1Shared momentum;
  TH1Shared enter_time;
  TH1Shared enter_time_zoom;
  TH1Shared enter_y;
  TH1Shared enter_x;
  TH1Shared enter_z;

  /**
 *  Intialize the histograms
 *  \param prefix A prefix to be added to each histogram name
 */
  void Initialize(const std::string &prefix, const geo::BoxBoundedGeo &detector_volume);
 
  /**
 * Fill the histograms with a single interaction
 * \param vertex_index The index of this vertex into the list of truth/reco interactions
 * \param is_truth Whether this interaction is true or reco
 * \param event The reco event object
 * \param core_truth The list of true interactions from sbncode core
 */
  void Fill(const std::vector<size_t> &cosmic_tracks, const std::map<size_t, numu::RecoTrack> &true_tracks);
};
  
  } // namespace SBNOSc
} // namespace ana
#endif
