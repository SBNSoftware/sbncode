#ifndef _sbnanalysis_CosmicHisto_hh
#define _sbnanalysis_CosmicHisto_hh

#include <string>
#include <vector>
#include <map>

#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "../Data/TrueParticle.h"
#include "../Data/Mode.h"

#include "HistoList.h"

class TH1D;
class TH2D;
namespace ana {
 namespace SBNOsc {

/**
 *  Histograms associated with neutrino interactions. Filled for the list of all true
 *  and reco vertices. These histograms are constructed per interaction mode per cut.
 */
struct CosmicHistos : public HistoList {

  TH1D *momentum;
  TH1D *enter_time;
  TH1D *enter_time_zoom;
  TH1D *enter_y;
  TH1D *enter_x;
  TH1D *enter_z;

  /**
 *  Intialize the histograms
 *  \param prefix A prefix to be added to each histogram name
 */
  void Initialize(const std::string &prefix, const geo::BoxBoundedGeo &detector_volume);
 
  /**
 * Fill the histograms with all cosmic muons in a single event
 *
 * \param true_particles The map of true particles in the event
 */
  void Fill(const std::map<size_t, numu::TrueParticle> &true_particles);
};
  
  } // namespace SBNOSc
} // namespace ana
#endif
