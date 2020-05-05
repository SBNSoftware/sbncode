#ifndef _sbnanalysis_TrueParticleHisto_hh
#define _sbnanalysis_TrueParticleHisto_hh

#include <string>

#include "../Data/TrueParticle.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "HistoList.h"

class TH1D;
class TH2D;
namespace ana {
 namespace SBNOsc {

/**
 * Histograms associated with truth information
 * for tracks
 */
struct TrueParticleHistos : public HistoList {

  TH1D *momentum;
  TH1D *deposited_energy;
  TH1D *length;
  TH1D *kinetic_energy;
  TH1D *frac_deposited;
  TH1D *theta;

  TH1D *vertex_x;
  TH1D *vertex_y;
  TH1D *vertex_z;

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
  void Fill(const numu::TrueParticle &particle);
};
  
  } // namespace SBNOSc
} // namespace ana
#endif
