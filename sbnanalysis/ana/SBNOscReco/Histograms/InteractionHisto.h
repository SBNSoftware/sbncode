#ifndef _sbnanalysis_InteractionHisto_hh
#define _sbnanalysis_InteractionHisto_hh

#include <string>
#include <vector>

#include <core/Event.hh>
#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "../Data/RecoEvent.h"
#include "../Data/RecoTrack.h"
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
struct InteractionHistos : public HistoList {

  TH1D *track_length; //!< Length of the reconstructed primary track
  TH1D *nuE; //!< Neutrino energy
  TH1D *track_p; //!< Primary track momentum
  TH1D *true_deposited_energy; //!< Deposited energy of true track
  TH1D *beam_center_distance; //!< Distance of the neutrino interaction to the beam center
  TH1D *Q2; //!< Q2 of the interaction
  TH1D *true_contained_length; //!< True contained length of primary track
  TH1D *true_track_multiplicity; //!< True particle multiplicity of the interaction
  TH1D *crosses_tpc; //!< Whether the primary track crosses a TPC boundary
  TH1D *dist_to_match; //!< Distance from this vertex to the closest matching vertex reco->truth and truth->reco
  TH1D *primary_track_completion; //!< Completion of the primary track
  TH1D *n_reco_vertices; //!< Number of reconstructed vertices in the event with this vertex
  TH1D *maxpe_crt_intime_hit; //!< Maximum number of PE's in a single CRT hit in time with the beam
  TH1D *crt_hit_times;
  TH1D *closest_crt_hit_time;
  TH2D *intime_crt_hits_xy;
  TH2D *intime_crt_hits_xz;
  TH2D *intime_crt_hits_yz;
  TH2D *vertex_xy;
  TH2D *vertex_yz;
  TH2D *vertex_xz;
  TH2D *crthit_xy;
  TH2D *crthit_yz;
  TH2D *crthit_xz;
  TH2D *light_trigger;

  /**
 *  Intialize the histograms
 *  \param prefix A prefix to be added to each histogram name
 *  \param mode The mode of interaction for these histograms
 *  \param index The cut index for these histograms.
 */
  void Initialize(const std::string &prefix, const geo::BoxBoundedGeo &detector_volume, const std::vector<double> &tagger_volume);
 
  /**
 * Fill the histograms with a single interaction
 * \param vertex_index The index of this vertex into the list of truth/reco interactions
 * \param is_truth Whether this interaction is true or reco
 * \param event The reco event object
 * \param core_truth The list of true interactions from sbncode core
 */
  void Fill(
    unsigned vertex_index,
    bool is_truth,
    const numu::RecoEvent &event,
    const std::vector<event::Interaction> &core_truth);
};
  
  } // namespace SBNOSc
} // namespace ana
#endif
