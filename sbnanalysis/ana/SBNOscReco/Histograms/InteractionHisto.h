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
class TProfile;

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
  TH1D *crt_pes;
  TH1D *fmatch_score;
  TH2D *fmatch_score_true_time;
  TH2D *fmatch_score_true_time_zoom;
  TH1D *fmatch_score_outtime;
  TH1D *fmatch_score_intime;
  TH2D *fmatch_time_true_time_zoom;
  TH1D *fmatch_time;
  TH1D *fmatch_time_real_time;
  TH2D *intime_crt_hits_xy;
  TH2D *intime_crt_hits_xz;
  TH2D *intime_crt_hits_yz;
  TH2D *vertex_xy;
  TH2D *vertex_yz;
  TH2D *vertex_xz;
  TH1D *vertex_x;
  TH1D *vertex_y;
  TH1D *vertex_z;
  TH2D *light_trigger;
  TH1D *nu_score;
  TProfile *nu_score_nu_energy;
  TProfile *nu_score_muon_length;
  TH2D *finalstate_tr;

  /**
 *  Intialize the histograms
 *  \param prefix A prefix to be added to each histogram name
 *  \param mode The mode of interaction for these histograms
 *  \param index The cut index for these histograms.
 */
  void Initialize(const std::string &prefix, const geo::BoxBoundedGeo &detector_volume, const std::vector<double> &tagger_volume);
 
  /**
 * Fill the histograms with a single reconstructed neutrino candidate
 * \param vertex The reconstructed vertex being filled
 * \param event The reco event object
 * \param core_truth The list of true interactions from sbncode core
 */
  void Fill(
    const numu::RecoInteraction &vertex,
    const numu::RecoEvent &event,
    const std::vector<event::Interaction> &core_truth);

  /**
 * Fill the histograms with a single true interaction
 * \param interaction The true interaction
 * \param mctruth_id ID of the interaction (index into the list of Interactions)
 * \param event The reco event object
 */
  void Fill(const event::Interaction &interaction, unsigned mctruth_id, const numu::RecoEvent &event);

private:
  void FillEvent(const numu::RecoEvent &event);
};
  
  } // namespace SBNOSc
} // namespace ana
#endif
