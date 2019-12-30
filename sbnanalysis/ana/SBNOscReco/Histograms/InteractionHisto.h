#ifndef _sbnanalysis_InteractionHisto_hh
#define _sbnanalysis_InteractionHisto_hh

#include <string>
#include <vector>

#include <core/Event.hh>
#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "../Data/RecoEvent.h"
#include "../Data/RecoTrack.h"
#include "../Data/Mode.h"
#include "../MultiThread/THShared.h"

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

  TH1Shared track_length; //!< Length of the reconstructed primary track
  TH1Shared nuE; //!< Neutrino energy
  TH1Shared track_p; //!< Primary track momentum
  TH1Shared true_deposited_energy; //!< Deposited energy of true track
  TH1Shared beam_center_distance; //!< Distance of the neutrino interaction to the beam center
  TH1Shared Q2; //!< Q2 of the interaction
  TH1Shared true_contained_length; //!< True contained length of primary track
  TH1Shared true_track_multiplicity; //!< True particle multiplicity of the interaction
  TH1Shared crosses_tpc; //!< Whether the primary track crosses a TPC boundary
  TH1Shared dist_to_match; //!< Distance from this vertex to the closest matching vertex reco->truth and truth->reco
  TH1Shared primary_track_completion; //!< Completion of the primary track
  TH1Shared n_reco_vertices; //!< Number of reconstructed vertices in the event with this vertex
  TH1Shared maxpe_crt_intime_hit; //!< Maximum number of PE's in a single CRT hit in time with the beam
  TH1Shared crt_hit_times;
  TH1Shared closest_crt_hit_time;
  TH1Shared crt_pes;
  TH1Shared fmatch_score;
  TH2Shared fmatch_score_true_time;
  TH2Shared fmatch_score_true_time_zoom;
  TH1Shared fmatch_score_outtime;
  TH1Shared fmatch_score_intime;
  TH2Shared fmatch_time_true_time_zoom;
  TH1Shared fmatch_time;
  TH1Shared fmatch_time_real_time;
  TH2Shared intime_crt_hits_xy;
  TH2Shared intime_crt_hits_xz;
  TH2Shared intime_crt_hits_yz;
  TH2Shared vertex_xy;
  TH2Shared vertex_yz;
  TH2Shared vertex_xz;
  TH2Shared light_trigger;

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
