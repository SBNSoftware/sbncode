#ifndef _sbnumurecodata_RecoTrack_hh
#define _sbnumurecodata_RecoTrack_hh

#include <vector>

#include "TVector3.h"

#include "ana/SBNOscReco/Data/CRTMatch.h"
#include "ana/SBNOscReco/Data/FlashMatch.h"
#include "ana/SBNOscReco/Data/TruthMatch.h"

namespace numu {

/**
* Information of TPC track objects
*/
struct RecoTrack {
  float deposited_energy; //!< Energy deposited in the TPC
  float deposited_energy_max; //!< Maximum of deposited energy across the 3 planes
  float deposited_energy_avg; //!< Average of deposited energy across the 3 planes
  float deposited_energy_med; //!< Median of deposited energy across the 3 planes
  float range_momentum; //!< Momentum calculation of track using range using the PID hypothesis for this track (muon/proton) [GeV].
  float range_momentum_muon; //!< Momentum calculation of the track using range under the assumption it is a muon [GeV].
  
  float fwd_mcs_momentum; //!< MCS momentum calculation under hypothesis track is forward using PID hypotheis for this track [GeV].
  float fwd_mcs_momentum_muon; //!< MCS momentum calculation under hypotheis is forward assuming track is muon [GeV].
  float fwd_mcs_momentum_err; //!< MCS momentum calculation fit error under hypothesis track is forward using PID hypothsesis [GeV].
  float fwd_mcs_momentum_muon_err; //!< MCS momentum calculation fit error under hypotheis track is forward assuming track is muon [GeV].
  float bwd_mcs_momentum; //!< MCS momentum calculation under hypothesis track is backward using PID hypotheis for this track [GeV].
  float bwd_mcs_momentum_muon; //!< MCS momentum calculation under hypotheis is backward assuming track is muon [GeV].
  float bwd_mcs_momentum_err; //!< MCS momentum calculation fit error under hypothesis track is backward using PID hypothsesis [GeV].
  float bwd_mcs_momentum_muon_err; //!< MCS momentum calculation fit error under hypotheis track is backward assuming track is muon [GeV].
  bool mcs_is_backward; //!< Whether the MCS fit calculation believes the track is backwards
  
  float momentum; //!< Best guess at momentum [GeV]. Currently unfilled
  float energy; //!< Best guess at energy [GeV]. Currently unfilled
  
  float chi2_proton; //!< Chi2 of dE/dx to proton hypothesis. Combined against all planes.
  float chi2_kaon; //!< Chi2 of dE/dx to kaon hypotheis. Combined against all planes.
  float chi2_pion; //!< Chi2 of dE/dx to pion hypotheis. Combined against all planes.
  float chi2_muon; //!< Chi2 of dE/dx to muon hypotheis. Combined agaisnt all planes.
  float min_chi2; //!< Minimum chi2 value across all hypothesis
  int pid_n_dof; //!< Number of d.o.f. in chi2 fit
  int pdgid; //!< Particle ID that minimizes chi2 
  
  bool is_muon; //!< Whether the particle ID is a muon
  float length; //!< Length og track
  float costh; //!< cosine of angle to z axis
  bool contained_in_cryo; //!< is it contained a single cryostat?
  bool contained_in_tpc; //!< is it contained in a single TPC?
  bool crosses_tpc; //!< does it cross a tpc?
  bool is_contained; //!< is it contained in the "containment volume"?
  TVector3 start; //!< start position of track
  TVector3 end; //!< end position of track
  float dist_to_vertex; //!< Distance of track start point to interaction vertex (if it exists) [cm]
  TrackTruthMatch match; //!< Truth matching information
  
  std::vector<CRTMatch> crt_match; //!< Optional CRTMatch -- vector has 1 or 0 entries
  std::vector<FlashMatch> flash_match; //!< Optional FlashMatch -- vector has 1 or 0 entries
  int ID; //!< ID/index of this track. Does not necessarily correspond to the Pandora index
  
  float stopping_chisq_start; //!< Chi2 fraction of stopping vs. not-stopping hypothesis to track start points
  float stopping_chisq_finish; //!< Chi2 fraction of stopping vs. not-stopping hypotheis to track end point.
  std::vector<float> tpc_t0s; //!< List of T0 values assigned to track by pandora [us].
  
  RecoTrack():
    deposited_energy(-1),
    deposited_energy_max(-1),
    deposited_energy_avg(-1),
    deposited_energy_med(-1),
    range_momentum(-1),
    range_momentum_muon(-1),
    fwd_mcs_momentum(-1),
    fwd_mcs_momentum_muon(-1),
    fwd_mcs_momentum_err(-1),
    fwd_mcs_momentum_muon_err(-1),
    bwd_mcs_momentum(-1),
    bwd_mcs_momentum_muon(-1),
    bwd_mcs_momentum_err(-1),
    bwd_mcs_momentum_muon_err(-1),
    mcs_is_backward(false),
    momentum(-1),
    energy(-1),
    chi2_proton(-1),
    chi2_kaon(-1),
    chi2_pion(-1),
    chi2_muon(-1),
    min_chi2(-1.5),
    pid_n_dof(-1),
    pdgid(-1),
    is_muon(false),
    length(-1),
    costh(-999),
    contained_in_cryo(false),
    contained_in_tpc(false),
    crosses_tpc(false),
    is_contained(false),
    start(-999, -999, -999),
    end(-999, -999, -999),
    dist_to_vertex(-1),
    crt_match({}),
    flash_match({}),
    ID(-1),
    stopping_chisq_start(-1),
    stopping_chisq_finish(-1),
    tpc_t0s()
    {}
};
} // namespace numu
#endif
