#ifndef _sbnumurecodata_RecoTrack_hh
#define _sbnumurecodata_RecoTrack_hh

#include <vector>

#include "TVector3.h"

#include "ana/SBNOscReco/Data/CRTMatch.h"
#include "ana/SBNOscReco/Data/FlashMatch.h"
#include "ana/SBNOscReco/Data/TruthMatch.h"

namespace numu {

enum Wall {
  wNone=0,
  wTop=1,
  wBottom=2,
  wLeft=3,
  wRight=4,
  wFront=5,
  wBack=6
};

struct MCSFitResult {
  float fwd_mcs_momentum; //!< MCS momentum calculation under hypothesis track is forward
  float fwd_mcs_momentum_err; //!< MCS momentum calculation fit error under hypothesis track is forward
  float bwd_mcs_momentum; //!< MCS momentum calculation under hypothesis track is backward
  float bwd_mcs_momentum_err; //!< MCS momentum calculation fit error under hypothesis track is backward
  bool mcs_is_backward; //!< Whether the MCS fit calculation believes the track is backwards 

  MCSFitResult():
    fwd_mcs_momentum(-1),
    fwd_mcs_momentum_err(-1),
    bwd_mcs_momentum(-1),
    bwd_mcs_momentum_err(-1),
    mcs_is_backward(false)
  {}
};

/**
* Information of TPC track objects
*/
struct RecoTrack {
  float deposited_energy; //!< Energy deposited in the TPC
  float deposited_energy_max; //!< Maximum of deposited energy across the 3 planes
  float deposited_energy_avg; //!< Average of deposited energy across the 3 planes
  float deposited_energy_med; //!< Median of deposited energy across the 3 planes
  float range_momentum_muon; //!< Range momentum calculation of the track using range under the assumption it is a muon [GeV].
  float range_momentum_proton; //!< Range momentum calculation of track using range using the assumption it is a proton [GeV].
  MCSFitResult mcs_muon; //!< MCS calculation result for Muon PID hypothesis
  MCSFitResult mcs_pion; //!< MCS calculation result for Pion PID hypothesis
  MCSFitResult mcs_proton; //!< MCS calculation result for Proton PID hypothesis
  MCSFitResult mcs_kaon; //!< MCS calculation result for Kaon PID hypothesis
  
  float mean_trucated_dQdx; //!< Mean of dQdx values inside the standard deviation
  
  float momentum; //!< Best guess at momentum [GeV]. 
  float range_momentum; //!< momentum calculated from range method
  float mcs_momentum; //!< momentum calculated from mcs method
  float energy; //!< Best guess at energy [GeV]. 
  
  float chi2_proton; //!< Chi2 of dE/dx to proton hypothesis. Combined against all planes.
  float chi2_kaon; //!< Chi2 of dE/dx to kaon hypotheis. Combined against all planes.
  float chi2_pion; //!< Chi2 of dE/dx to pion hypotheis. Combined against all planes.
  float chi2_muon; //!< Chi2 of dE/dx to muon hypotheis. Combined agaisnt all planes.
  float min_chi2; //!< Minimum chi2 value across all hypothesis
  int pid_n_dof; //!< Number of d.o.f. in chi2 fit
  int pdgid; //!< Particle ID that minimizes chi2 
  
  Wall wall_enter; //!< the face of the TPC that the track crosses on enter
  Wall wall_exit; //!< the face of the TPC that the track crosses on exit
  bool is_muon; //!< Whether the particle ID is a muon
  float length; //!< Length of track
  float artlength; //!< Test to check default length from pandora
  float costh; //!< cosine of angle to z axis
  bool contained_in_cryo; //!< is it contained a single cryostat?
  bool contained_in_tpc; //!< is it contained in a single TPC?
  bool crosses_tpc; //!< does it cross a tpc?
  bool is_contained; //!< is it contained in the "containment volume"?
  TVector3 start; //!< start position of track
  float start_time; //!< start time of track
  float end_time; //!< end time of track 
  TVector3 end; //!< end position of track
  float dist_to_vertex; //!< Distance of track start point to interaction vertex (if it exists) [cm]
  TrackTruthMatch match; //!< Truth matching information
  
  CRTMatch crt_match; //!< CRTMatch
  FlashMatch flash_match; //!< Flash matching info
  int ID; //!< ID/index of this track. Does not necessarily correspond to the Pandora index
  
  float stopping_chisq_start; //!< Chi2 fraction of stopping vs. not-stopping hypothesis to track start points
  float stopping_chisq_finish; //!< Chi2 fraction of stopping vs. not-stopping hypotheis to track end point.
  std::vector<float> tpc_t0s; //!< List of T0 values assigned to track by pandora [us].
  
  RecoTrack():
    deposited_energy(-1),
    deposited_energy_max(-1),
    deposited_energy_avg(-1),
    deposited_energy_med(-1),
    range_momentum_muon(-1),
    range_momentum_proton(-1),
    mean_trucated_dQdx(-9999.),
    momentum(-1),
    range_momentum(-1),
    mcs_momentum(-1),
    energy(-1),
    chi2_proton(-1),
    chi2_kaon(-1),
    chi2_pion(-1),
    chi2_muon(-1),
    min_chi2(-1.5),
    pid_n_dof(-1),
    pdgid(-1),
    wall_enter(numu::wNone),
    wall_exit(numu::wNone),
    is_muon(false),
    length(-1),
    artlength(-1),
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
