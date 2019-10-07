#ifndef _sbnumurecodata_RecoTrack_hh
#define _sbnumurecodata_RecoTrack_hh

#include <vector>

#include "TVector3.h"

#include "CRTMatch.h"
#include "FlashMatch.h"
#include "TruthMatch.h"

namespace numu {

struct RecoTrack {
  // track specific info
  float deposited_energy_max;
  float deposited_energy_avg;
  float deposited_energy_med;
  float range_momentum;
  float range_momentum_muon;
  
  float fwd_mcs_momentum;
  float fwd_mcs_momentum_muon;
  float fwd_mcs_momentum_err;
  float fwd_mcs_momentum_muon_err;
  float bwd_mcs_momentum;
  float bwd_mcs_momentum_muon;
  float bwd_mcs_momentum_err;
  float bwd_mcs_momentum_muon_err;
  bool mcs_is_backward;
  
  // buest guess at momentum 
  float momentum;
  float energy;
  
  float chi2_proton;
  float chi2_kaon;
  float chi2_pion;
  float chi2_muon;
  float min_chi2;
  int pid_n_dof;
  int pdgid; //!< particle id
  
  bool is_muon;
  float length;
  float costh; //!< cosine of angle to z axis
  bool contained_in_cryo; //!< is it contained a single cryostat?
  bool contained_in_tpc; //!< is it contained in a single TPC?
  bool crosses_tpc; //!< does it cross a tpc?
  bool is_contained; //!< is it contained in the "containment volume"?
  TVector3 start; //!< start position of track
  TVector3 end; //!< end position of track
  float dist_to_vertex;
  TrackTruthMatch match;
  
  std::vector<CRTMatch> crt_match; 
  std::vector<FlashMatch> flash_match;
  int ID;
  
  float stopping_chisq_start;
  float stopping_chisq_finish;
  std::vector<float> tpc_t0s;
  
  RecoTrack():
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
    is_muon(-1),
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
    
    // More involved info -- need this later???
    // std::vector<TLorentzVector> trajectory;
    // std::vector<float> calo_dEdx;
    // std::vector<float> calo_extent;
};
} // namespace numu
#endif
