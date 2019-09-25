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
  double deposited_energy_max;
  double deposited_energy_avg;
  double deposited_energy_med;
  double range_momentum;
  double range_momentum_muon;
  
  double fwd_mcs_momentum;
  double fwd_mcs_momentum_muon;
  double fwd_mcs_momentum_err;
  double fwd_mcs_momentum_muon_err;
  double bwd_mcs_momentum;
  double bwd_mcs_momentum_muon;
  double bwd_mcs_momentum_err;
  double bwd_mcs_momentum_muon_err;
  bool mcs_is_backward;
  
  // buest guess at momentum 
  double momentum;
  double energy;
  
  double chi2_proton;
  double chi2_kaon;
  double chi2_pion;
  double chi2_muon;
  double min_chi2;
  int pid_n_dof;
  int pdgid; //!< particle id
  
  bool is_muon;
  double length;
  double costh; //!< cosine of angle to z axis
  bool contained_in_cryo; //!< is it contained a single cryostat?
  bool contained_in_tpc; //!< is it contained in a single TPC?
  bool crosses_tpc; //!< does it cross a tpc?
  bool is_contained; //!< is it contained in the "containment volume"?
  TVector3 start; //!< start position of track
  TVector3 end; //!< end position of track
  double dist_to_vertex;
  TrackTruthMatch match;
  
  std::vector<CRTMatch> crt_match; 
  std::vector<FlashMatch> flash_match;
  int ID;
  
  double stopping_chisq_start;
  double stopping_chisq_finish;
  std::vector<double> tpc_t0s;
  
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
    // std::vector<double> calo_dEdx;
    // std::vector<double> calo_extent;
};
} // namespace numu
#endif
