#ifndef _sbnumurecodata_RecoTrack_hh
#define _sbnumurecodata_RecoTrack_hh

#include <vector>

#include "TVector3.h"

#include "ana/SBNOscReco/Data/CRTMatch.h"
#include "ana/SBNOscReco/Data/TruthMatch.h"
#include "ana/SBNOscReco/Data/CaloEnergy.h"

namespace numu {

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
  float range_momentum_muon; //!< Range momentum calculation of the track using range under the assumption it is a muon [GeV].
  float range_momentum_proton; //!< Range momentum calculation of track using range using the assumption it is a proton [GeV].
  MCSFitResult mcs_muon; //!< MCS calculation result for Muon PID hypothesis
  MCSFitResult mcs_pion; //!< MCS calculation result for Pion PID hypothesis
  MCSFitResult mcs_proton; //!< MCS calculation result for Proton PID hypothesis
  MCSFitResult mcs_kaon; //!< MCS calculation result for Kaon PID hypothesis
  
  float chi2_proton; //!< Chi2 of dE/dx to proton hypothesis. Combined against all planes.
  float chi2_kaon; //!< Chi2 of dE/dx to kaon hypotheis. Combined against all planes.
  float chi2_pion; //!< Chi2 of dE/dx to pion hypotheis. Combined against all planes.
  float chi2_muon; //!< Chi2 of dE/dx to muon hypotheis. Combined agaisnt all planes.
  float min_chi2; //!< Minimum chi2 value across all hypotheses
  int pid_n_dof; //!< Number of d.o.f. in chi2 fit
  int pdgid; //!< Particle ID that minimizes chi2 
  
  float length; //!< Length of track
  float theta; //!< angle to z axis
  float phi; //!< angle about z axis

  bool crosses_tpc; //!< does it cross a tpc?
  bool is_contained; //!< is it contained in the "containment volume"?

  TVector3 start; //!< start position of track
  TVector3 end; //!< end position of track

  TrackTruth truth; //!< Truth information on track
  
  CRTMatch crt_match; //!< CRTMatch

  int ID; //!< ID/index of this track. Does not necessarily correspond to the Pandora index
  
  float stopping_chisq_start; //!< Chi2 fraction of stopping vs. not-stopping hypothesis to track start points
  float stopping_chisq_finish; //!< Chi2 fraction of stopping vs. not-stopping hypotheis to track end point.
  float proton_muon_score; // Score between 0 (proton-like) and 1 (muon-like)

  CaloEnergy calo_energy;
  CaloEnergy calE_plus_daughters;
  
  RecoTrack():
    range_momentum_muon(-1),
    range_momentum_proton(-1),
    chi2_proton(-1),
    chi2_kaon(-1),
    chi2_pion(-1),
    chi2_muon(-1),
    min_chi2(-1.5),
    pid_n_dof(-1),
    pdgid(-1),
    length(-1),
    theta(-999),
    phi(-999),
    crosses_tpc(false),
    is_contained(false),
    start(-999, -999, -999),
    end(-999, -999, -999),
    crt_match({}),
    ID(-1),
    stopping_chisq_start(-1),
    stopping_chisq_finish(-1)
    {}
};
} // namespace numu
#endif
