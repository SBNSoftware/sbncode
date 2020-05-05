#ifndef _sbnumurecodata_PIDData_hh
#define _sbnumurecodata_PIDData_hh

namespace numu::flat {
struct PIDData {
  int true_pdg;
  float trueP;
  float purity;
  float completion;
  float length;
  float theta;
  float phi;
  float chi2_muon;
  float chi2_proton;
  float chi2_pion;
  float crt_hit_distance;
  bool contained;
  int n_trk_daughters;
  int n_shr_daughters;
};

} // end namespace
#endif
