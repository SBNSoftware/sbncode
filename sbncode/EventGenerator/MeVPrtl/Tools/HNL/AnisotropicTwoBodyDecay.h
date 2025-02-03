#ifndef AnisotropicTwoBodyDecay_h
#define AnisotropicTwoBodyDecay_h

#include "TLorentzVector.h"
#include "TRandom3.h"


//Implementation of HNL Two Body Decays Anisotropies (N -> m+l-)
//Valid as long as the HNL is decaying into a charged lepton and a charged pseudoscalar meson
//See arXiv:1905.00284 for more details @LuisPelegrina

namespace evgen {
namespace ldm {
namespace AnTwoBD {

  TVector3 RandomUnitVector(TRandom3 *rand);
  void RotateToHNLPolFrame(TVector3 p_HNL, TVector3 &p_l, TVector3 &p_m);
  double two_body_momentum(double parent_mass, double childA_mass, double childB_mass);
  double lambda(double a, double b, double c);
  double I_1_plus(double x,double y, double cos_theta);
  double I_1_minus(double x,double y, double cos_theta);
  double GetMaxIntegral(int LeptonPDG, double Pol, double eps_m, double eps_l);
  void AnisotropicTwoBodyDist( TLorentzVector pHNL, TLorentzVector &pl, TLorentzVector &pm, double m_HNL, int LeptonPDG, int MesonPDG, double Pol);    

} // namespace AnTwoBD
} // namespace ldm
} // namespace evgen

#endif
