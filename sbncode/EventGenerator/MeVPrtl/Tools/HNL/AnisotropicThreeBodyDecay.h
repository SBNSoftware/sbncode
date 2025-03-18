#ifndef AnisotropicThreeBodyDecay_h
#define AnisotropicThreeBodyDecay_h

#include "TLorentzVector.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlStage.h"

//Implementation of HNL Three Body Decays Anisotropies
//Valid as long as the HNL is decaying into a neutrino and identical final-state charged leptons
//See arXiv:2104.05719 for more details @LuisPelegrina


namespace evgen {
namespace ldm {
namespace AnThreeBD {
  
  double GL(bool CC);
  double GR();
  void CMgLgR(double (&CM)[6], double Ue4, double Umu4, double Ut4, int LeptonPDG);
  void CDAgLgR(double (&CD)[6], double Ue4, double Umu4, double Ut4, int LeptonPDG);
  void CDNgLgR(double (&CD)[6], double Ue4,double Umu4, double Ut4, int LeptonPDG);
  void KinDep(double (&K)[6], double Pol, double z_num, double z_ll, double ct_ll, double gam_ll, double m_HNL, double m_p, double m_m);
  bool IsDalitzAllowed(double z_ll, double s, double d, double z_num);
  double MSqDM(double (&C)[6], double z_num, double z_ll, double ct_ll, double gam_ll, double m_HNL, double m_p, double m_m, double Pol);
  double MaxMSqDM(double m_HNL, int LeptonPDG, double Ue4, double Umu4, double Ut4, double m_p, double m_m, double Pol, double (&C)[6], bool Majorana, CLHEP::HepRandomEngine* fEngine);
  void RF4vecs(TLorentzVector &pnu, TLorentzVector &pm, TLorentzVector &pp, double z_num, double z_ll, double ct_ll, double gam_ll, double m_HNL, double m_p, double m_m, CLHEP::HepRandomEngine* fEngine);
  void AnisotropicThreeBodyDist(TLorentzVector pHNL,TLorentzVector &pnu, TLorentzVector &pm, TLorentzVector &pp, double m_HNL, double Ue4, double Umu4, double Ut4, int LeptonPDG, bool Majorana, bool AntiHNL, double Pol, CLHEP::HepRandomEngine* fEngine);
  void RotateToHNLPolFrame(TLorentzVector pHNL, TLorentzVector &pnu, TLorentzVector &pm, TLorentzVector &pp);

  double GetRandom(CLHEP::HepRandomEngine* fEngine);

} // namespace AnThreeBD
} // namespace ldm
} // namespace evgen

#endif
