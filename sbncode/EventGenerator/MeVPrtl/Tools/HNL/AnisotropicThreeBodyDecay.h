#ifndef AnisotropicThreeBodyDecay_h
#define AnisotropicThreeBodyDecay_h

#include "TLorentzVector.h"


//Implementation of HNL Three Body Decays Anisotropies
//Valid as long as the HNL is decaying into a neutrino and identical final-state charged leptons
//See arXiv:2104.05719 for more details @LuisPelegrina


namespace evgen {
namespace ldm {
namespace AnThreeBD {
  
  double GL(bool CC);
  double GR();
  void CMgLgR(double (&CM)[6],double Ue4,double Umu4,double Ut4, int LeptonPDG);
  void CDAgLgR(double (&CD)[6],double Ue4,double Umu4,double Ut4, int LeptonPDG);
  void CDNgLgR(double (&CD)[6],double Ue4,double Umu4,double Ut4, int LeptonPDG);
  void KinDep(double (&K)[6],double Pol,double z_num,double z_ll, double ct_ll,double gam_ll,double m_HNL,double m_p,double m_m);
  bool IsDalitzAllowed(double z_ll, double s, double d,double z_num);
  double MSqDM(double (&C)[6],double z_num, double z_ll,  double ct_ll, double gam_ll,double m_HNL,double m_p,double m_m,double Pol);
  double MaxMSqDM(double m_HNL,int  LeptonPDG,double Ue4,double Umu4,double Ut4,double m_p,double m_m,double Pol,double (&C)[6],bool Majorana);
  void RF4vecs(TLorentzVector &pnu,TLorentzVector &pm,TLorentzVector &pp,double z_num,double z_ll, double ct_ll,double gam_ll,double m_HNL,double m_p,double m_m);
  //void AnisotropicThreeBodyDist((double (&pnu)[4],(double (&pm)[4],(double (&pp)[4],double m_HNL,double Ue4,double Umu4,double Ut4,int LeptonPDG,bool Majorana,double Pol);
  void AnisotropicThreeBodyDist(TLorentzVector &pnu,TLorentzVector &pm,TLorentzVector &pp,double m_HNL,double Ue4,double Umu4,double Ut4,int LeptonPDG,bool Majorana,bool AntiHNL,double Pol);
  
  double PolHNL(double m_HNL, double m_l, double m_M);
    
} // namespace AnThreeBD
} // namespace ldm
} // namespace evgen

#endif
