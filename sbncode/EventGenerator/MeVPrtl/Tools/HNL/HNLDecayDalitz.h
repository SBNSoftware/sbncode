#ifndef HNLDecayDalitz_h
#define HNLDecayDalitz_h

#include "TLorentzVector.h"

namespace evgen {
namespace ldm {

  double HNLNuDiLepLNVDalitz(TLorentzVector K, TLorentzVector LA, TLorentzVector LB, TLorentzVector LBD, TLorentzVector NU); // TODO: implement

  double HNLLepPiLNVDalitz(TLorentzVector K, TLorentzVector LA, TLorentzVector N, TLorentzVector PI, TLorentzVector LB);
  double HNLLepPiLNCDalitz(TLorentzVector K, TLorentzVector LA, TLorentzVector N, TLorentzVector PI, TLorentzVector LB);

  double HNLLepPiDalitzMax(double mK, double mA, double mN, double mP, double mB);
  double twobody_momentum(double parent_mass, double childA_mass, double childB_mass);

} // namespace ldm
} // namespace evgen

#endif

