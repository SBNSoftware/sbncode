#include "HNLDecayDalitz.h"


// STUB
double evgen::ldm::HNLNuDiLepLNVDalitz(TLorentzVector K, TLorentzVector LA, TLorentzVector LB, TLorentzVector LBD, TLorentzVector NU) {

  auto pcross = [K, LA](TLorentzVector F) {
    return 2. * LA.Dot(K) * F.Dot(K) - LA.Dot(F) * K.Dot(K);
  };

  double A1 = 0.;
  double A2 = 0.;
  double A3 = 0.;

  return A1*2*NU.Dot(LB)*pcross(LBD) + A2*2*NU.Dot(LBD)*pcross(LB) + A3*LB.Dot(LB)*pcross(NU); 
}

double evgen::ldm::twobody_momentum(double parent_mass, double childA_mass, double childB_mass) {
  if (parent_mass < childA_mass + childB_mass) return -1.;

  return sqrt(parent_mass * parent_mass * parent_mass * parent_mass 
    -2 * parent_mass * parent_mass * childA_mass * childA_mass
    -2 * parent_mass * parent_mass * childB_mass * childB_mass
       + childA_mass * childA_mass * childA_mass * childA_mass 
       + childB_mass * childB_mass * childB_mass * childB_mass 
    -2 * childA_mass * childA_mass * childB_mass * childB_mass) / ( 2 * parent_mass );

}

double evgen::ldm::HNLLepPiLNCDalitz(TLorentzVector K, TLorentzVector LA, TLorentzVector N, TLorentzVector PI, TLorentzVector LB) {
  return 8 * K.Dot(LA) * PI.Dot(LB) * K.Dot(N) * PI.Dot(N) \
        -4 * K.Dot(LA) * PI.Dot(LB) * K.Dot(PI) * N.Dot(N) \
        -4 * K.Dot(K) * PI.Dot(LB) * LA.Dot(N) * PI.Dot(N) \
        +2 * K.Dot(K) * PI.Dot(LB) * LA.Dot(PI) * N.Dot(N) \
        -4 * K.Dot(LA) * PI.Dot(PI) * K.Dot(N) * LB.Dot(N) \
        +2 * K.Dot(LA) * PI.Dot(PI) * K.Dot(LB) * N.Dot(N) \
        +2 * K.Dot(K) * PI.Dot(PI) * LA.Dot(N) * LB.Dot(N) \
        -1 * K.Dot(K) * PI.Dot(PI) * LA.Dot(LB) * N.Dot(N);
}

double evgen::ldm::HNLLepPiLNVDalitz(TLorentzVector K, TLorentzVector LA, TLorentzVector N, TLorentzVector PI, TLorentzVector LB) {
  return N.Dot(N) * (
    4 * K.Dot(LA) * K.Dot(PI) * LB.Dot(PI) \
   -2 * K.Dot(K) * LA.Dot(PI) * LB.Dot(PI) \
   -2 * K.Dot(LA) * K.Dot(LB) * PI.Dot(PI) \
   +1 * K.Dot(K) * LA.Dot(LB) * PI.Dot(PI)
  );
}

// Work in the X frame
//              >
//          lB -
//            -
//     lA  N -   K
//    <-----o<-----
//         -
//        - pi
//       -
//      <
//
double evgen::ldm::HNLLepPiDalitzMax(double mK, double mA, double mN, double mP, double mB) {
  double pNA = evgen::ldm::twobody_momentum(mK, mA, mN);
  double pPB = evgen::ldm::twobody_momentum(mN, mP, mB);

  TLorentzVector K(TVector3(0, 0, pNA * mK / mN), mK * sqrt(1 + pNA * pNA / (mN*mN)));
  TLorentzVector LA(TVector3(0, 0, pNA * mK / mN), sqrt(mA*mA + pNA * pNA * mK * mK / (mN * mN)));
  TLorentzVector N(TVector3(0, 0, 0), mN);
  // In the LNC case, if mN > mA, then the direction of lB wants to be aligned with lA (+z)
  // (because of angular momentum conservation)
  // Flipping the mass ordering or LNC/LNV flips the sign
  double sign = mN > mA ? 1 : -1;
  TLorentzVector LB(TVector3(0, 0, sign * pPB), sqrt(mB * mB + pPB * pPB));
  TLorentzVector PI(TVector3(0, 0, -sign * pPB), sqrt(mP * mP + pPB * pPB));

  return HNLLepPiLNCDalitz(K, LA, N, PI, LB);
}
