#ifndef IMeVPrtlConstants_h
#define IMeVPrtlConstants_h

#include "TVector3.h"

namespace evgen {
namespace ldm {

// NOTE: all values here use units of GeV/ns/cm

// Masses

// P.A. Zylaet al.(Particle Data Group), Prog. Theor. Exp. Phys.2020, 083C01 (2020)

static const double elec_mass = 0.0005109989461; // GeV https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
static const double muon_mass = 0.1056583745; // GeV https://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf
static const double piplus_mass = 0.13957039; // GeV https://pdg.lbl.gov/2020/tables/rpp2020-tab-mesons-light.pdf
static const double pizero_mass = 0.1349768; // GeV https://pdg.lbl.gov/2020/tables/rpp2020-tab-mesons-light.pdf
static const double pionp_mass = piplus_mass; // GeV https://pdg.lbl.gov/2020/tables/rpp2020-tab-mesons-light.pdf
static const double pplus_mass = piplus_mass; // GeV
static const double pzero_mass = pizero_mass; // GeV
static const double kaonp_mass = 0.493677; // GeV https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf 
static const double kplus_mass = kaonp_mass; // GeV https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf 
static const double klong_mass = 0.497611; // GeV https://pdg.lbl.gov/2020/listings/rpp2020-list-K-zero.pdf
static const double tquark_mass = 172.76; // GeV https://pdg.lbl.gov/2020/tables/rpp2020-sum-quarks.pdf (direct measurements)

// Couplings
static const double Gfermi = 1.166379e-5; // 1/GeV^2 https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf 
static const double higgs_vev = 1. / sqrt(sqrt(2)*Gfermi); // GeV (246.22)
static const double sin2thetaW = 0.2312; // electroweak mixing angle https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
static const double gL = -0.5 + sin2thetaW;
static const double gR = sin2thetaW;
static const double fpion = 0.1302; // Pion decay constant [GeV] https://pdg.lbl.gov/2020/reviews/rpp2020-rev-pseudoscalar-meson-decay-cons.pdf (FLAG 19 average)

// unit conversion
static const double hbar = 6.582119569e-16; // GeV*ns or eV*s https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
static const double c_cm_per_ns = 29.9792458; // cm / ns https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf

// kaon lifetimes
static const double kplus_lifetime = 1.238e1; // ns https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
static const double klong_lifetime = 5.116e1; // ns https://pdg.lbl.gov/2020/listings/rpp2020-list-K-zero-L.pdf (FIT)

// Kaon decay branching ratios
static const double kaonp_mup_numu = 0.6356; // From PDG: https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
static const double kaonp_ep_nue = 1.582e-5; // From PDG: https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf

// CKM matrix
static const double abs_Vud_squared = 0.97370 * 0.97370; // https://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf (12.7)

// Computed using the Wolfenstein parameterization, where:
// Vts = -A \lambda^2
// Vtd = A \lambda^3 (1 - \rho - I\eta)
//
// With, from https://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf:
// A = 0.790
// \lambda = 0.2265
// \rho = 0.141
// \eta = 0.357
static const double abs_VtsVtd_squared = 1.19777e-7;
static const double rel_VtsVtd_squared = 1.02136e-7;
// (OLD: 1.0185e-07)
  
// Useful computations
double twobody_momentum(double parent_mass, double childA_mass, double childB_mass);
int calcPrtlRayWgt(double rest_frame_p, double M, TVector3 dir, TVector3 boost, double rand,
                     double& lab_frame_p_out, double& costh_rest_out, double& wgt);
double forwardPrtlEnergy(double parentM, double secM, double prtlM, double parentE);
double secPDG2Mass(int pdg);

// Minimum possible cos theta for a given decay
double minKinematicCosTheta(double parentM, double secM, double prtlM, double parentE);

  }
}
#endif
