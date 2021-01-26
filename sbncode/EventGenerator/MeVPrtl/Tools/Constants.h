#ifndef IMeVPrtlConstants_h
#define IMeVPrtlConstants_h

#include "TVector3.h"

namespace evgen {
namespace ldm {

// NOTE: all values here use units of GeV/ns/cm

// masses
static const double elec_mass = 0.000511; // GeV
static const double muon_mass = 0.105658; // GeV
static const double piplus_mass = 0.139570; // GeV
static const double pizero_mass = 0.134977; // GeV
static const double kaonp_mass = 0.493677; // GeV
static const double pionp_mass = 0.139570; // GeV
static const double kplus_mass = 0.493677; // GeV
static const double klong_mass = 0.497611; // GeV
static const double pplus_mass = 0.139570; // GeV
static const double pzero_mass = 0.134977; // GeV
static const double tquark_mass = 172.76; // GeV

// Couplings
static const double Gfermi = 1.16638e-5; // 1/GeV^2
static const double higgs_vev = 246.22; // GeV (=1/sqrt(Gfermi*sqrt(2)))
static const double sin2thetaW = 0.2229; // electroweak mixing angle
static const double gL = -0.5 + sin2thetaW;
static const double gR = sin2thetaW;
static const double fpion = 0.1302; // Pion decay constant [GeV]

// unit conversion
static const double hbar = 6.582119569e-16; // GeV*ns or eV*s
static const double c_cm_per_ns = 29.9792; // cm / ns

// kaon lifetimes
static const double kplus_lifetime = 1.238e1; // ns
static const double klong_lifetime = 5.116e1; // ns

// Kaon decay branching ratios
static const double kaonp_mup_numu = 0.6356; // From PDG: https://pdg.lbl.gov/2018/listings/rpp2018-list-K-plus-minus.pdf 
static const double kaonp_ep_nue = 1.582e-5; // From PDG

// CKM matrix
static const double abs_VtsVtd_squared = 1.0185e-07;
static const double rel_VtsVtd_squared = 1.0185e-07;
static const double abs_Vud_squared = 0.97420 * 0.97420;
  
// Useful computations
double twobody_momentum(double parent_mass, double childA_mass, double childB_mass);
int calcPrtlRayWgt(double rest_frame_p, double M, TVector3 dir, TVector3 boost, double rand,
                     double& lab_frame_p_out, double& costh_rest_out, double& wgt);
double forwardPrtlEnergy(double parentM, double secM, double prtlM, double parentE);
double secPDG2Mass(int pdg);


}

}
#endif
