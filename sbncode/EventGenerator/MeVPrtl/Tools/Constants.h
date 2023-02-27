#ifndef IMeVPrtlConstants_h
#define IMeVPrtlConstants_h

#include "TVector3.h"

#include "fhiclcpp/ParameterSet.h"

namespace evgen {
namespace ldm {

// NOTE: all values here use units of GeV/ns/cm

class Constants {
private:
  // Singleton
  static Constants &InstanceMut() {
    static Constants instance;

    return instance;
  }

public:
  // Masses
  double elec_mass;
  double muon_mass;
  double piplus_mass;
  double pizero_mass;
  double kplus_mass;
  double klong_mass;
  double tquark_mass;
  double tau_mass;
  double eta_mass;
  double rho_mass;
  double etap_mass;

  // Couplings
  double fine_structure_constant;
  double Gfermi;
  double higgs_vev;
  double sin2thetaW;
  double gL;
  double gR;
  double fpion;
  double feta;
  double fetap;
  double grho;

  // unit conversion
  double hbar;
  double c_cm_per_ns;

  // kaon lifetimes
  double kplus_lifetime;
  double klong_lifetime;

  // other lifetimes
  double tau_lifetime;

  // and widths
  double rho_width;

  // kaon decay branching ratios
  double kaonp_mup_numu;
  double kaonp_ep_nue;

  // CKM matrix
  double abs_Vud_squared;
  double abs_VtsVtd_squared;
  double rel_VtsVtd_squared;


  Constants();

  // Delete bad functions
  Constants(Constants const&) = delete;
  void operator=(Constants const&)  = delete;

  // For configuration
  static void Configure(const fhicl::ParameterSet &p);

  // Public access
  static const Constants &Instance() { return InstanceMut(); }
};

// Useful computations
double twobody_momentum(double parent_mass, double childA_mass, double childB_mass);
int calcPrtlRayWgt(double rest_frame_p, double M, TVector3 dir, TVector3 boost, double rand,
                     double& lab_frame_p_out, double& costh_rest_out, double& wgt);
double forwardPrtlEnergy(double parentM, double secM, double prtlM, double parentE);
double PDG2Mass(int pdg);

// Minimum possible cos theta for a given decay
double minKinematicCosTheta(double parentM, double secM, double prtlM, double parentE);

  }
}
#endif
