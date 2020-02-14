////////////////////////////////////////////////////////////////////////
/// \class PMNSOpt
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework.
///
/// Two optimizations are relevant:
/// The construction of the Hamiltonian follows DocDB-XXXX (to be posted)
/// The eigensystem determination is based on the following reference:
///
///......................................................................
///
/// Int. J. Mod. Phys. C       VOLUME 19, NUMBER 03            MARCH 2008
///
///     Efficient numerical diagonalization of hermitian 3x3 matrices
///
///                            Joachim Kopp
///                  Max–Planck–Institut für Kernphysik
///             Postfach 10 39 80, 69029 Heidelberg, Germany
///                    (Received 19 October 2007)
///
///                                523
///......................................................................
///
/// The code structure follows the implementation written by M. Messier
/// in the PMNS class.
///
/// \version $Id: PMNSOpt.h,v 1.1 2013/01/19 16:09:57 jcoelho Exp $
///
/// @author joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////
#ifndef PMNSOPT_H
#define PMNSOPT_H
#include <list>
#include <complex>

namespace osc {
  // Some useful complex numbers
  static std::complex<double> zero(0.0,0.0);
  static std::complex<double> one (1.0,0.0);

  // Unit conversion constants
  static const double kKm2eV  = 5.06773103202e+09; ///< km to eV^-1
  static const double kK2     = 4.62711492217e-09; ///< mole/GeV^2/cm^3 to eV
  static const double kGeV2eV = 1.0e+09;           ///< GeV to eV

  //G_F in units of GeV^-2
  static const double kGf     = 1.166371e-5;

  /// Optimized version of \ref PMNS
  class PMNSOpt {
  public:
    PMNSOpt();
    virtual ~PMNSOpt();

    /// Set the parameters of the PMNS matrix
    /// @param th12    - The angle theta_12 in radians
    /// @param th23    - The angle theta_23 in radians
    /// @param th13    - The angle theta_13 in radians
    /// @param deltacp - The CPV phase delta_CP in radians
    virtual void SetMix(double th12, double th23, double th13, double deltacp);

    /// Set the mass-splittings
    /// @param dm21 - m2^2-m1^2 in eV^2
    /// @param dm32 - m3^2-m2^2 in eV^2
    virtual void SetDeltaMsqrs(double dm21, double dm32);

    /// Propagate a neutrino through a slab of matter
    /// @param L - length of slab (km)
    /// @param E - neutrino energy in GeV
    /// @param Ne - electron number density of matter in mole/cm^3
    /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
    virtual void PropMatter(double L, double E, double Ne, int anti=1);
    virtual void PropMatter(const std::list<double>& L,
                    double                   E,
                    const std::list<double>& Ne,
                    int anti);

    /// Propagate a neutrino through vacuum
    /// @param L - length of slab (km)
    /// @param E - neutrino energy in GeV
    /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
    virtual void PropVacuum(double L, double E, int anti=1);

    /// Return the probability of final state in flavour flv
    /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)
    virtual double P(int flv) const;

    /// Erase memory of neutrino propagate and reset neutrino
    /// to pure flavour flv. Preserves values of mixing and mass-splittings
    /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)
    virtual void ResetToFlavour(int flv=1);

  protected:
    // A shorthand...
    typedef std::complex<double> complex;

    /// Build H*lv, where H is the Hamiltonian in vacuum on flavour basis
    /// and lv is the oscillation length
    virtual void BuildHlv();

    /// Solve the full Hamiltonian for eigenvectors and eigenvalues
    /// @param E - neutrino energy in GeV
    /// @param Ne - electron number density of matter in mole/cm^3
    /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
    virtual void SolveHam(double E, double Ne, int anti);

    /// Set the eigensystem to the analytic solution of the vacuum Hamiltonian
    /// @param E - neutrino energy in GeV
    /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
    virtual void SetVacuumEigensystem(double E, int anti);

    double  fDm21;          ///< m^2_2 - m^2_1 in vacuum
    double  fDm31;          ///< m^2_3 - m^2_1 in vacuum
    double  fTheta12;       ///< theta12 mixing angle
    double  fTheta23;       ///< theta23 mixing angle
    double  fTheta13;       ///< theta13 mixing angle
    double  fDeltaCP;       ///< CP violating phase
    complex fHlv[3][3];     ///< dimensionless matrix H*lv
    complex fEvec[3][3];    ///< Eigenvectors of the Hamiltonian
    double  fEval[3];       ///< Eigenvalues of the Hamiltonian
    complex fNuState[3];    ///< The neutrino current state
    double  fCachedNe;      ///< Cached electron density
    double  fCachedE;       ///< Cached neutrino energy
    int     fCachedAnti;    ///< Cached nu/nubar selector
    bool    fBuiltHlv;      ///< Tag to avoid rebuilding Hlv
  };
}
#endif
////////////////////////////////////////////////////////////////////////
