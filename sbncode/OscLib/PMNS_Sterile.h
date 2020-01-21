////////////////////////////////////////////////////////////////////////
/// \class PMNS_Sterile
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        n-neutrino framework.
///
/// \version $Id: PMNS_Sterile.h,v 1.3 2014/04/17 21:01:58 jcoelho Exp $
///
/// @author joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////
#ifndef PMNS_STERILE_H
#define PMNS_STERILE_H
#include <list>
#include <complex>
#include <vector>

namespace osc {
  class PMNS_Sterile {
    public:

    PMNS_Sterile(int NumNus);
    virtual ~PMNS_Sterile();
    
    /// Set the parameters of the PMNS matrix
    /// @param th    - The angle theta_ij in radians
    virtual void SetAngle(int i, int j, double th);
    virtual void SetDelta(int i, int j, double delta);
    
    /// Set the mass-splittings
    /// @param dmi1 - mi^2-m1^2 in eV^2
    virtual void SetDm(int i, double dm);

    /// Set standard 3-flavor parameters
    virtual void SetStdPars();
    
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

    /// Return the probability of final state in flavour flv
    /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)
    virtual double P(int flv) const;
    
    /// Erase memory of neutrino propagate and reset neutrino
    /// to pure flavour flv. Preserves values of mixing and mass-splittings
    /// @param flv - final flavor (0,1,2) = (nue,numu,nutau)
    virtual void ResetToFlavour(int flv=1);
    
    /// Getters
    virtual int    GetNFlavors() const { return fNumNus; }
    virtual double GetDm(int i) const { return fDm[i-1]; }
    virtual double GetAngle(int i, int j) const { return fTheta[i-1][j-1]; }
    virtual double GetDelta(int i, int j) const { return fDelta[i-1][j-1]; }
    
  protected:
    
    // A shorthand...
    typedef std::complex<double> complex;

    virtual void InitializeVectors();

    /// Rotate the Hamiltonian by theta_ij and delta_ij
    virtual void RotateH(int i,int j);

    /// Build Hms = H*2E, where H is the Hamiltonian in vacuum on flavour basis   
    /// and E is the neutrino energy. This is effectively the matrix of masses squared.
    virtual void BuildHms();

    /// Solve the full Hamiltonian for eigenvectors and eigenvalues
    /// @param E - neutrino energy in GeV
    /// @param Ne - electron number density of matter in mole/cm^3
    /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
    virtual void SolveHam(double E, double Ne, int anti);

    int fNumNus;
    
    std::vector<double>                 fDm;      ///< m^2_i - m^2_1 in vacuum
    std::vector< std::vector<double> >  fTheta;   ///< theta[i][j] mixing angle
    std::vector< std::vector<double> >  fDelta;   ///< delta[i][j] CP violating phase

    std::vector<complex>                fNuState; ///< The neutrino current state
    std::vector< std::vector<complex> > fHms;     ///< matrix H*2E in eV^2

    double  fCachedNe;      ///< Cached electron density
    double  fCachedE;       ///< Cached neutrino energy
    int     fCachedAnti;    ///< Cached nu/nubar selector
    bool    fBuiltHms;      ///< Tag to avoid rebuilding Hms

    struct Priv;
    Priv* d;

  };
}
#endif
////////////////////////////////////////////////////////////////////////
