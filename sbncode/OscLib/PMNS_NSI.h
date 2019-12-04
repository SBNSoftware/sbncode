////////////////////////////////////////////////////////////////////////
/// \class PMNS_NSI
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with NSI. 
///
/// This class inherits from the PMNSOpt class
///
/// \version $Id: PMNS_NSI.h,v 1.2 2013/04/03 19:59:31 jcoelho Exp $
///
/// @author joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////
#ifndef PMNS_NSI_H
#define PMNS_NSI_H
#include <list>
#include <complex>

#include "OscLib/PMNSOpt.h"

namespace osc {
  class PMNS_NSI : public PMNSOpt {
  public:
    PMNS_NSI();
    virtual ~PMNS_NSI();
    
    void SetNSI(double eps_ee,      double eps_emu,      double eps_etau,
                double eps_mumu,    double eps_mutau,    double eps_tautau,
                double delta_emu=0, double delta_etau=0, double delta_mutau=0);

  protected:
    /// Solve the full Hamiltonian for eigenvectors and eigenvalues
    /// @param E - neutrino energy in GeV
    /// @param Ne - electron number density of matter in mole/cm^3
    /// @param anti - +1 = neutrino case, -1 = anti-neutrino case
    virtual void SolveHam(double E, double Ne, int anti);
    
    double  fEps_ee;        ///< NSI parameter ee
    double  fEps_mumu;      ///< NSI parameter mumu
    double  fEps_tautau;    ///< NSI parameter tautau
    complex fEps_emu;       ///< NSI parameter emu
    complex fEps_etau;      ///< NSI parameter etau
    complex fEps_mutau;     ///< NSI parameter mutau
    bool    fResetNSI;      ///< True when NSI parameters are changed
  };
}
#endif
////////////////////////////////////////////////////////////////////////
