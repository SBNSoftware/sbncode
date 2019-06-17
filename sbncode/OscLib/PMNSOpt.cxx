////////////////////////////////////////////////////////////////////////
// $Id: PMNSOpt.cxx,v 1.1 2013/01/19 16:09:57 jcoelho Exp $
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework. Two optimizations are relevant:
// The construction of the Hamiltonian follows DocDB-XXXX (to be posted)
// The eigensystem determination is based on the following reference:
//
//......................................................................
//
// Int. J. Mod. Phys. C       VOLUME 19, NUMBER 03            MARCH 2008
//
//     Efficient numerical diagonalization of hermitian 3x3 matrices
//
//                            Joachim Kopp
//                  Max–Planck–Institut für Kernphysik
//             Postfach 10 39 80, 69029 Heidelberg, Germany
//                    (Received 19 October 2007)
//
//                                523
//......................................................................
//
// Throughout I have taken:
//   - L to be the neutrino flight distance in km
//   - E to be the neutrino energy in GeV
//   - Dmij to be the differences between the mass-squares in eV^2
//   - Ne to be the electron number density in mole/cm^3
//   - theta12,theta23,theta13,deltaCP to be in radians
//
// The code structure follows the implementation written by M. Messier
// in the PMNS class.
//
// joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////

#include "OscLib/PMNSOpt.h"

// Just pull in all the cxx files from MatrixDecomp. This way we don't have
// another library to worry about building and linking everywhere.
#include "OscLib/MatrixDecomp/zhetrd3.cxx"
#include "OscLib/MatrixDecomp/zheevc3.cxx"
#include "OscLib/MatrixDecomp/zheevh3.cxx"
#include "OscLib/MatrixDecomp/zheevq3.cxx"

#include <cstdlib>
#include <cassert>
#include <math.h>

using namespace osc;

//......................................................................

PMNSOpt::PMNSOpt()
{
  this->SetMix(0.,0.,0.,0.);
  this->SetDeltaMsqrs(0.,0.);
  this->ResetToFlavour(1);
  fCachedNe = 0.0;
  fCachedE =  1.0;
  fCachedAnti = 1;
  fBuiltHlv = false;
}

PMNSOpt::~PMNSOpt(){
}

//......................................................................

void PMNSOpt::SetMix(double th12, double th23, double th13, double deltacp)
{

  fTheta12 = th12;
  fTheta23 = th23;
  fTheta13 = th13;
  fDeltaCP = deltacp;

  fBuiltHlv = false;

}

//......................................................................
///
/// Set the mass-splittings. These are m_2^2-m_1^2, m_3^2-m_2^2
/// and m_3^2-m_1^2 in eV^2
///
void PMNSOpt::SetDeltaMsqrs(double dm21, double dm32)
{

  fDm21 = dm21;
  fDm31 = dm32 + dm21;

  if(fDm31==0) fDm31 = 1.0e-12;

  fBuiltHlv = false;

}

//......................................................................
///
/// Build H*lv, where H is the Hamiltonian in vacuum on flavour basis
/// and lv is the oscillation length
///
/// This is a dimentionless hermitian matrix, so only the
/// upper triangular part needs to be filled
///
/// The construction of the Hamiltonian avoids computing terms that
/// are simply zero. This has a big impact in the computation time.
/// This construction is described in DocDB-XXXX (to be posted)
///
void PMNSOpt::BuildHlv()
{

  // Check if anything changed
  if(fBuiltHlv) return;

  // Create temp variables
  double sij, cij, h00, h11, h01;
  complex expCP, h02, h12;
  
  // Hamiltonian in mass base. Only one entry is variable.
  h11 = fDm21 / fDm31;
  
  // Rotate over theta12
  sij = sin(fTheta12);
  cij = cos(fTheta12);
 
  // There are 3 non-zero entries after rephasing so that h22 = 0
  h00 = h11 * sij * sij - 1;
  h01 = h11 * sij * cij;
  h11 = h11 * cij * cij - 1;

  // Rotate over theta13 with deltaCP
  sij = sin(fTheta13);
  cij = cos(fTheta13);
  expCP = complex(cos(fDeltaCP), -sin(fDeltaCP));
  
  // There are 5 non-zero entries after rephasing so that h22 = 0
  h02 = (-h00 * sij * cij) * expCP;
  h12 = (-h01 * sij) * expCP;
  h11 -= h00 * sij * sij;                                         
  h00 *= cij * cij  -  sij * sij;
  h01 *= cij;

  // Finally, rotate over theta23
  sij = sin(fTheta23);
  cij = cos(fTheta23);

  // Fill the Hamiltonian rephased so that h22 = -h11
  fHlv[0][0] = h00 - 0.5 * h11;
  fHlv[1][1] = 0.5 * h11 * (cij * cij - sij * sij)  +  2 * real(h12) * cij * sij;
  fHlv[2][2] = -fHlv[1][1];

  fHlv[0][1] = h02 * sij  +  h01 * cij;
  fHlv[0][2] = h02 * cij  -  h01 * sij;
  fHlv[1][2] = h12 - (h11 * cij + 2 * real(h12) * sij) * sij;

  // Tag as built
  fBuiltHlv = true;

}

//......................................................................
///
/// Solve the full Hamiltonian for eigenvectors and eigenvalues
/// This is using a method from the GLoBES software available at
/// http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
/// We should cite them accordingly
///
void PMNSOpt::SolveHam(double E, double Ne, int anti)
{

  // Check if anything has changed before recalculating
  if(Ne!=fCachedNe || E!=fCachedE || anti!=fCachedAnti || !fBuiltHlv ){
    fCachedNe = Ne;
    fCachedE = E;
    fCachedAnti = anti;
    this->BuildHlv();
  }
  else return;

  double lv = 2 * kGeV2eV*E / fDm31;  // Osc. length in eV^-1
  double kr2GNe = kK2*M_SQRT2*kGf*Ne; // Matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV
  complex A[3][3];
  for(int i=0;i<3;i++){
    A[i][i] = fHlv[i][i]/lv;
    for(int j=i+1;j<3;j++){
      if(anti>0) A[i][j] = fHlv[i][j]/lv;
      else       A[i][j] = conj(fHlv[i][j])/lv;
    }
  }
  if(anti>0) A[0][0] += kr2GNe;
  else       A[0][0] -= kr2GNe;

  // Solve Hamiltonian for eigensystem using the GLoBES method
  zheevh3(A,fEvec,fEval);

}

///.....................................................................
///
/// Propagate the current neutrino state over a distance L in km
/// with an energy E in GeV through constant matter of density
/// Ne in mole/cm^3.
/// @param anti - +1 = neutrino case, -1 = anti-neutrino case
///
void PMNSOpt::PropMatter(double L, double E, double Ne, int anti)
{

  // Solve Hamiltonian
  this->SolveHam(E, Ne, anti);

  // Store coefficients of propagation eigenstates
  complex nuComp[3];

  for(int i=0;i<3;i++){
    nuComp[i] = 0;
    for(int j=0;j<3;j++){
      nuComp[i] += fNuState[j] * conj(fEvec[j][i]);
    }
  }

  for(int i=0;i<3;i++)fNuState[i] = 0;

  // Propagate neutrino state
  for(int j=0;j<3;j++){
    double s, c;

#ifdef DARWINBUILD
    s = std::sin(-fEval[j] * kKm2eV * L);
    c = std::cos(-fEval[j] * kKm2eV * L);
#else
    sincos(-fEval[j] * kKm2eV*L, &s, &c);
#endif
    
    complex jPart = complex(c, s) * nuComp[j];

    for(int i=0;i<3;i++){
      fNuState[i] += jPart * fEvec[i][j];
    }
  }

}

//......................................................................
///
/// Do several layers in a row. L and Ne must have the same length
///
void PMNSOpt::PropMatter(const std::list<double>& L,
                           double                   E,
                           const std::list<double>& Ne,
                           int anti)
{
  if (L.size()!=Ne.size()) abort();
  std::list<double>::const_iterator Li  (L.begin());
  std::list<double>::const_iterator Lend(L.end());
  std::list<double>::const_iterator Ni  (Ne.begin());
  for (; Li!=Lend; ++Li, ++Ni) {
    // For very low densities, use vacumm
    static const double kRhoCutoff = 1.0E-6; // Ne in moles/cm^3
    if (*Ni<kRhoCutoff) this->PropVacuum(*Li, E, anti);
    else                this->PropMatter(*Li, E, *Ni, anti);
  }
}


///.....................................................................
///
/// We know the vacuum eigensystem, so just write it explicitly
/// The eigenvalues depend on energy, so E needs to be provided in GeV
/// @param anti - +1 = neutrino case, -1 = anti-neutrino case
///
void PMNSOpt::SetVacuumEigensystem(double E, int anti)
{

  double  s12, s23, s13, c12, c23, c13;
  complex expidelta(cos(fDeltaCP), anti * sin(fDeltaCP));

  s12 = sin(fTheta12);  s23 = sin(fTheta23);  s13 = sin(fTheta13);
  c12 = cos(fTheta12);  c23 = cos(fTheta23);  c13 = cos(fTheta13);

  fEvec[0][0] =  c12*c13;
  fEvec[0][1] =  s12*c13;
  fEvec[0][2] =  s13*conj(expidelta);

  fEvec[1][0] = -s12*c23-c12*s23*s13*expidelta;
  fEvec[1][1] =  c12*c23-s12*s23*s13*expidelta;
  fEvec[1][2] =  s23*c13;

  fEvec[2][0] =  s12*s23-c12*c23*s13*expidelta;
  fEvec[2][1] = -c12*s23-s12*c23*s13*expidelta;
  fEvec[2][2] =  c23*c13;

  fEval[0] = 0;
  fEval[1] = fDm21 / (2 * kGeV2eV*E);
  fEval[2] = fDm31 / (2 * kGeV2eV*E);

}

///.....................................................................
///
/// Propagate the current neutrino state over a distance L in km
/// with an energy E in GeV through vacuum
/// @param anti - +1 = neutrino case, -1 = anti-neutrino case
///
void PMNSOpt::PropVacuum(double L, double E, int anti)
{

  this->SetVacuumEigensystem(E, anti);

  complex nuComp[3];

  for(int i=0;i<3;i++){
    nuComp[i] = 0;
    for(int j=0;j<3;j++){
      nuComp[i] += fNuState[j] * conj(fEvec[j][i]);
    }
  }

  for(int i=0;i<3;i++){
    fNuState[i] = 0;
    for(int j=0;j<3;j++){
      complex iEval(0.0,fEval[j]);
      fNuState[i] +=  exp(-iEval * kKm2eV*L) * nuComp[j] * fEvec[i][j];
    }
  }

}

//......................................................................
///
/// Reset the neutrino state back to a pure flavour where
/// it starts
///
void PMNSOpt::ResetToFlavour(int flv)
{
  int i;
  for (i=0; i<3; ++i){
    if (i==flv) fNuState[i] = one;
    else        fNuState[i] = zero;
  }
}

//......................................................................
///
/// Compute oscillation probability of flavour flv
///
/// 0 = nue, 1 = numu, 2 = nutau
///
double PMNSOpt::P(int flv) const
{
  assert(flv>=0 && flv<3);
  return norm(fNuState[flv]);
}

////////////////////////////////////////////////////////////////////////
