////////////////////////////////////////////////////////////////////////
// $Id: PMNS_Sterile.cxx,v 1.7 2014/04/17 21:34:38 jcoelho Exp $
//
// Implementation of oscillations of neutrinos in matter in a
// n-neutrino framework. 
//
//......................................................................
//
// Throughout I have taken:
//   - L to be the neutrino flight distance in km
//   - E to be the neutrino energy in GeV
//   - Dmij to be the differences between the mass-squares in eV^2
//   - Ne to be the electron number density in mole/cm^3
//   - theta_ij and delta_ij to be in radians
//
// joao.coelho@tufts.edu
////////////////////////////////////////////////////////////////////////

#include "OscLib/PMNS_Sterile.h"

#include <iostream>
#include <cassert>
#include <stdlib.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_cgsm.h>

namespace osc
{
  struct PMNS_Sterile::Priv
  {
    Priv(int numNu)
    {
      // Allocate memory for the GSL objects
      fEval = gsl_vector_alloc(numNu);
      fEvec = gsl_matrix_complex_alloc(numNu, numNu);
      H_GSL = gsl_matrix_complex_alloc(numNu, numNu);
      W_GSL = gsl_eigen_hermv_alloc(numNu);
    }

    Priv(const Priv& p)
    {
      const int numNu = p.fEval->size;
      fEval = gsl_vector_alloc(numNu);
      fEvec = gsl_matrix_complex_alloc(numNu, numNu);
      H_GSL = gsl_matrix_complex_alloc(numNu, numNu);
      W_GSL = gsl_eigen_hermv_alloc(numNu);

      gsl_vector_memcpy(fEval, p.fEval);
      gsl_matrix_complex_memcpy(fEvec, p.fEvec);
      gsl_matrix_complex_memcpy(H_GSL, p.H_GSL);
      // No way to copy. Hope leaving it empty is OK
    }

    ~Priv()
    {
      // Free memory from GSL objects
      gsl_matrix_complex_free(H_GSL);
      gsl_eigen_hermv_free(W_GSL);
      gsl_matrix_complex_free(fEvec);
      gsl_vector_free(fEval);
    }

    gsl_vector* fEval;
    gsl_matrix_complex* fEvec;
    gsl_matrix_complex* H_GSL;
    gsl_eigen_hermv_workspace* W_GSL;
  };
}


// Some useful complex numbers
static std::complex<double> zero(0.0,0.0);
static std::complex<double> one (1.0,0.0);

// Unit conversion constants
static const double kKm2eV  = 5.06773103202e+09; ///< km to eV^-1
static const double kK2     = 4.62711492217e-09; ///< mole/GeV^2/cm^3 to eV
static const double kGeV2eV = 1.0e+09;           ///< GeV to eV

//G_F in units of GeV^-2
static const double kGf     = 1.166371E-5;

using namespace std;

using namespace osc;

//......................................................................

PMNS_Sterile::PMNS_Sterile(int NumNus) 
  : d(new PMNS_Sterile::Priv(NumNus))
{

  fNumNus = NumNus;
  this->SetStdPars();
  this->ResetToFlavour(1);
  fCachedNe = 0.0;
  fCachedE =  1.0;
  fCachedAnti = 1;
}

PMNS_Sterile::~PMNS_Sterile()
{
  delete d;
}

void PMNS_Sterile::InitializeVectors()
{

  fDm    = vector<double>(fNumNus, 0);
  fTheta = vector< vector<double> >(fNumNus, vector<double>(fNumNus,0));
  fDelta = vector< vector<double> >(fNumNus, vector<double>(fNumNus,0));

  fNuState = vector<complex>(fNumNus, zero);
  fHms     = vector< vector<complex> >(fNumNus, vector<complex>(fNumNus,zero));

}

void PMNS_Sterile::SetStdPars()
{

  this->InitializeVectors();

  if(fNumNus>2){
    this->SetAngle(1,2,0.6);
    this->SetAngle(1,3,0.16);
    this->SetAngle(2,3,0.7);
    this->SetDm(2,7.6e-5);
    this->SetDm(3,2.4e-3);
  }
  else if(fNumNus==2){
    this->SetAngle(1,2,0.7);
    this->SetDm(2,2.4e-3);
  }
  
}

//......................................................................
///
/// Set mixing angles in radians. Requires i < j.
/// 
void PMNS_Sterile::SetAngle(int i, int j, double th) 
{

  if(i>j){
    cout << "First argument should be smaller than second argument" << endl;
    cout << "Setting reverse order (Theta" << j << i << "). " << endl;
    int temp = i;
    i = j;
    j = temp;
  }
  if(i<1 || i>fNumNus-1 || j<2 || j>fNumNus){
    cout << "Theta" << i << j << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }

  fTheta[i-1][j-1] = th;

  fBuiltHms = false;

}

//......................................................................
///
/// Set CP violating phases in radians. Requires i+1 < j.
/// 
void PMNS_Sterile::SetDelta(int i, int j, double delta) 
{

  if(i>j){
    cout << "First argument should be smaller than second argument" << endl;
    cout << "Setting reverse order (Delta" << j << i << "). " << endl;
    int temp = i;
    i = j;
    j = temp;
  }
  if(i<1 || i>fNumNus-1 || j<2 || j>fNumNus){
    cout << "Delta" << i << j << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }
  if(i+1==j){
    cout << "Rotation " << i << j << " is real. Doing nothing." << endl;
    return;
  }

  fDelta[i-1][j-1] = delta;

  fBuiltHms = false;

}

//......................................................................
///
/// Set the mass-splittings. These are m_j^2-m_1^2
///
void PMNS_Sterile::SetDm(int j, double dm) 
{

  if(j<2 || j>fNumNus){
    cout << "Dm" << j << "1 not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }

  fDm[j-1] = dm;

  fBuiltHms = false;

}

//......................................................................
///
/// Rotate the Hamiltonian by the angle theta_ij and phase delta_ij.
/// The rotations assume all off-diagonal elements with i > j are zero.
/// This is correct if the order of rotations is chosen appropriately
/// and it speeds up computation by skipping null terms
///
void PMNS_Sterile::RotateH(int i,int j){

  // Do nothing if angle is zero
  if(fTheta[i][j]==0) return;

  double fSinBuffer = sin(fTheta[i][j]);
  double fCosBuffer = cos(fTheta[i][j]);

  double  fHmsBufferD;
  complex fHmsBufferC;

  // With Delta
  if(i+1<j){
    complex fExpBuffer = complex(cos(fDelta[i][j]), -sin(fDelta[i][j]));

    // General case
    if(i>0){
      // Top columns
      for(int k=0; k<i; k++){
        fHmsBufferC = fHms[k][i];

        fHms[k][i] *= fCosBuffer;
        fHms[k][i] += fHms[k][j] * fSinBuffer * conj(fExpBuffer);

        fHms[k][j] *= fCosBuffer;
        fHms[k][j] -= fHmsBufferC * fSinBuffer * fExpBuffer;
      }

      // Middle row and column
      for(int k=i+1; k<j; k++){
        fHmsBufferC = fHms[k][j];
    
        fHms[k][j] *= fCosBuffer;
        fHms[k][j] -= conj(fHms[i][k]) * fSinBuffer * fExpBuffer;

        fHms[i][k] *= fCosBuffer;
        fHms[i][k] += fSinBuffer * fExpBuffer * conj(fHmsBufferC);
      }

      // Nodes ij
      fHmsBufferC = fHms[i][i];
      fHmsBufferD = real(fHms[j][j]);

      fHms[i][i] *= fCosBuffer * fCosBuffer;
      fHms[i][i] += 2 * fSinBuffer * fCosBuffer * real(fHms[i][j] * conj(fExpBuffer));
      fHms[i][i] += fSinBuffer * fHms[j][j] * fSinBuffer;

      fHms[j][j] *= fCosBuffer * fCosBuffer;
      fHms[j][j] += fSinBuffer * fHmsBufferC * fSinBuffer;
      fHms[j][j] -= 2 * fSinBuffer * fCosBuffer * real(fHms[i][j] * conj(fExpBuffer));
  
      fHms[i][j] -= 2 * fSinBuffer * real(fHms[i][j] * conj(fExpBuffer)) * fSinBuffer * fExpBuffer;
      fHms[i][j] += fSinBuffer * fCosBuffer * (fHmsBufferD - fHmsBufferC) * fExpBuffer;

    }
    // First rotation on j (No top columns)
    else{
      // Middle rows and columns
      for(int k=i+1; k<j; k++){
        fHms[k][j] = -conj(fHms[i][k]) * fSinBuffer * fExpBuffer;

        fHms[i][k] *= fCosBuffer;
      }

      // Nodes ij
      fHmsBufferD = real(fHms[i][i]);

      fHms[i][j] = fSinBuffer * fCosBuffer * (fHms[j][j] - fHmsBufferD) * fExpBuffer;

      fHms[i][i] *= fCosBuffer * fCosBuffer;
      fHms[i][i] += fSinBuffer * fHms[j][j] * fSinBuffer;

      fHms[j][j] *= fCosBuffer * fCosBuffer;
      fHms[j][j] += fSinBuffer * fHmsBufferD * fSinBuffer;
    }
  
  }
  // Without Delta (No middle rows or columns: j = i+1)
  else{
    // General case
    if(i>0){
      // Top columns
      for(int k=0; k<i; k++){
        fHmsBufferC = fHms[k][i];

        fHms[k][i] *= fCosBuffer;
        fHms[k][i] += fHms[k][j] * fSinBuffer;

        fHms[k][j] *= fCosBuffer;
        fHms[k][j] -= fHmsBufferC * fSinBuffer;
      }

      // Nodes ij
      fHmsBufferC = fHms[i][i];
      fHmsBufferD = real(fHms[j][j]);

      fHms[i][i] *= fCosBuffer * fCosBuffer;
      fHms[i][i] += 2 * fSinBuffer * fCosBuffer * real(fHms[i][j]);
      fHms[i][i] += fSinBuffer * fHms[j][j] * fSinBuffer;

      fHms[j][j] *= fCosBuffer * fCosBuffer;
      fHms[j][j] += fSinBuffer * fHmsBufferC * fSinBuffer;
      fHms[j][j] -= 2 * fSinBuffer * fCosBuffer * real(fHms[i][j]);
  
      fHms[i][j] -= 2 * fSinBuffer * real(fHms[i][j]) * fSinBuffer;
      fHms[i][j] += fSinBuffer * fCosBuffer * (fHmsBufferD - fHmsBufferC);

    }
    // First rotation (theta12)
    else{

      fHms[i][j] = fSinBuffer * fCosBuffer * fHms[j][j];

      fHms[i][i] = fSinBuffer * fHms[j][j] * fSinBuffer;

      fHms[j][j] *= fCosBuffer * fCosBuffer;
  
    }
  }

}

//......................................................................
///
/// Build Hms = H*2E, where H is the Hamiltonian in vacuum on flavour basis
/// and E is the neutrino energy in eV. Hms is effectively the matrix of 
/// masses squared.
///
/// This is a hermitian matrix, so only the
/// upper triangular part needs to be filled
///
/// The construction of the Hamiltonian avoids computing terms that
/// are simply zero. This has a big impact in the computation time.
/// This construction is described in DocDB-XXXX (to be posted)
///
void PMNS_Sterile::BuildHms() 
{

  // Check if anything changed
  if(fBuiltHms) return;

  for(int j=0; j<fNumNus; j++){
    // Set mass splitting
    fHms[j][j] = fDm[j];
    // Reset off-diagonal elements
    for(int i=0; i<j; i++){
      fHms[i][j] = 0;
    }
    // Rotate j neutrinos
    for(int i=0; i<j; i++){
      this->RotateH(i,j);
    }
  }
  
  // Tag as built
  fBuiltHms = true;

}

//......................................................................
///
/// Solve the full Hamiltonian for eigenvectors and eigenvalues
/// This uses GSL to solve the eigensystem.
///
/// The matter effect assumes # of neutrons = # of electrons.
///
void PMNS_Sterile::SolveHam(double E, double Ne, int anti)
{

  // Check if anything has changed before recalculating
  if(Ne!=fCachedNe || E!=fCachedE || anti!=fCachedAnti || !fBuiltHms ){
    fCachedNe = Ne;
    fCachedE = E;
    fCachedAnti = anti;
    this->BuildHms();
  }
  else return;

  double lv = 2 * kGeV2eV*E;          // 2E in eV 
  double kr2GNe = kK2*M_SQRT2*kGf*Ne; // Matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV

  for(size_t i=0;i<size_t(fNumNus);i++){
    complex buf = fHms[i][i]/lv;
    *gsl_matrix_complex_ptr(d->H_GSL, i, i) = gsl_complex_rect(real(buf),0);
    for(size_t j=i+1;j<size_t(fNumNus);j++){
      buf = fHms[i][j]/lv;
      if(anti>0) *gsl_matrix_complex_ptr(d->H_GSL, j, i) = gsl_complex_rect(real(buf),-imag(buf));
      else       *gsl_matrix_complex_ptr(d->H_GSL, j, i) = gsl_complex_rect(real(buf), imag(buf));
    }
    if(i>2){
      // Subtract NC coherent forward scattering from sterile neutrinos. See arXiv:hep-ph/0606054v3, eq. 3.15, for example. 
      if(anti>0) *gsl_matrix_complex_ptr(d->H_GSL, i, i) = gsl_complex_add_real(gsl_matrix_complex_get(d->H_GSL,i,i) ,  kr2GNe/2);
      else       *gsl_matrix_complex_ptr(d->H_GSL, i, i) = gsl_complex_add_real(gsl_matrix_complex_get(d->H_GSL,i,i) , -kr2GNe/2);;
    }
  }
  // Add nue CC coherent forward scattering from sterile neutrinos. 
  if(anti>0) *gsl_matrix_complex_ptr(d->H_GSL, 0, 0) = gsl_complex_add_real(gsl_matrix_complex_get(d->H_GSL,0,0) ,  kr2GNe);
  else       *gsl_matrix_complex_ptr(d->H_GSL, 0, 0) = gsl_complex_add_real(gsl_matrix_complex_get(d->H_GSL,0,0) , -kr2GNe);

  // Solve Hamiltonian for eigensystem
  gsl_eigen_hermv(d->H_GSL, d->fEval, d->fEvec, d->W_GSL);
  
}

///.....................................................................
///
/// Propagate the current neutrino state over a distance L in km
/// with an energy E in GeV through constant matter of density
/// Ne in mole/cm^3. 
/// @param anti - +1 = neutrino case, -1 = anti-neutrino case
///
void PMNS_Sterile::PropMatter(double L, double E, double Ne, int anti) 
{

  // Solve Hamiltonian
  this->SolveHam(E, Ne, anti);

  // Store coefficients of propagation eigenstates
  vector<complex> nuComp(fNumNus, zero);
  for(int i=0;i<fNumNus;i++){
    nuComp[i] = 0;
    for(int j=0;j<fNumNus;j++){
      gsl_complex buf = gsl_matrix_complex_get(d->fEvec,j,i);
      complex evecji = complex( GSL_REAL(buf), GSL_IMAG(buf) );
      nuComp[i] += fNuState[j] * conj( evecji );
    }
  }

  // Propagate neutrino state
  for(int i=0;i<fNumNus;i++){
    fNuState[i] = 0;
    for(int j=0;j<fNumNus;j++){
      gsl_complex buf = gsl_matrix_complex_get(d->fEvec,i,j);
      complex evecij = complex( GSL_REAL(buf), GSL_IMAG(buf) );
      complex iEval(0.0,gsl_vector_get(d->fEval,j));
      fNuState[i] +=  exp(-iEval * kKm2eV*L) * nuComp[j] * evecij;
    }
  }

}

//......................................................................
///
/// Do several layers in a row. L and Ne must have the same length
///
void PMNS_Sterile::PropMatter(const std::list<double>& L,
                              double                   E,
                              const std::list<double>& Ne,
                              int anti)
{
  if (L.size()!=Ne.size()) abort();
  std::list<double>::const_iterator Li  (L.begin());
  std::list<double>::const_iterator Lend(L.end());
  std::list<double>::const_iterator Ni  (Ne.begin());
  for (; Li!=Lend; ++Li, ++Ni) {
    //cout << *Li << " km and " << *Ni << " mol/cm^3" << endl;
    this->PropMatter(*Li, E, *Ni, anti);
  }
}


//......................................................................
///
/// Reset the neutrino state back to a pure flavour where
/// it starts
///
void PMNS_Sterile::ResetToFlavour(int flv) 
{
  int i;
  for (i=0; i<fNumNus; ++i){
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
double PMNS_Sterile::P(int flv) const
{
  assert(flv>=0 && flv<fNumNus);
  return norm(fNuState[flv]);
}

////////////////////////////////////////////////////////////////////////
