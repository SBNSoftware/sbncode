#include <string>
#include <vector>
#include <cmath>
#include "TVector3.h"
#include "Plane.hh"

namespace selection{

  // Plane constructor
  Plane::Plane(const TVector3 &V, const TVector3 &A, const TVector3 &B) :
    m_V(V),
    m_A(A),
    m_B(B){
      m_a = ((1/m_A.Mag()) * m_A);
      m_b = ((1/m_B.Mag()) * m_B);
      m_n = m_a.Cross(m_b);  
      m_alpha = m_A.Mag();
      m_beta  = m_B.Mag(); 
  }

  // Plane object getters
  TVector3 Plane::GetV() const {return m_V;}
  TVector3 Plane::GetA() const {return m_A;}
  TVector3 Plane::GetB() const {return m_B;}
  TVector3 Plane::GetUnitA() const {return m_a;}
  TVector3 Plane::GetUnitB() const {return m_b;} 
  TVector3 Plane::GetUnitN() const {return m_n;}
  float Plane::GetAlpha() const {return m_alpha;}
  float Plane::GetBeta()  const {return m_beta;}
} // Selection
