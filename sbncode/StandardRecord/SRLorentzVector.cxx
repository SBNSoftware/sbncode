////////////////////////////////////////////////////////////////////////
// \file    SRLorentzVector.cxx
// \author  Christopher Backhouse - bckhouse@caltech.edu
////////////////////////////////////////////////////////////////////////
 
#include "SRLorentzVector.h"
 
namespace caf
{
  SRLorentzVector::SRLorentzVector() :
    E (std::numeric_limits<float>::signaling_NaN()),
    px(std::numeric_limits<float>::signaling_NaN()),
    py(std::numeric_limits<float>::signaling_NaN()),
    pz(std::numeric_limits<float>::signaling_NaN())
  {
  }

  SRLorentzVector::SRLorentzVector(const TLorentzVector& v)
    : E(v.E()), px(v.X()), py(v.Y()), pz(v.Z())
  {
  }

  SRLorentzVector::~SRLorentzVector()
  {
  }

  SRLorentzVector::operator TLorentzVector() const
  {
   return TLorentzVector(px, py, pz, E);
  }
} // end namespace caf
////////////////////////////////////////////////////////////////////////
