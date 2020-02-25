////////////////////////////////////////////////////////////////////////
// \file    SRLorentzVector.h
// \author  Christopher Backhouse - bckhouse@caltech.edu
////////////////////////////////////////////////////////////////////////

#ifndef SRLORENTZVECTOR_H
#define SRLORENTZVECTOR_H
 
#ifndef __GCCXML__
#include "TLorentzVector.h"
#include "TVector3.h"
#endif
 
namespace caf
{
  /// 4-vector with more efficient storage than TLorentzVector
  class SRLorentzVector
  {
  public:
    SRLorentzVector();
    virtual ~SRLorentzVector();

#ifndef __GCCXML__
    SRLorentzVector(const TLorentzVector& v);

    /// Recommend users convert back to TLorentzVector for boosts etc
    operator TLorentzVector() const;

    // For access as a position vector. For momentum use the member variables
    // directly.
    float T() const {return E;}
    float X() const {return px;}
    float Y() const {return py;}
    float Z() const {return pz;}
    float Mag() const {return TMath::Sqrt(px*px + py*py + pz*pz);}
    float Beta() const {return Mag()/E;}
    float Gamma() const {return 1.0/TMath::Sqrt(1-Beta()*Beta());}
 
   TVector3 Vect() const {return TVector3(px, py, pz);}
 #endif
 
    float E;
    float px;
    float py;
    float pz;
  };

 } // end namespace

#endif // SRLORENTZVECTOR_H
//////////////////////////////////////////////////////////////////////////////
