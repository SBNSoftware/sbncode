////////////////////////////////////////////////////////////////////////
// \file    SRVector3D.cxx
// \brief   TVector3 have three doubles plus ROOT IDs. This adds up. 
//          SRVector3D contain only the 3 floats we need.
// \author  Christopher Backhouse - bckhouse@caltech.edu
////////////////////////////////////////////////////////////////////////

 #include "SRVector3D.h"
 
 namespace caf
 {
  SRVector3D::SRVector3D() :
    x(std::numeric_limits<float>::signaling_NaN()),
    y(std::numeric_limits<float>::signaling_NaN()),
    z(std::numeric_limits<float>::signaling_NaN())
   {
   }

   SRVector3D::SRVector3D(float _x, float _y, float _z) :
    x(_x), y(_y), z(_z)
  {
  }

  SRVector3D::SRVector3D(const TVector3& v) :
    x(v.X()), y(v.Y()), z(v.Z())
   {
   }
 
   SRVector3D::~SRVector3D()
   {
   }

   void SRVector3D::SetXYZ(float _x, float _y, float _z)
   {
     x = _x;
     y = _y;
     z = _z;
   }
 
   SRVector3D::operator TVector3() const
   {
     return TVector3(x, y, z);
   }
 
 } // end namespace caf
////////////////////////////////////////////////////////////////////////
