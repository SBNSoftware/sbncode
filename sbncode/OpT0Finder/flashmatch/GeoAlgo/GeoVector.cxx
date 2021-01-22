#ifndef BASICTOOL_GEOVECTOR_CXX
#define BASICTOOL_GEOVECTOR_CXX

#include "GeoVector.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <cmath>
#include <limits>

namespace geoalgo {

  Vector::Vector(const double x, const double y) : Vector(2)
  { (*this)[0] = x; (*this)[1] = y; }
  
  Vector::Vector(const double x, const double y, const double z) : Vector(3)
  { (*this)[0] = x; (*this)[1] = y; (*this)[2] = z; }
  
  Vector::Vector(const TVector3 &pt) : Vector(3)
  { (*this)[0] = pt[0]; (*this)[1] = pt[1]; (*this)[2] = pt[2]; }
  
  Vector::Vector(const TLorentzVector &pt) : Vector(3)
  { (*this)[0] = pt[0]; (*this)[1] = pt[1]; (*this)[2] = pt[2]; }

  bool Vector::IsValid() const {
    
    for (auto const &v : (*this)){
      // if any point is different from kINVALID_DOUBLE
      // then the point is valid
      if (v != kINVALID_DOUBLE)
	return true;
    }
    
    return false;
  }

  double Vector::SqLength() const {
    double res=0;
    for(auto const &v : (*this)) res += v*v;
    return res;
  }

  double Vector::Length() const
  { return sqrt(SqLength()); }
  
  double Vector::SqDist(const Vector &obj) const {
    compat(obj);
    return _SqDist_(obj);
  }

  double Vector::Dist(const Vector& obj) const 
  { return sqrt(SqDist(obj)); }

  double Vector::Dot(const Vector &obj) const {
    compat(obj);
    return _Dot_(obj);
  }
    
  Vector Vector::Cross(const Vector &obj) const {
    
    if(size()!=3 || obj.size()!=3)
      
      throw GeoAlgoException("<<Cross>> only possible for 3-dimensional vectors!");
    
    return _Cross_(obj);
  }

  double Vector::Phi() const {
    return (*this)[0] == 0.0 && (*this)[1] == 0.0 ? 0.0 : atan2((*this)[1],(*this)[0]);
  }
  
  double Vector::Theta() const {
    if ( size() != 3 )
      throw GeoAlgoException("<<Theta>> Only possible for 3-dimensional vectors!");

    return (*this).Length() == 0.0 ? 0.0 : acos( (*this)[2] / (*this).Length() );
  }

  double Vector::Angle(const Vector &obj) const {
    compat(obj);
    if(size()!=2 && size()!=3)
      throw GeoAlgoException("<<Angle>> only possible for 2 or 3-dimensional vectors!");
    return _Angle_(obj);
  }
  
  TLorentzVector Vector::ToTLorentzVector() const {
    if(size()!=3)
      throw GeoAlgoException("<<ToTLorentsVector>> only possible for 3-dimensional vectors!");
    return TLorentzVector((*this)[0],(*this)[1],(*this)[2],0.);
  }

  void Vector::Normalize() { (*this) /= this->Length(); }

  Vector Vector::Dir() const {
    Vector res(*this);
    res /= res.Length();
    return res;
  }
      
  void Vector::compat(const Vector& obj) const {
    if(size() != obj.size()) {
      std::ostringstream msg;
	msg << "<<" << __FUNCTION__ << ">>" 
	    << " size mismatch: "
	    << size() << " != " << obj.size()
	    << std::endl;
	throw GeoAlgoException(msg.str());
    }
  }

  double Vector::_SqDist_(const Vector& obj) const
  {
    double dist = 0;
    for(size_t i=0; i<size(); ++i) dist += ((*this)[i] - obj[i]) * ((*this)[i] - obj[i]);
      return dist;
  }

  double Vector::_Dist_(const Vector& obj) const
  { return sqrt(_SqDist_(obj)); }

  double Vector::_Dot_(const Vector& obj) const
    { return (*this) * obj; }

  Vector Vector::_Cross_(const Vector& obj) const
  {
    Vector res(3);
    res[0] = (*this)[1] * obj[2] - obj[1] * (*this)[2];
    res[1] = (*this)[2] * obj[0] - obj[2] * (*this)[0];
    res[2] = (*this)[0] * obj[1] - obj[0] * (*this)[1];
    return res;
  }    

  double Vector::_Angle_(const Vector& obj) const
  { return acos( _Dot_(obj) / Length() / obj.Length() ); }
  

  void Vector::RotateX(const double& theta)
  {

    double c = cos(theta);
    double s = sin(theta);
    
    double ynew = (*this)[1] * c - (*this)[2] * s;
    double znew = (*this)[1] * s + (*this)[2] * c;

    (*this)[1] = ynew;
    (*this)[2] = znew;
    
    return;
  }


  void Vector::RotateY(const double& theta)
  {

    double c = cos(theta);
    double s = sin(theta);
    
    double xnew =   (*this)[0] * c + (*this)[2] * s;
    double znew = - (*this)[0] * s + (*this)[2] * c;
    
    (*this)[0] = xnew;
    (*this)[2] = znew;

    return;
  }


  void Vector::RotateZ(const double& theta)
  {

    double c = cos(theta);
    double s = sin(theta);
    
    double xnew = (*this)[0] * c - (*this)[1] * s;
    double ynew = (*this)[0] * s + (*this)[1] * c;

    (*this)[0] = xnew;
    (*this)[1] = ynew;
    
    return;
  }

  std::string Vector::dump() const {
    
    std::string msg="Pt (";
    for(size_t i=0; i<this->size(); ++i)
      msg += (std::to_string((*this)[i]) + ( (i+1) == this->size()  ? ")" : ","));
    return msg;
  }

}

#endif


