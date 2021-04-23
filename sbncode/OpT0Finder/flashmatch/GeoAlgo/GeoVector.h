/**
 * \file GeoObjects.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class Point and Vector
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOVECTOR_H
#define BASICTOOL_GEOVECTOR_H

#include "GeoAlgoConstants.h"
#include "GeoAlgoException.h"
#include <TVector3.h>
#include <TLorentzVector.h>

namespace geoalgo {

  // Forward declaration (don't worry)
  class Trajectory;
  class LineSegment;
  class HalfLine;
  class Sphere;
  class GeoAlgo;

  /**
     \class Vector
     This class represents an n-dimensional vector
  */
  class Vector : public std::vector<double> {
    friend class Trajectory;
    friend class HalfLine;
    friend class LineSegment;
    friend class Sphere;
    friend class GeoAlgo;
  public:
    /// Default ctor
    Vector() : std::vector<double>()
    {}

    /// Ctor to instantiate with invalid value
    Vector(size_t n) : std::vector<double>(n,kINVALID_DOUBLE)
    {}

    /// Default ctor w/ a bare std::vector<double>
    Vector(const std::vector<double> &obj) : std::vector<double>(obj)
    {}


    Vector(const double x, const double y);                 ///< ctor w/ x & y
    Vector(const double x, const double y, const double z); ///< ctor w/ x, y & z
    Vector(const TVector3 &pt);                             ///< ctor w/ TVector3
    Vector(const TLorentzVector &pt);                       ///< ctor w/ TLorentzVector
    
    virtual ~Vector(){} ///< Default dtor

    void   Normalize(); ///< Normalize itself
    bool   IsValid () const; ///< Check if point is valid    
    double SqLength() const; ///< Compute the squared length of the vector
    double Length  () const; ///< Compute the length of the vector
    Vector Dir     () const; ///< Return a direction unit vector
    double Phi	   () const; ///< Compute the angle Phi
    double Theta   () const; ///< Compute the angle theta

    double SqDist(const Vector &obj) const; ///< Compute the squared distance to another vector
    double Dist  (const Vector& obj) const; ///< Compute the distance to another vector
    double Dot   (const Vector &obj) const; /// Compute a dot product of two vectors
    Vector Cross (const Vector &obj) const; /// Compute a cross product of two vectors
    double Angle (const Vector &obj) const; /// Compute an opening angle w.r.t. the given vector

    TLorentzVector ToTLorentzVector() const; ///< Convert geovector to TLorentzVector (with 4th element set equal to 0)

    /// Dimensional check for a compatibility
    void compat(const Vector& obj) const;

    /// rotation operations
    void RotateX(const double& theta);
    void RotateY(const double& theta);
    void RotateZ(const double& theta);

    std::string dump() const;

  protected:

    /// Compute the squared-distance to another vector w/o dimension check
    double _SqDist_(const Vector& obj) const;
    /// Compute the distance to another vector w/o dimension check
    double _Dist_(const Vector& obj) const;
    /// Compute a dot product w/o dimention check.
    double _Dot_(const Vector& obj) const;
    /// Compute a cross product w/o dimension check.
    Vector _Cross_(const Vector& obj) const;
    /// Compute the angle in degrees between 2 vectors w/o dimension check.
    double _Angle_(const Vector& obj) const;

  public:
    //
    // binary/uniry operators
    //
    inline Vector& operator+=(const Vector& rhs) { 
      for(size_t i=0; i<size(); ++i) (*this)[i] += rhs[i];
      return *this;
    }

    inline Vector& operator-=(const Vector& rhs) {
      for(size_t i=0; i<size(); ++i) (*this)[i] -= rhs[i];
      return *this;
    }

    inline Vector& operator*=(const double rhs) {
      for(auto &v : *this) v *= rhs;
      return *this;
    }

    inline Vector& operator/=(const double rhs) {
      for(auto &v : *this) v /= rhs;
      return *this;
    }

    inline Vector& operator=(const Vector& rhs) {
      this->resize(rhs.size());
      for(size_t i=0; i<rhs.size(); ++i) (*this)[i]=rhs[i];
      return (*this);
    }

    inline Vector operator+(const Vector& rhs) const
    { 
      Vector res((*this));
      res+=rhs;
      return res;
    }

    inline Vector operator-(const Vector& rhs) const
    { 
      Vector res((*this));
      res-=rhs;
      return res;
    }

    inline double operator*(const Vector& rhs) const
    { 
      double res=0;
      for(size_t i=0; i<size(); ++i) res += (*this)[i] * rhs[i];
      return res;
    }

    inline Vector operator*(const double& rhs) const
    {
      Vector res((*this));
      res *= rhs;
      return res;
    }

    inline Vector operator/(const double& rhs) const
    {
      Vector res((*this));
      res /= rhs;
      return res;
    }

    inline bool operator< ( const Vector& rhs ) const 
    { 
      compat(rhs);
      for(size_t i=0; i<size(); ++i)
	if((*this)[i] < rhs[i]) return true;
      return false;
    }

    inline bool operator< ( const double& rhs) const 
    { return Length() < rhs; }

    inline bool operator== ( const Vector& rhs) const
    { 
      compat(rhs); 
      for(size_t i=0; i<size(); ++i) 
	if((*this)[i] != rhs[i]) return false;
      return true;
    }

    inline bool operator!= ( const Vector& rhs) const
    { return !(*this == rhs); }

    /// Streamer
    #ifndef __CINT__
    friend std::ostream& operator << (std::ostream &o, ::geoalgo::Vector const& a)
    { o << a.dump(); return o; }
    #endif

  };

  /// Point has same feature as Vector
  typedef Vector Vector_t;
  typedef Vector Point_t;
}

// Define a pointer comparison                                                                                                                     
#ifndef __GCCXML__
namespace std {
  template <>
  class less<geoalgo::Vector*>
  {
  public:
    bool operator()( const geoalgo::Vector* lhs, const geoalgo::Vector* rhs )
    { return (*lhs) < (*rhs); }
  };
}
#endif

#endif
/** @} */ // end of doxygen group 

