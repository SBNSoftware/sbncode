#ifndef BASICTOOL_GEOTRAJECTORY_CXX
#define BASICTOOL_GEOTRAJECTORY_CXX

#include "GeoTrajectory.h"
#include <sstream>
namespace geoalgo {

  Trajectory::Trajectory(size_t npoints, size_t ndimension) 
    : std::vector<geoalgo::Point_t>(npoints, Point_t(ndimension))
  {}

  Trajectory::Trajectory(const std::vector<std::vector<double> > &obj)
  {
    this->reserve(obj.size());
    for(auto const& p : obj) this->push_back(Point_t(p));
  }

  Trajectory::Trajectory(const std::vector<geoalgo::Point_t> &obj)
  {
    this->reserve(obj.size());
    for(auto const& p : obj) this->push_back(p);
  }

  double Trajectory::Length(size_t start_step,size_t end_step) const {

    if(end_step == 0) end_step = size() - 1; // By default end_step is 0. Then consider the whole trajectory()

    // Sanity checks
    if(start_step >= end_step) throw GeoAlgoException("Cannot have start step >= end step!");

    if(end_step >= size()) throw GeoAlgoException("Requested step index bigger than size!");

    // if length < 2, no length
    if(size()<2) return 0;
    
    double length = 0;
    for(size_t i=start_step; i<end_step; ++i)
      
      length += (*this)[i]._Dist_((*this)[i+1]);
    
    return length;
  }

  bool Trajectory::IsLonger(double ref) const {

    if(size()<2) return false;

    double length = 0;
    for(size_t i=0; i<size()-1; ++i) {

      length += (*this)[i]._Dist_((*this)[i+1]);

      if(length > ref) return true;
    }

    return false;    
  }

  void Trajectory::push_back(const Point_t& obj) {
    compat(obj); 
    if (!(size() && obj == (*rbegin())))
      std::vector<geoalgo::Point_t>::push_back(obj);
  }

  void Trajectory::compat(const Point_t& obj) const {
    
    if(!size()) return;
    if( (*(this->begin())).size() != obj.size() ) {
      
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
	  << " size mismatch: "
	  << (*(this->begin())).size() << " != " << obj.size() << std::endl;
      throw GeoAlgoException(msg.str());
    }
  }

  void Trajectory::compat(const Trajectory &obj) const {
    
    if(!size() || !(obj.size())) return;
    
    if( (*(this->begin())).size() != (*obj.begin()).size() ) {
      
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
	  << " size mismatch: "
	  << (*(this->begin())).size() << " != " << (*obj.begin()).size() << std::endl;
      throw GeoAlgoException(msg.str());
      
    }
  }

  Vector Trajectory::Dir(size_t i) const {
    
    if(size() < (i+2)) {
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
	  << " length=" << size() << " is too short to find a direction @ index=" << i << std::endl;
      throw GeoAlgoException(msg.str());
    }
    return _Dir_(i);
  }
  
  Vector Trajectory::_Dir_(size_t i) const {

    return ((*this)[i+1] - (*this)[i]);
    
  }
}

#endif


