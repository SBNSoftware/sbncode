/**
 * \file TimeRange.h
 *
 * \ingroup Algorithm
 * 
 * \brief Class def header for a class TimeRange
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithm

    @{*/
#ifndef OPT0FINDER_TIMERANGE_H
#define OPT0FINDER_TIMERANGE_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>

namespace flashmatch {
 

  /**
     \class TimeRange
     User defined class TimeRange ... these comments are used to generate
     doxygen documentation!
  */
  class TimeRange {
    
  public:
    
    /// Default constructor
    TimeRange(double start=0,double end=0)
    { SetRange(start,end); }

    /// Default destructor
    ~TimeRange(){}

    void SetRange(double start,double end)
    {
      if(start >= end) throw std::exception();
      _start=start;
      _end=end;
    }

    double Start() const { return _start; }

    double End() const { return _end; }

    inline bool operator<(const TimeRange& rhs) const {
      if(_end < rhs.Start()) return true;
      if(rhs.End() < _start) return false;
      return false;
    }

    inline bool operator<(const double& rhs) const {
      if(_end < rhs) return true;
      if(rhs < _start ) return false;
      return false;
    }

  private:

    double _start;
    double _end;
    
  };

  class TimeRangeSet {

  public:
    TimeRangeSet(){}

    ~TimeRangeSet(){}

    //int RangeIndex(double time);

    void Print() const
    {
      for(auto const& r : _time_range_v)
	std::cout<<r.Start()<<" => "<<r.End()<<std::endl;
    }
    
    bool Overlap(const flashmatch::TimeRange& range) const
    {
      auto low = std::lower_bound(_time_range_v.begin(),
				  _time_range_v.end(),
				  range);
      if(low == _time_range_v.end()) return false;
      if(range.End() >= (*low).Start()) return true;
      return false;
    }

    bool Overlap(double time) const
    { auto low = std::lower_bound(_time_range_v.begin(),
				  _time_range_v.end(),
				  time);
      if(low == _time_range_v.end()) return false;
      return ((*low).Start() <= time && time <= (*low).End());
    }

    void Insert(const flashmatch::TimeRange& range)
    {
      auto low = std::lower_bound(_time_range_v.begin(),
				  _time_range_v.end(),
				  range);

      if(low == _time_range_v.end()) {
	_time_range_v.push_back(range);
	return;
      }

      // Check if found one is overlapping
      if(range.End() < (*low).Start()){

	double next_start=0;
	double next_end  =0;
	double start=range.Start();
	double end=range.End();
	while(low != _time_range_v.end()) {
	  next_start = (*low).Start();
	  next_end   = (*low).End();
	  (*low).SetRange(start,end);
	  start = next_start;
	  end = next_end;
	  ++low;
	}
	_time_range_v.push_back(TimeRange(start,end));
      }else{

	if((*low).Start() <= range.End()) {
	  
	  (*low).SetRange(std::min((*low).Start(),range.Start()),
			  std::max((*low).End(),range.End()));

	  // If this overlaps with previous window, set it to previous window's max boundary
	  size_t dist = std::distance(_time_range_v.begin(),low);
	  if(dist) {
	    auto prev = low-1;
	    if((*prev).End() > (*low).Start()) (*low).SetRange((*prev).End(),(*low).End());
	  }
	  
	  // If this overlaps with next window, set it to next window's min boundary
	  if( (dist+1) < _time_range_v.size() ) {
	    auto next = low+1;
	    if((*low).End() > (*next).Start()) (*low).SetRange((*low).Start(),(*next).Start());
	  }
	}
      }
    }
    
  private:

    std::vector<flashmatch::TimeRange> _time_range_v;
    
  };
}

namespace std {
  template <>
  class less<flashmatch::TimeRange*>
  {
  public:
    bool operator()( const flashmatch::TimeRange* lhs, const flashmatch::TimeRange* rhs )
    { return (*lhs) < (*rhs); }
  };
}

#endif
/** @} */ // end of doxygen group 

