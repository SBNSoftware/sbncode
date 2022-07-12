/**
 * \file TimeCompatMatch.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class TimeCompatMatch
 *
 * @author david caratelli
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_TIMECOMPATMATCH_H
#define OPT0FINDER_TIMECOMPATMATCH_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0

#include "flashmatch/Base/BaseProhibitAlgo.h"
#include "flashmatch/Base/FlashProhibitFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseProhibitAlgo.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashProhibitFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif

#include <cmath>
#include <sstream>


namespace flashmatch {
  
  /**
     \class TimeCompatMatch
     Simple flash matching algorithm. Based on absolute time of flash and track
     w.r.t. trigger time, if the two objects are incompatible (because the time-difference
     is larger than a full drift window) the match is not allowed
  */
  class TimeCompatMatch : public BaseProhibitAlgo {
    
  public:
    
    /// Default constructor
    TimeCompatMatch(const std::string name="TimeCompatMatch");
    
    /// Default destructor
    ~TimeCompatMatch(){}

    bool MatchCompatible(const QCluster_t& clus, const Flash_t& flash);

  protected:

    void _Configure_(const Config_t &pset);

  private:

    /// Buffer time to allow some uncertainty [us]
    double _time_window;
    /// Shift in beam timing w.r.t. flash's time reference =0
    double _time_shift;

  };

  /**
     \class flashmatch::TimeCompatMatchFactory
  */
  class TimeCompatMatchFactory : public FlashProhibitFactoryBase {
  public:
    /// ctor
    TimeCompatMatchFactory() { FlashProhibitFactory::get().add_factory("TimeCompatMatch",this); }
    /// dtor
    ~TimeCompatMatchFactory() {}
    /// creation method
    BaseProhibitAlgo* create(const std::string instance_name) { return new TimeCompatMatch(instance_name); }
  };
  
}
#endif
/** @} */ // end of doxygen group 

