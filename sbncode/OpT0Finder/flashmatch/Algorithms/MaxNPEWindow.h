/**
 * \file MaxNPEWindow.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class MaxNPEWindow
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_MAXNPEWINDOW_H
#define OPT0FINDER_MAXNPEWINDOW_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseFlashFilter.h"
#include "flashmatch/Base/FlashFilterFactory.h"
#else
#include "sbndcode/OpT0Finder/flashmatch/Base/BaseFlashFilter.h"
#include "sbndcode/OpT0Finder/flashmatch/Base/FlashFilterFactory.h"
#endif

namespace flashmatch {
  /**
     \class MaxNPEWindow
     Type flashmatch::BaseFlashFilter algorithm class. It applies a very simple  \n
     filter based on time window to eliminate multiple flashes coming from the \n
     same interaction. The algorithm works in the following steps:             \n

     0) Order a list of provided flash using the amount photoelectrons \n
     (the bigger, the higher priority). \n

     1) Loop over the list of flashes in high=>low priority order. \n

     1.1) For each flash, attempt to open a time window T1 => T2 w.r.t. \n
     flash time where T1 & T2 are configured via user \n

     1.2) If a flash falls within any of previously-opened time window, \n
     ignore this flash. \n

     1.3) Terminate the loop when the number of photoelectrons drops below \n
     the threshold set by a user.
  */
  class MaxNPEWindow : public BaseFlashFilter{
    
  public:
    
    /// Default constructor
    MaxNPEWindow(const std::string name="MaxNPEWindow");
    
    /// Default destructor
    ~MaxNPEWindow(){}

    /// Implementation of a virtual function
    IDArray_t Filter(const FlashArray_t&);

  protected:

    void _Configure_(const Config_t &pset);

  private:

    double _time_lower_bound; ///< T1 [us]: the lower edge of the opened time window
    double _time_upper_bound; ///< T2 [us]: the upper edge of the opened time window
    double _npe_threshold;    ///< threshold [p.e.]: to ignore any flash below this value
    
  };

  /**
     \class flashmatch::MaxNPEWindowFactory
  */
  class MaxNPEWindowFactory : public FlashFilterFactoryBase {
  public:
    /// ctor
    MaxNPEWindowFactory() { FlashFilterFactory::get().add_factory("MaxNPEWindow",this); }
    /// dtor
    ~MaxNPEWindowFactory() {}
    /// creation method
    BaseFlashFilter* create(const std::string instance_name) { return new MaxNPEWindow(instance_name); }
  };
  
}

#endif
/** @} */ // end of doxygen group 

