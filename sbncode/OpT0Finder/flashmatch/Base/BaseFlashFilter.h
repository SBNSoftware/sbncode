/**
 * \file BaseFlashFilter.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseFlashFilter
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef OPT0FINDER_BASEFLASHFILTER_H
#define OPT0FINDER_BASEFLASHFILTER_H

#include "BaseAlgorithm.h"

namespace flashmatch {
  /**
     \class BaseFlashFilter
     Algorithm base class for filtering out optical flashes from a candidate list \n
     to be matched with TPC object.
  */
  class BaseFlashFilter : public BaseAlgorithm {
    
  public:
    
    /// Default constructor
    BaseFlashFilter(const std::string name="noname") : BaseAlgorithm(kFlashFilter,name)
    {}
    
    /// Default destructor
    virtual ~BaseFlashFilter(){}

    /**
       CORE FUNCTION: takes in a list of flash and returns a list of flash TO BE USED for matching. \n
     */
    virtual IDArray_t Filter(const FlashArray_t&) = 0;
    
  };
}
  
#endif
/** @} */ // end of doxygen group 

