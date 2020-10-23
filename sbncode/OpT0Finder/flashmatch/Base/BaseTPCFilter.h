/**
 * \file BaseTPCFilter.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseTPCFilter
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef OPT0FINDER_BASETPCFILTER_H
#define OPT0FINDER_BASETPCFILTER_H

#include "BaseAlgorithm.h"

namespace flashmatch {
  /**
     \class BaseTPCFilter
     Algorithm base class for filtering out optical TPC objects from a candidate list \n
     to be matched with TPC object.
  */
  class BaseTPCFilter : public BaseAlgorithm{
    
  public:
    
    /// Default constructor
    BaseTPCFilter(const std::string name="noname") : BaseAlgorithm(kTPCFilter,name)
    {}
 
    /// Default destructor
    virtual ~BaseTPCFilter(){}
    
    /**
       CORE FUNCTION: takes in a list of TPC objects and returns a list of TPC object TO BE USED for matching. \n
     */
    virtual IDArray_t Filter(const QClusterArray_t&) = 0;
    
  };
}

#endif
/** @} */ // end of doxygen group 

