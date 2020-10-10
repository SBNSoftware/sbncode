/**
 * \file FilterArray.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class FilterArray
 *
 * @author david caratelli
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_FILTERARRAY_H
#define OPT0FINDER_FILTERARRAY_H

#include "OpT0Finder/Base/BaseTPCFilter.h"
namespace flashana {
  /**
     \class FilterArray
     Implementation of flashana::BaseTPCFilter abstract algorithm class. \n
     It applies a _very_ simple filtering: excludes TPC objects (flashana::QCluster_t) \n
     that contains less than specified number of points. \n
  */
  class FilterArray : public BaseTPCFilter{
    
  public:
    
    /// Default constructor
    FilterArray();
    
    /// Default destructor
    ~FilterArray(){}

    /// Implementation of virtualfunction
    IDArray_t Filter(const QClusterArray_t&);

    /**
     * @brief append a filter algorithm to the array
     */						
    void AppendFilterAlgo(flashana::BaseTPCFilter* filter) { _filter_v.push_back(filter); }

  private:

    // array of filter algorithms to be applied
    std::vector<flashana::BaseTPCFilter*> _filter_v;
    
  };
}

#endif
/** @} */ // end of doxygen group 

