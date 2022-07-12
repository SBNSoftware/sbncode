/**
 * \file NPEFlashFilter.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class NPEFlashFilter
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_NPEFLASHFILTER_H
#define OPT0FINDER_NPEFLASHFILTER_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseFlashFilter.h"
#include "flashmatch/Base/FlashFilterFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseFlashFilter.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashFilterFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif

namespace flashmatch {
  /**
     \class NPEFlashFilter
     Type flashmatch::BaseFlashFilter algorithm class. It applies a very simple  \n
     filter based on number of photo-electron count.
  */
  class NPEFlashFilter : public BaseFlashFilter{
    
  public:
    
    /// Default constructor
    NPEFlashFilter(const std::string name="NPEFlashFilter");
    
    /// Default destructor
    ~NPEFlashFilter(){}

    /// Implementation of a virtual function
    IDArray_t Filter(const FlashArray_t&);

  protected:

    void _Configure_(const Config_t &pset);

  private:

    double _npe_threshold;    ///< threshold [p.e.]: to ignore any flash below this value
    
  };

  /**
     \class flashmatch::NPEFlashFilterFactory
  */
  class NPEFlashFilterFactory : public FlashFilterFactoryBase {
  public:
    /// ctor
    NPEFlashFilterFactory() { FlashFilterFactory::get().add_factory("NPEFlashFilter",this); }
    /// dtor
    ~NPEFlashFilterFactory() {}
    /// creation method
    BaseFlashFilter* create(const std::string instance_name) { return new NPEFlashFilter(instance_name); }
  };
  
}

#endif
/** @} */ // end of doxygen group 

