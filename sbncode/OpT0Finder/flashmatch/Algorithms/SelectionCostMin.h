/**
 * \file SelectionCostMin.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class SelectionCostMin
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_SelectionCostMin_H
#define OPT0FINDER_SelectionCostMin_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseMatchSelection.h"
#include "flashmatch/Base/MatchSelectionFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseMatchSelection.h"
#include "sbncode/OpT0Finder/flashmatch/Base/MatchSelectionFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif

namespace flashmatch {
  /**
     \class SelectionCostMin
  */
  class SelectionCostMin : public BaseMatchSelection{
    
  public:
    
    /// Default constructor
    SelectionCostMin(const std::string name="SelectionCostMin");
    
    /// Default destructor
    ~SelectionCostMin(){}

    std::vector<FlashMatch_t> Select(const std::vector<std::vector<FlashMatch_t> >& match_data);

  protected:

    void _Configure_(const Config_t &pset);

  private:
    
    size_t _min_num_pt; ///< mininum number of QPoint_t to pass the filter
    
  };

  /**
     \class flashmatch::SelectionCostMinFactory
  */
  class SelectionCostMinFactory : public MatchSelectionFactoryBase {
  public:
    /// ctor
    SelectionCostMinFactory() { MatchSelectionFactory::get().add_factory("SelectionCostMin",this); }
    /// dtor
    ~SelectionCostMinFactory() {}
    /// creation method
    BaseMatchSelection* create(const std::string instance_name) { return new SelectionCostMin(instance_name); }
  };
  
}

#endif
/** @} */ // end of doxygen group 

