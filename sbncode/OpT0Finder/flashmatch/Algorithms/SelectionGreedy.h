/**
 * \file SelectionGreedy.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class SelectionGreedy
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_SelectionGreedy_H
#define OPT0FINDER_SelectionGreedy_H

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
     \class SelectionGreedy
  */
  class SelectionGreedy : public BaseMatchSelection {
    
  public:
    
    /// Default constructor
    SelectionGreedy(const std::string name="SelectionGreedy");
    
    /// Default destructor
    ~SelectionGreedy(){}

   std::vector<FlashMatch_t> Select(const std::vector<std::vector<FlashMatch_t> >& match_data);

  protected:

    void _Configure_(const Config_t &pset);

  private:
    
    bool   _allow_reuse_flash;    ///< allow one flash to be matched against multiple tpc object
    double _score_max_threshold;  ///< ignore touch-matched pairs with a score larger than this threshold
    double _score_min_threshold;  ///< ignore flash-matched pairs with a score less than this threshold
    double _score_max_ceiling;    ///< treat matched pairs with score values higher than this ceiling same (reset to this value)
  };

  /**
     \class flashmatch::SelectionGreedyFactory
  */
  class SelectionGreedyFactory : public MatchSelectionFactoryBase {
  public:
    /// ctor
    SelectionGreedyFactory() { MatchSelectionFactory::get().add_factory("SelectionGreedy",this); }
    /// dtor
    ~SelectionGreedyFactory() {}
    /// creation method
    BaseMatchSelection* create(const std::string instance_name) { return new SelectionGreedy(instance_name); }
  };
  
}

#endif
/** @} */ // end of doxygen group 

