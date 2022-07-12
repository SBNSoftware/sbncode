/**
 * \file AnodeMatch.h
 *
 * \ingroup Algorithms
 *
 * \brief Class def header for a class AnodeMatch
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_AnodeMatch_H
#define OPT0FINDER_AnodeMatch_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseTouchMatch.h"
#include "flashmatch/Base/TouchMatchFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseTouchMatch.h"
#include "sbncode/OpT0Finder/flashmatch/Base/TouchMatchFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif

namespace flashmatch {

  /**
     \class AnodeMatch
  */
  class AnodeMatch : public BaseTouchMatch {

  public:

    /// Default constructor
    AnodeMatch(const std::string name="AnodeMatch");

    /// Default destructor
    ~AnodeMatch(){}

    void Match(const QCluster_t&, const Flash_t&, FlashMatch_t& match);

  protected:

    void _Configure_(const Config_t &pset);

  private:
    double _time_window; ///< maximum time-diff from the hypothesized timing from TPC to flash
    double _space_range; ///< spatial extent from the edge TPC track to inspect PMTs to see if they see enough PE to claim a match.
    double _pe_threshold; ///< minimum PE required for a PMT that is within the range of _space_range_squared from the track edge points. 
  };

  /**
     \class flashmatch::AnodeMatchFactory
  */
  class AnodeMatchFactory : public TouchMatchFactoryBase {
  public:
    /// ctor
    AnodeMatchFactory() { TouchMatchFactory::get().add_factory("AnodeMatch",this); }
    /// dtor
    ~AnodeMatchFactory() {}
    /// creation method
    BaseTouchMatch* create(const std::string instance_name) { return new AnodeMatch(instance_name); }
  };

}
#endif
/** @} */ // end of doxygen group
