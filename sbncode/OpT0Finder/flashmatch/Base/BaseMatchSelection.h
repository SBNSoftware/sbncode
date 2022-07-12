/**
 * \file BaseMatchSelection.h
 *
 * \ingroup Base
 *
 * \brief Class def header for a class BaseMatchSelection
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef OPT0FINDER_BaseMatchSelection_H
#define OPT0FINDER_BaseMatchSelection_H

#include "BaseAlgorithm.h"

namespace flashmatch {

  class FlashMatchManager;

  /**
     \class BaseMatchSelection
     Algorithm base class for matching flashmatch::QCluster_t (TPC object) and \n
     flashmatch::Flash_t (flash). It creates flashmatch::FlashMatch_t which contains \n
     matching infomration.
  */
  class BaseMatchSelection : public BaseAlgorithm{
    friend class FlashMatchManager;

  public:

    /// Default constructor
    BaseMatchSelection(const std::string name="noname") : BaseAlgorithm(kFlashMatch,name)
    {}

    /// Default destructor
    virtual ~BaseMatchSelection(){}

    /**
     * Given all matching data, return the selected correct matches
     */
    virtual std::vector<FlashMatch_t> Select(const std::vector<std::vector<FlashMatch_t> >& match_data) = 0;

  };
}

#endif
/** @} */ // end of doxygen group
