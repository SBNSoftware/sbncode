/**
 * \file BaseTouchMatch.h
 *
 * \ingroup Base
 *
 * \brief Class def header for a class BaseTouchMatch
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef OPT0FINDER_BaseTouchMATCH_H
#define OPT0FINDER_BaseTouchMATCH_H

#include "BaseAlgorithm.h"
namespace flashmatch {

  class FlashMatchManager;

  /**
     \class BaseTouchMatch
     Algorithm base class for matching flashmatch::QCluster_t (TPC object) and \n
     flashmatch::Flash_t (flash). It creates flashmatch::FlashMatch_t which contains \n
     matching infomration.
  */
  class BaseTouchMatch : public BaseAlgorithm{
    friend class FlashMatchManager;

  public:

    /// Default constructor
    BaseTouchMatch(const std::string name="noname") : BaseAlgorithm(kFlashMatch,name)
    {}

    /// Default destructor
    virtual ~BaseTouchMatch(){}

    /// Inspect geometry+timing for T0 tagged tracks
    virtual void Match(const QCluster_t&, const Flash_t&, FlashMatch_t& match) = 0;

  };
}

#endif
/** @} */ // end of doxygen group
