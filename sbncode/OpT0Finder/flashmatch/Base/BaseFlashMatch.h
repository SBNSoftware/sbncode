/**
 * \file BaseFlashMatch.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseFlashMatch
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef OPT0FINDER_BASEFLASHMATCH_H
#define OPT0FINDER_BASEFLASHMATCH_H

#include "BaseAlgorithm.h"
#include "BaseFlashHypothesis.h"
namespace flashmatch {

  class FlashMatchManager;

  /**
     \class BaseFlashMatch
     Algorithm base class for matching flashmatch::QCluster_t (TPC object) and \n
     flashmatch::Flash_t (flash). It creates flashmatch::FlashMatch_t which contains \n
     matching infomration.
  */
  class BaseFlashMatch : public BaseAlgorithm{
    friend class FlashMatchManager;
    
  public:
    
    /// Default constructor
    BaseFlashMatch(const std::string name="noname") : BaseAlgorithm(kFlashMatch,name)
    {}
    
    /// Default destructor
    virtual ~BaseFlashMatch(){}

    /**
       CORE FUNCTION: takes in flashmatch::QCluster_t (TPC object) and flashmatch::Flash_t (flash) \n
       and inspect if two are consistent (i.t. matched) or not. Returns flashmatch::FlashMatch_t \n
       which represents the compatibility of two inputs. In particular the algorithm needs to  \n
       set the "score" and "QPoint_t" of the return object. The former represents the goodness \n
       of a match with a value larger than 0. Negative value is interpreted as no match.       \n
       The latter represents the matched TPC 3D point. The "tpc_id and "flash_id" of the return \n
       object is set by flashmatch::FlashMatchManager, the caller of the algorithm, as it manages \n
       the overall collection of user input flash and TPC objects. \n
       \n
       Note it is flashmatch::FlashMatchManager configuration option to allow an assignment of the \n
       same flash to multiple TPC object or not. If not allowed, a match with a higher "score"  \n
       in the return object is chosen.
     */
    virtual FlashMatch_t Match(const QCluster_t&, const Flash_t&) = 0;

    /// Method to call flash hypothesis
    Flash_t GetEstimate(const QCluster_t&) const;

    /// Method to simply fill provided reference of flashmatch::Flash_t
    void FillEstimate(const QCluster_t&, Flash_t&) const;

    /// Sets the TPC and Cryo numbers
    virtual void SetTPCCryo(int tpc, int cryo) = 0;

  private:

    void SetFlashHypothesis(flashmatch::BaseFlashHypothesis*);

    flashmatch::BaseFlashHypothesis* _flash_hypothesis;

  protected:

    int _tpc = 0; ///< The TPC number to use
    int _cryo = 0; ///< The Cryostat number to use

  };
}

#endif
/** @} */ // end of doxygen group

