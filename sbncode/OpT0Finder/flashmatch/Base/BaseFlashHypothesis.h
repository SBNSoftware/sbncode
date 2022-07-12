/**
 * \file FlashHypothesis.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class FlashHypothesis
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef BASEFLASHHYPOTHESIS_H
#define BASEFLASHHYPOTHESIS_H

#include "BaseAlgorithm.h"

namespace flashmatch {
  /**
     \class FlashHypothesis
     User defined class FlashHypothesis ... these comments are used to generate
     doxygen documentation!
  */
  class BaseFlashHypothesis : public flashmatch::BaseAlgorithm {
    
  public:

    /// Default constructor
    BaseFlashHypothesis(const std::string name="noname");
    
    /// Default destructor
    ~BaseFlashHypothesis(){}

    /// Method to create flashmatch::Flash_t object and return
    Flash_t GetEstimate(const QCluster_t&) const;

    /// Method to simply fill provided reference of flashmatch::Flash_t
    virtual void FillEstimate(const QCluster_t&, Flash_t&) const = 0;

    /// Sets the channels to use
    void SetChannelMask(std::vector<size_t> ch_mask);

    /// Sets the channels sensitive to visible light
    void SetUncoatedPMTs(std::vector<size_t> ch_uncoated);

  protected:

    std::vector<bool> _channel_mask; ///< The list of channels to use
    std::vector<bool> _uncoated_pmt_list; ///< A list of opdet sensitive to visible (reflected) light

  };
}
#endif
/** @} */ // end of doxygen group 

