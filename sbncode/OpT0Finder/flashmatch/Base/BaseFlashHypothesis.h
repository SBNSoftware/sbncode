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

#if USING_LARSOFT == 1
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#endif

namespace flashmatch {
  /**
     \class FlashHypothesis
     User defined class FlashHypothesis ... these comments are used to generate
     doxygen documentation!
  */
  class BaseFlashHypothesis : public flashmatch::BaseAlgorithm {
    
  public:

    /// Default constructor
    BaseFlashHypothesis(const std::string name="noname")
      : flashmatch::BaseAlgorithm(flashmatch::kFlashHypothesis,name)
    {}
    
    /// Default destructor
    ~BaseFlashHypothesis(){}

    /// Method to create flashmatch::Flash_t object and return
    Flash_t GetEstimate(const QCluster_t&) const;

    /// Method to simply fill provided reference of flashmatch::Flash_t
    virtual void FillEstimate(const QCluster_t&, Flash_t&) const = 0;

    /// Sets the channels to use
    void SetChannelMask(std::vector<int> ch_mask) { _channel_mask = ch_mask; }

    /// Sets the channel type (pmt vs. xarapuca)
    void SetChannelType(std::vector<int> ch_type) { _channel_type = ch_type; }

    /// Sets the channels sensitive to visible light
    void SetUncoatedPMTs(std::vector<int> ch_uncoated) { _uncoated_pmt_list = ch_uncoated; }

    #if USING_LARSOFT == 1
    /// Sets the semi analytical model
    void SetSemiAnalyticalModel(std::unique_ptr<phot::SemiAnalyticalModel> model) { _semi_model = std::move(model); }
    #endif

  protected:

    std::vector<int> _channel_mask; ///< The list of channels to use
    std::vector<int> _uncoated_pmt_list; ///< A list of opdet sensitive to visible (reflected) light
    std::vector<int> _channel_type; 
    // std::vector<int> _visible_ch_list; 

    #if USING_LARSOFT == 1
    std::unique_ptr<phot::SemiAnalyticalModel> _semi_model;
    #endif

  };
}
#endif
/** @} */ // end of doxygen group 

