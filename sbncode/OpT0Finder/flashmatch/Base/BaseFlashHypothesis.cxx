#ifndef BASEFLASHHYPOTHESIS_CXX
#define BASEFLASHHYPOTHESIS_CXX

#include "BaseFlashHypothesis.h"
#include "OpT0FinderException.h"
namespace flashmatch {

  BaseFlashHypothesis::BaseFlashHypothesis(std::string name)
    : BaseAlgorithm(kFlashHypothesis,name)
    , _channel_mask(DetectorSpecs::GetME().NOpDets(),false)
    , _uncoated_pmt_list(DetectorSpecs::GetME().NOpDets(),false)
  {}


  Flash_t BaseFlashHypothesis::GetEstimate(const QCluster_t& tpc) const
  {
    Flash_t res;
    //res.pe_v.resize(OpDetXArray().size());

    FillEstimate(tpc,res);
    return res;
  }

  void BaseFlashHypothesis::SetChannelMask(std::vector<size_t> ch_mask)
  {
    for(auto const& v : ch_mask) {
      if(v >= _channel_mask.size()) {
        FLASH_CRITICAL() << "ChannelMask array contains ID " << v
                         << " >= number of opdets (" << DetectorSpecs::GetME().NOpDets() << ")!" << std::endl;
        throw OpT0FinderException();
      }
      _channel_mask[v] = true;
    }
  }

  void BaseFlashHypothesis::SetUncoatedPMTs(std::vector<size_t> ch_uncoated)
  {
    for(auto const& v : ch_uncoated) {
      if(v >= _uncoated_pmt_list.size()) {
        FLASH_CRITICAL() << "Uncoated PMT List array contains ID " << v
                         << " >= number of opdets (" << DetectorSpecs::GetME().NOpDets() << ")!" << std::endl;
        throw OpT0FinderException();
      }
      _uncoated_pmt_list[v] = true;
    }
  }

}
#endif
