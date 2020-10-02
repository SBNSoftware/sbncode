/**
 * \class CRTData
 *
 * \ingroup crt
 *
 * \brief CRT Hit Info
 *
 * \author $Author: David Lorca $
 *
 */

#ifndef FlashTriggerPrimitive_SBN_hh_
#define FlashTriggerPrimitive_SBN_hh_

#include <cstdint>
#include <vector>
#include <map>
#include <string>
#include <utility>

namespace sbn {
  struct FlashTriggerPrimitive{
    unsigned channel;
    struct Trig {
      int adc;
      int tdc;
    };
    std::vector<Trig> triggers;
  };

} // namespace sbnd

#endif
