////////////////////////////////////////////////////////////////////////
// Class:       PMTFlashTriggerFilter
// Plugin Type: filter (art v3_03_01)
// File:        PMTFlashTriggerFilter_module.cc
//
// Generated at Fri Apr  3 11:02:06 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include "sbncode/LArRecoProducer/Products/FlashTriggerPrimitive.hh"

#include <memory>
#include <iostream>
#include <vector>

namespace sbn {
  class PMTFlashTriggerFilter;
  bool HasTrigger(const std::vector<FlashTriggerPrimitive> &primitives, int threshold, unsigned n_above_threshold);
}


class sbn::PMTFlashTriggerFilter : public art::EDFilter {
public:
  explicit PMTFlashTriggerFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTFlashTriggerFilter(PMTFlashTriggerFilter const&) = delete;
  PMTFlashTriggerFilter(PMTFlashTriggerFilter&&) = delete;
  PMTFlashTriggerFilter& operator=(PMTFlashTriggerFilter const&) = delete;
  PMTFlashTriggerFilter& operator=(PMTFlashTriggerFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  art::InputTag fFlashTriggerPrimitiveLabel;
  unsigned fNPMTAboveThreshold;
  int fPMTTriggerThreshold;
  bool fStoreDataProduct;
};


sbn::PMTFlashTriggerFilter::PMTFlashTriggerFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fFlashTriggerPrimitiveLabel(p.get<std::string>("FlashTriggerPrimitiveLabel")),
    fNPMTAboveThreshold(p.get<unsigned>("NPMTAboveThreshold")),
    fPMTTriggerThreshold(p.get<int>("PMTTriggerThreshold")),
    fStoreDataProduct(p.get<bool>("StoreDataProduct", false))
  // More initializers here.
{

  if (fStoreDataProduct) produces<bool>();
}

bool sbn::PMTFlashTriggerFilter::filter(art::Event& e)
{
  art::Handle<std::vector<sbn::FlashTriggerPrimitive>> flashtrig_handle;
  e.getByLabel(fFlashTriggerPrimitiveLabel, flashtrig_handle);

  bool ret = false;
  if (flashtrig_handle.isValid()) {
    ret = sbn::HasTrigger(*flashtrig_handle, fPMTTriggerThreshold, fNPMTAboveThreshold);
  }

  if (fStoreDataProduct) {
    std::unique_ptr<bool> store(new bool);
    *store = ret;
    e.put(std::move(store));
    return true;
  }
  return ret;
}

bool sbn::HasTrigger(const std::vector<sbn::FlashTriggerPrimitive> &primitives, int threshold, unsigned n_above_threshold) {
  if (n_above_threshold == 0) return true;

  std::map<int, std::vector<unsigned>> above_threshold;

  for (const sbn::FlashTriggerPrimitive &primitive: primitives) {
    for (const sbn::FlashTriggerPrimitive::Trig &trig: primitive.triggers) {
      if (trig.adc <= threshold) {
        above_threshold[trig.tdc].push_back(primitive.channel);
      }
    }
  }

  for (auto const &pair: above_threshold) {
    if (pair.second.size() >= n_above_threshold) {
      return true;
    }
  }

  return false;
}


DEFINE_ART_MODULE(sbn::PMTFlashTriggerFilter)
