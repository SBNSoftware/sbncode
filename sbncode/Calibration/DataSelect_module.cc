////////////////////////////////////////////////////////////////////////
// Class:       DataSelect
// Plugin Type: filter (art v3_05_01)
// File:        DataSelect_module.cc
//
// Module for selecting events based on the presence of data products.
//
// Useful for processing data events, which sometimes have missing 
// data products.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Hit.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

namespace sbn {
  namespace util {
    class DataSelect;
  }
}


class sbn::util::DataSelect : public art::EDFilter {
public:
  explicit DataSelect(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DataSelect(DataSelect const&) = delete;
  DataSelect(DataSelect&&) = delete;
  DataSelect& operator=(DataSelect const&) = delete;
  DataSelect& operator=(DataSelect&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  std::string fHitLabel;
};


sbn::util::DataSelect::DataSelect(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fHitLabel(p.get<std::string>("HitLabel"))
{}

bool sbn::util::DataSelect::filter(art::Event& e)
{

  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitLabel, hitHandle);

  return hitHandle.isValid();
}

DEFINE_ART_MODULE(sbn::util::DataSelect)
