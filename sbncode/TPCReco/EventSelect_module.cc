////////////////////////////////////////////////////////////////////////
// Class:       EventSelect
// Plugin Type: filter (art v3_05_01)
// File:        EventSelect_module.cc
//
// Generated at Mon Oct  5 13:03:33 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_10_00.
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

namespace sbn {
  namespace util {
    class EventSelect;
  }
}


class sbn::util::EventSelect : public art::EDFilter {
public:
  explicit EventSelect(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EventSelect(EventSelect const&) = delete;
  EventSelect(EventSelect&&) = delete;
  EventSelect& operator=(EventSelect const&) = delete;
  EventSelect& operator=(EventSelect&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  void respondToOpenInputFile(const art::FileBlock& fb) { fFileNo ++; }

private:

  // Declare member data here.
  std::vector<std::string> fConfig;
  std::vector<int> fFileNoSelect;
  std::vector<int> fEvtSelect;
  int fFileNo;
  unsigned fNAccept;
};


sbn::util::EventSelect::EventSelect(fhicl::ParameterSet const& p)
  : EDFilter{p},
    fConfig(p.get<std::vector<std::string>>("Select"))
{
  for (const std::string &s: fConfig) {
    fFileNoSelect.push_back(-1);
    fEvtSelect.push_back(-1);
    sscanf(s.c_str(), "%d:%d", &fFileNoSelect.back(), &fEvtSelect.back());
  }
  fFileNo = -1;
  fNAccept = 0;
}

bool sbn::util::EventSelect::filter(art::Event& e)
{
  //if (fNAccept == fFileNoSelect.size()) {
  //  throw cet::exception("EventSelect Filter") << "Finished selecting events from list. ";
  //}

  for (unsigned i = 0; i < fFileNoSelect.size(); i++) {
    if (fFileNoSelect[i] == fFileNo && fEvtSelect[i] == (int)e.event()) {
      fNAccept ++;
      return true;
    }
  }
  return false;
}

DEFINE_ART_MODULE(sbn::util::EventSelect)
