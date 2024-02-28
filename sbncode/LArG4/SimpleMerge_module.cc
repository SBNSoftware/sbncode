////////////////////////////////////////////////////////////////////////
// Class:       SimpleMerge
// Plugin Type: producer (Unknown Unknown)
// File:        producer_module.cc
//
// Generated at Tue Dec 19 08:26:07 2023 by Brinden Carlson using cetskelgen
// from  version .
//
// Producer module for merging multiple input sources into a single output
// without trackID offsets or other complications from merge sim sources.

////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
//#include <utility>

#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcorealg/CoreUtils/enumerate.h"

namespace sbn {
  class SimpleMerge;
}


class sbn::SimpleMerge : public art::EDProducer {
public:
  explicit SimpleMerge(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleMerge(SimpleMerge const&) = delete;
  SimpleMerge(SimpleMerge&&) = delete;
  SimpleMerge& operator=(SimpleMerge const&) = delete;
  SimpleMerge& operator=(SimpleMerge&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::vector<art::InputTag> const fInputSourcesLabels;
  std::vector<int> fResetMotherIDs;
  bool fResetMotherID;

};


sbn::SimpleMerge::SimpleMerge(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fInputSourcesLabels(p.get<std::vector<art::InputTag>>("InputSourcesLabels"))
  , fResetMotherIDs(p.get<std::vector<int>>("ResetMotherIDs", std::vector<int>{}))
  // Anything else you want to merge can be added here...
{
  produces<std::vector<simb::MCParticle>>();
  for (art::InputTag const& tag : fInputSourcesLabels) {
    consumes<std::vector<simb::MCParticle>>(tag);
  }
  fResetMotherID = fResetMotherIDs.size() > 0;
}

void sbn::SimpleMerge::produce(art::Event& e)
{
  auto partCol = std::make_unique<std::vector<simb::MCParticle>>();
  for (auto const& [i_source,input_label] : util::enumerate(fInputSourcesLabels)) {
    auto const& parts = e.getProduct<std::vector<simb::MCParticle>>(input_label);
    for (auto part : parts) {
      if (fResetMotherID){
        for (auto const& offset : fResetMotherIDs) {
          if (part.Mother() == offset) {
            part.SetMother(0);
            break;
          } // if part.Mother
        } // for offsets
      } // if fResetMotherID
      partCol->push_back(part);
    }// for parts
  }// for input_label
  e.put(std::move(partCol));
}// produce

DEFINE_ART_MODULE(sbn::SimpleMerge)
