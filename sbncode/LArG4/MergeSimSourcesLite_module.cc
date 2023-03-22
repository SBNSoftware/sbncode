////////////////////////////////////////////////////////////////////////
// Class:       MergeSimSourcesLite
// Plugin Type: producer (Unknown Unknown)
// File:        MergeSimSourcesLite_module.cc
//
// Generated at Wed Mar 22 00:14:48 2023 by Laura Mai Hien Domine using cetskelgen
// from  version .
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

//#include "MergeSimSourcesLiteUtility.h"
#include "lardataobj/MCBase/MCParticleLite.h"
//#include "larsim/MergeSimSources/MergeSimSources_module.cc"

#include <memory>

namespace sbn {
  class MergeSimSourcesLite;
}


class sbn::MergeSimSourcesLite : public sim::MergeSimSources {
public:
  explicit MergeSimSourcesLite(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MergeSimSourcesLite(MergeSimSourcesLite const&) = delete;
  MergeSimSourcesLite(MergeSimSourcesLite&&) = delete;
  MergeSimSourcesLite& operator=(MergeSimSourcesLite const&) = delete;
  MergeSimSourcesLite& operator=(MergeSimSourcesLite&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

};


sbn::MergeSimSourcesLite::MergeSimSourcesLite(fhicl::ParameterSet const& p)
  : sim::MergeSimSources{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  for (art::InputTag const& tag : fInputSourcesLabels) {
    if (fFillMCParticles) {
      consumes<std::vector<sim::MCParticleLite>>(tag);
      produces<std::vector<sim::MCParticleLite>>();
    }
  }
}

void sbn::MergeSimSourcesLite::produce(art::Event& e)
{
  // Run the parent's produce
  sim::MergeSimSources::produce(e);

  // Deal with MCParticleLite
  auto partLiteCol = std::make_unique<std::vector<sim::MCParticleLite>>();

  MergeSimSourcesLiteUtility MergeUtility{fTrackIDOffsets};

  for (auto const& [i_source, input_label] : util::enumerate(fInputSourcesLabels)) {
    if (fFillMCParticles) {
      art::PtrMaker<sim::MCParticleLite> const makePartPtr{e};
      auto const input_partLiteCol = e.getValidHandle<sim::MCParticleLite>>(input_label);
      MergeUtility.MergeMCParticleLites(*partLiteCol, *input_partLiteCol, i_source);
    }
  }

  if (fFillMCParticles) {
    e.put(std::move(partLiteCol));
  }
}

DEFINE_ART_MODULE(sbn::MergeSimSourcesLite)
