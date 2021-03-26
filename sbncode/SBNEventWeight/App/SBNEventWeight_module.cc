////////////////////////////////////////////////////////////////////////
// Class: SBNEventWeight
// Module Type: producer
// File: SBNEventWeight_module.cc
//
// Generated at Fri Mar 20 09:36:11 2015 by Zarko Pavlovic using artmod
// from cetpkgsupport v1_08_04.
//
// Ported from uboonecode to larsim on Feb 14 2018 by Marco Del Tutto
//
// Ported from larsim to sbncode on Dec 22 2020 by A. Mastbaum
////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "fhiclcpp/ParameterSet.h"

#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightParameterSet.h"
#include "sbncode/SBNEventWeight/Base/WeightManager.h"

namespace sbn {
  namespace evwgh {

class SBNEventWeight : public art::EDProducer {
public:
  explicit SBNEventWeight(fhicl::ParameterSet const& p);

  SBNEventWeight(SBNEventWeight const &) = delete;
  SBNEventWeight(SBNEventWeight &&) = delete;
  SBNEventWeight& operator = (SBNEventWeight const&) = delete;
  SBNEventWeight& operator = (SBNEventWeight&&) = delete;

private:
  void produce(art::Event& e) override;
  void beginRun(art::Run& run) override;

private:
  WeightManager fWeightManager;
  std::string fGenieModuleLabel;
};


SBNEventWeight::SBNEventWeight(fhicl::ParameterSet const& p) : EDProducer{p} {
  fGenieModuleLabel = p.get<std::string>("genie_module_label", "generator");

  const size_t n_func = fWeightManager.Configure(p, *this);
  if (n_func > 0) {
    produces<std::vector<sbn::evwgh::EventWeightMap> >();
    produces<art::Assns<simb::MCTruth, sbn::evwgh::EventWeightMap> >();
    produces<std::vector<sbn::evwgh::EventWeightParameterSet>, art::InRun>();
  }
}


void SBNEventWeight::produce(art::Event& e) {
  auto mcwghvec = std::make_unique<std::vector<EventWeightMap> >();
  auto wghassns = std::make_unique<art::Assns<simb::MCTruth, sbn::evwgh::EventWeightMap> >();

  art::PtrMaker<sbn::evwgh::EventWeightMap> makeWeightPtr(e);

  // Get the MC generator information
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  auto const mcTruthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fGenieModuleLabel);
  art::fill_ptr_vector(mclist, mcTruthHandle);

  // Loop over all truth objects (e.g. neutrinos) in this event
  for (size_t i=0; i<mclist.size(); i++) {
    const EventWeightMap mcwgh = fWeightManager.Run(e, i);
    mcwghvec->push_back(std::move(mcwgh));

    art::Ptr<sbn::evwgh::EventWeightMap> wghPtr = makeWeightPtr(mcwghvec->size() - 1);
    wghassns->addSingle(mclist.at(i), wghPtr);
  }

  e.put(std::move(mcwghvec));
  e.put(std::move(wghassns));
}


void SBNEventWeight::beginRun(art::Run& run) {
  auto p = std::make_unique<std::vector<EventWeightParameterSet> >();
  
  for (auto const& it : fWeightManager.GetWeightCalcMap()) {
    p->push_back(it.second->fParameterSet);
  }

  run.put(std::move(p));
}

  }  // namespace evwgh
}  // namespace sbn

DEFINE_ART_MODULE(sbn::evwgh::SBNEventWeight)

