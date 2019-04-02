////////////////////////////////////////////////////////////////////////
// Class:       SBNAnalysis
// Plugin Type: analyzer (art v3_01_02)
// File:        SBNAnalysis_module.cc
//
// Generated at Tue Apr  2 04:48:23 2019 by Andrew Mastbaum using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "sbncode/sbnanalysis/core/Processor/ProcessorBase.hh"
#include "sbncode/sbnanalysis/core/Processor/ProcessorFactory.hh"

class SBNAnalysis : public art::EDAnalyzer {
public:
  explicit SBNAnalysis(fhicl::ParameterSet const& p);

  SBNAnalysis(SBNAnalysis const&) = delete;
  SBNAnalysis(SBNAnalysis&&) = delete;
  SBNAnalysis& operator=(SBNAnalysis const&) = delete;
  SBNAnalysis& operator=(SBNAnalysis&&) = delete;

  void beginJob() override;
  void analyze(art::Event const& e) override;
  void endJob() override;

protected:
  /** Processors and their configurations. */
  std::vector<std::pair<sbnanalysis::ProcessorBase*, fhicl::ParameterSet> > fProcessors;
};


SBNAnalysis::SBNAnalysis(fhicl::ParameterSet const& p) : EDAnalyzer{p} {
  std::string label = p.get<std::string>("module_label");
  std::vector<std::string> procnames = p.get<std::vector<std::string> >("processors");

  for (auto const& procname : procnames) {
    fhicl::ParameterSet pset = p.get<fhicl::ParameterSet>(procname);
    std::string type = pset.get<std::string>("type");
    sbnanalysis::ProcessorBase* proc = sbnanalysis::ProcessorFactory::Create(procname);
    if (!proc) {
      throw cet::exception(__FUNCTION__)
        << "Unknown processor " << procname << std::endl;
    }

    fProcessors.push_back({ proc, pset });
  }
}


void SBNAnalysis::beginJob() {
  for (auto it : fProcessors) {
    it.first->Setup(&it.second);
    it.first->Initialize(&it.second);
  }
}


void SBNAnalysis::analyze(art::Event const& ev) {
  for (auto it : fProcessors) {
    it.first->BuildEventTree(ev);

    bool accept = \
      it.first->ProcessEvent(ev, it.first->GetEvent()->truth, *it.first->fReco);

    if (accept) {
      it.first->FillTree();
    }

    it.first->EventCleanup();
  }
}


void SBNAnalysis::endJob() {
  for (auto it : fProcessors) {
    it.first->Finalize();
    it.first->Teardown();
  }

  for (auto it : fProcessors) {
    delete it.first;
  }
}

DEFINE_ART_MODULE(SBNAnalysis)

