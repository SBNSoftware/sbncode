#include <string>
#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "ProcessorBase.hh"
#include "ProcessorBlock.hh"

namespace core {

ProcessorBlock::ProcessorBlock() {}


ProcessorBlock::~ProcessorBlock() {}


void ProcessorBlock::AddProcessor(ProcessorBase* processor,
                                  fhicl::ParameterSet* config) {
  fProcessors.push_back({processor, config});
}


void ProcessorBlock::ProcessFiles(std::vector<std::string> filenames) {
  // Setup
  for (auto it : fProcessors) {
    it.first->Setup(it.second);
    it.first->Initialize(it.second);
  }

  // Event loop
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    for (auto it : fProcessors) {

      ProcessorBase* p = it.first;
      p->BuildEventTree(ev);
      p->SetupServices(ev);

      bool accept = p->ProcessEvent(ev, p->fEvent->truth, *p->fReco);

      // Set reco index
      for (size_t i=0; i<p->fReco->size(); i++) {
        p->fReco->at(i).index = i;
      }
      p->fEvent->nreco = p->fReco->size();

      if (accept) {
        p->FillTree();

        // For each reco event, fill the reco output tree
        p->fRecoEvent->experiment = p->fEvent->experiment;
        p->fRecoEvent->metadata = p->fEvent->metadata;

        for (auto const& reco : *p->fReco) {
          p->fRecoEvent->reco = reco;

          // Copy the associated truth interaction, if any
          p->fRecoEvent->truth.resize(0);
          int truthidx = p->fRecoEvent->reco.truth_index;
          if (truthidx >= 0 && truthidx < p->fEvent->truth.size()) {
            p->fRecoEvent->truth.push_back(p->fEvent->truth.at(truthidx));
          }

          p->FillRecoTree();
        }
      }

      p->EventCleanup();
    }
  }

  // Finalize
  for (auto it : fProcessors) {
    it.first->Finalize();
    it.first->Teardown();
  }
}


void ProcessorBlock::DeleteProcessors() {
  for (auto it : fProcessors) {
    delete it.first;
  }
}

}  // namespace core

