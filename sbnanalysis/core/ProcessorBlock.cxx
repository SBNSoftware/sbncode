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
      it.first->BuildEventTree(ev);

      bool accept = it.first->ProcessEvent(ev, it.first->fEvent->truth, *it.first->fReco);

      if (accept) {
        it.first->FillTree();

        // For each reco event, fill the reco output tree
        it.first->fRecoEvent->experiment = it.first->fEvent->experiment;
        it.first->fRecoEvent->metadata = it.first->fEvent->metadata;
        for (auto const& reco : *it.first->fReco) {
          it.first->fRecoEvent->reco = reco;
          it.first->FillRecoTree();
        }
      }

      it.first->EventCleanup();
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

