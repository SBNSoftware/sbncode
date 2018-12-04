#include <string>
#include <vector>
#include <json/json.h>
#include "ProcessorBase.hh"
#include "ProcessorBlock.hh"

namespace core {

ProcessorBlock::ProcessorBlock() {}


ProcessorBlock::~ProcessorBlock() {}


void ProcessorBlock::AddProcessor(ProcessorBase* processor,
                                  Json::Value* config) {
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

