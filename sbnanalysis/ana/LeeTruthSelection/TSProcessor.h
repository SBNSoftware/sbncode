#ifndef SBNANA_TSPROCESSOR_H
#define SBNANA_TSPROCESSOR_H

#include "core/SelectionBase.hh"

namespace ana {
namespace lee_truth_selection {

class TSProcessor: core::SelectionBase {
public:
  TSProcessor() : SelectionBase(), _selections() {}
  
  void Initialize(Json::Value *config);
  void ProcessEvent(gallery::Event& ev);
  void Finalize();

protected:
  std::vector<TSSelection> _selections;
};
}


}
#endif
