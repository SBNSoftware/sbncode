#ifndef __sbnanalysis_ProcessorFactory__
#define __sbnanalysis_ProcessorFactory__

#include <map>
#include <string>
#include "sbncode/sbnanalysis/core/Processor/ProcessorBase.hh"
#include "sbncode/sbnanalysis/core/Processor/ProcessorCreator.hh"

namespace sbnanalysis {

class ProcessorFactory {
public:
  static ProcessorBase* Create(const std::string& name);
  static void Register(const std::string& name, ProcessorCreator* creator);

private:
  static std::map<std::string, ProcessorCreator*>& GetTable();
};

}  // namespace sbnanalysis

#endif  // __sbnanalysis_ProcessorFactory__

