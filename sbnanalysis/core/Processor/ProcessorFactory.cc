#include "ProcessorFactory.hh"
#include "ProcessorCreator.hh"

namespace sbnanalysis {

ProcessorBase* ProcessorFactory::Create(const std::string& name) {
  std::map<std::string, ProcessorCreator*>::iterator it = \
    GetTable().find(name);

  if (it == GetTable().end()) {
    return nullptr;
  }

  return it->second->Create();
}

void ProcessorFactory::Register(const std::string& name,
                                ProcessorCreator* creator) {
  GetTable()[name] = creator;
}

std::map<std::string, ProcessorCreator*>& ProcessorFactory::GetTable() {
  static std::map<std::string, ProcessorCreator*> table;
  return table;
}

}  // namespace sbnanalysis

