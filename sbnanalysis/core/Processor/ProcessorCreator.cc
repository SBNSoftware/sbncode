#include "ProcessorCreator.hh"
#include "ProcessorFactory.hh"

namespace sbnanalysis {

ProcessorCreator::ProcessorCreator(const std::string& name) {
  ProcessorFactory::Register(name, this);
}

}  // namespace sbnanalysis

