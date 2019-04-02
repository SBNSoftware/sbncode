#ifndef __sbnanalysis_ProcessorCreator__
#define __sbnanalysis_ProcessorCreator__

#include <string>
#include "sbncode/sbnanalysis/core/Processor/ProcessorBase.hh"

namespace sbnanalysis {

class ProcessorCreator {
public:
  ProcessorCreator(const std::string& name);
  virtual ProcessorBase* Create() = 0;
};


template <class T>
class ProcessorImpl : public ProcessorCreator {
public:
  ProcessorImpl<T>(const std::string& name) : ProcessorCreator(name) {}
  virtual ProcessorBase* Create() { return new T; }
};


#define DECLARE_PROCESSOR(proc) \
  private: \
    static const ProcessorImpl<proc> creator;


#define REGISTER_PROCESSOR(proc) \
  const ProcessorImpl<proc> proc::creator(#proc);

}  // namespace sbnanalysis

#endif  // __sbnanalysis_ProcessorCreator__

