#ifndef __sbnanalysis_core_ProcessorBase__
#define __sbnanalysis_core_ProcessorBase__

/**
 * A generic processor that writes an sbnanalysis tree.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/25
 */

#include <string>
#include <vector>
#include "gallery/Event.h"

class TBranch;
class TFile;
class TTree;
class Event;

namespace Json {
  class Value;
}

namespace core {

class ProcessorBase {
public:
  ProcessorBase();
  virtual ~ProcessorBase();

  // Core processor base functions
  virtual void Setup(Json::Value* config=NULL);
  virtual void Teardown();

  virtual void ProcessFiles(std::vector<std::string> filenames,
                            Json::Value* config=NULL);

  template<class T>
  TBranch* AddBranch(std::string name, T* obj) {
    return fTree->Branch(name.c_str(), obj);
  }

  // User-defined
  virtual void ProcessEvent(gallery::Event& ev) = 0;
  virtual void Initialize(Json::Value* config=NULL) = 0;
  virtual void Finalize() = 0;

protected:
  void FillEventTree(gallery::Event& ev);

  unsigned long fEventIndex;
  std::string fOutputFilename;
  TFile* fOutputFile;
  TTree* fTree;
  Event* fEvent;
};

}  // namespace core


/** Macro to create plugin library for user-defined Processors. */
#define DECLARE_SBN_PROCESSOR(classname) \
extern "C" classname* CreateObject() { return new classname; } \
extern "C" void DestroyObject(classname* o) { delete o; }

#endif  // __sbnanalysis_core_ProcessorBase__

