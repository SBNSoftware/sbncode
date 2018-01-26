#ifndef __sbnanalysis_core_ProcessorBase__
#define __sbnanalysis_core_ProcessorBase__

/**
 * \file ProcessorBase.hh
 *
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

/** Core framework functionality. */
namespace core {

/**
 * \class core::ProcessorBase
 * \brief A generic tree-writing event-by-event processor.
 */
class ProcessorBase {
public:
  /** Constructor */
  ProcessorBase();

  /** Destructor */
  virtual ~ProcessorBase();

  /**
   * Process a set of files.
   *
   * \param filenames A list of art ROOT files to process
   * \param config A configuration as a JSON object
   */
  virtual void ProcessFiles(std::vector<std::string> filenames,
                            Json::Value* config=NULL);

  /**
   * Add a branch to the output tree.
   *
   * Called in user subclasses to augment the default event tree.
   * This mirrors the TTree::Branch API.
   *
   * \param name The branch name
   * \param obj A pointer to the object
   * \returns A pointer to the created TBranch (we retain ownership)
   */
  template<class T>
  TBranch* AddBranch(std::string name, T* obj) {
    return fTree->Branch(name.c_str(), obj);
  }

  /** User-defined functions. */

  /**
   * Process one event.
   *
   * \param ev The event, as a gallery::Event
   */
  virtual void ProcessEvent(gallery::Event& ev) = 0;

  /**
   * Perform user-level initialization.
   *
   * \param config A configuration, as a JSON object.
   */
  virtual void Initialize(Json::Value* config=NULL) = 0;

  /** Perform user-level finalization. */
  virtual void Finalize() = 0;

protected:
  /**
   * Perform framework-level initialization.
   *
   * \param config A configuration as a JSON object
   */
  virtual void Setup(Json::Value* config=NULL);

  /** Perform framework-level finalization. */
  virtual void Teardown();

  /**
   * Populate the default event tree variables.
   *
   * \param ev The current gallery event
  */
  void FillEventTree(gallery::Event& ev);

  unsigned long fEventIndex;  //!< An incrementing index
  std::string fOutputFilename;  //!< The output filename
  TFile* fOutputFile;  //!< The output ROOT file
  TTree* fTree;  //!< The output ROOT tree
  Event* fEvent;  //!< The standard output event data structure
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
};

}  // namespace core


/** Macro to create plugin library for user-defined Processors. */
#define DECLARE_SBN_PROCESSOR(classname) \
extern "C" classname* CreateObject() { return new classname; } \
extern "C" void DestroyObject(classname* o) { delete o; }

#endif  // __sbnanalysis_core_ProcessorBase__

