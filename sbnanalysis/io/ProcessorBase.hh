#ifndef __sbnanalysis_io_ProcessorBase__
#define __sbnanalysis_io_ProcessorBase__

/**
 * A generic processor that writes an sbnanalysis tree.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/25
 */

#include <string>
#include <vector>
#include "gallery/Event.h"

class TFile;
class TTree;
class Event;

namespace io {

class ProcessorBase {
public:
  ProcessorBase();

  virtual ~ProcessorBase();

  virtual void ProcessFile(std::vector<std::string> filenames);

  virtual void Initialize();

  virtual void Finalize();

  virtual void AddBranch() {}

  virtual void ProcessEvent(gallery::Event& ev) = 0;

protected:
  void FillEventTree(gallery::Event& ev);

  unsigned long fEventIndex;
  std::string fOutputFilename;
  TFile* fOutputFile;
  TTree* fTree;
  Event* fEvent;
};

}  // namespace io

#endif  // __sbnanalysis_io_ProcessorBase__

