#ifndef __sbnanalysis_core_ProcessorBlock__
#define __sbnanalysis_core_ProcessorBlock__

/**
 * \file ProcessorBlock.hh
 *
 * Processor management.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/30
 */

#include <string>
#include <vector>
#include "ProcessorBase.hh"

namespace Json {
  class Value;
}

namespace core {

/**
 * \class core::ProcessorBlock
 * \brief A set of Processors
 */
class ProcessorBlock {
public:
  /** Constructor */
  ProcessorBlock();

  /** Destructor */
  virtual ~ProcessorBlock();

  /**
   * Add a processor to the block.
   *
   * \param processor The processor
   * \param config The configuration, if any
   */
  virtual void AddProcessor(ProcessorBase* processor, Json::Value* config);

  /**
   * Process a set of files.
   *
   * \param filenames A list of art ROOT files to process
   */
  virtual void ProcessFiles(std::vector<std::string> filenames);

protected:
  /** Processors and their configurations. */
  std::vector<std::pair<ProcessorBase*, Json::Value*> > fProcessors;
};

}  // namespace core

#endif  // __sbnanalysis_core_ProcessorBlock__

