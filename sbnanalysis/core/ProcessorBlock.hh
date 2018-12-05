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

namespace fhicl {
  class ParameterSet;
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
   * Note that the ProcessorBlock takes ownership of the Processor.
   *
   * \param processor The processor
   * \param config The configuration, if any
   */
  virtual void AddProcessor(ProcessorBase* processor, fhicl::ParameterSet* config);

  /**
   * Process a set of files.
   *
   * \param filenames A list of art ROOT files to process
   */
  virtual void ProcessFiles(std::vector<std::string> filenames);

  /** Delete all processors owned by the block. */
  virtual void DeleteProcessors();

protected:
  /** Processors and their configurations. */
  std::vector<std::pair<ProcessorBase*, fhicl::ParameterSet*> > fProcessors;
};

}  // namespace core

#endif  // __sbnanalysis_core_ProcessorBlock__

