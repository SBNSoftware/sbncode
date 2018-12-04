#ifndef __sbnanalysis_core_Loader__
#define __sbnanalysis_core_Loader__

/**
 * \file Loader.hh
 *
 * Loading Processors from shared libraries.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/02/05
 */

#include <string>
#include <vector>

namespace Json {
  class Value;
}

namespace core {

class ProcessorBase;
class PostProcessorBase;

/**
 * \struct export_table
 * \brief Struct containing (macro defined) creation/deletion operations
 */
struct export_table {
  /** Function to create an Processor instance. */
  ProcessorBase* (*create)(void);

  /** Function to delete an Processor instance. */
  void (*destroy)(ProcessorBase*);
};

/**
 * \struct export_table_postprocess
 * \brief Struct containing (macro defined) creation/deletion operations
 */
struct export_table_postprocess {
  /** Function to create an PostProcessor instance. */
  PostProcessorBase* (*create)(void);

  /** Function to delete an PostProcessor instance. */
  void (*destroy)(PostProcessorBase*);
};

/**
 * Load a processor from a shared library.
 *
 * \param libpath Path to the shared library
 * \returns A table of exported hooks for creating/deleting instances
 */
export_table* LoadProcessor(char* libname);

/**
 * Load a post-processor from a shared library.
 *
 * \param libpath Path to the shared library
 * \returns A table of exported hooks for creating/deleting instances
 */
export_table_postprocess* LoadPostProcessor(char* libname);


/**
 * Load configuration from JSON file to object.
 *
 * \param config Path to the JSON file
 * \returns Configuration as a JSON object
 */
Json::Value* LoadConfig(char* configfile);


/**
 * \class ProcessorDef
 * \brief Definition for a dynamically-loaded Processor
 */
struct ProcessorDef {
  char* name;  //!< Name of the Processor
  struct export_table* exports;  //!< Table for ctor/dtor access
  Json::Value* config;  //!< Configuration as an object
  core::ProcessorBase* proc;  //!< Pointer to the processor
};

}  // namespace core

#endif  // __sbnanalysis_core_Loader__

