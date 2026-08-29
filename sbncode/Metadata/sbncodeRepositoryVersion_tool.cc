/**
 * @file   RepositoryVersion_sbncode_tool.cc
 * @brief  _art_ tool reporting the version of `sbncode`-related packages.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 18, 2025
 * 
 */

// SBN libraries
#include "sbncode/Metadata/RepositoryVersionReporter.h"
#include "sbncode/Metadata/RepositoryVersion_sbncode.h"
#include "sbnobj/Metadata/RepositoryVersion_sbnobj.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard libraries
#include <string>
#include <vector>


// -----------------------------------------------------------------------------
namespace sbn { struct sbncodeRepositoryVersion; }
/**
 * @brief Implements the `sbn::RepositoryVersionReporter` interface for
 * `sbncode`.
 * 
 * It collects information from the following repositories: `sbncode` and
 * `sbnobj`.
 * 
 */
struct sbn::sbncodeRepositoryVersion: public sbn::RepositoryVersionReporter {
  
  struct Config {};
  
  using Parameters = art::ToolConfigTable<Config>;
  
  sbncodeRepositoryVersion(Parameters const&);
  
}; // sbn::sbncodeRepositoryVersion()


// -----------------------------------------------------------------------------
// ---  implementation
// -----------------------------------------------------------------------------
sbn::sbncodeRepositoryVersion::sbncodeRepositoryVersion(Parameters const&) {
  
  packageVersions.items.emplace_back("sbnobj", ::RepositoryVersion_sbnobj);
  packageVersions.items.emplace_back("sbncode", ::RepositoryVersion_sbncode);
  
  packageVersions.finish();
  
}


// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(sbn::sbncodeRepositoryVersion)


// -----------------------------------------------------------------------------
