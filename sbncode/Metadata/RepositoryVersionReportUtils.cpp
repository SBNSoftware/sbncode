/**
 * @file   sbncode/Metadata/RepositoryVersionReportUtils.cxx
 * @brief  Interface for _art_ tool reporting the version of packages.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 18, 2025
 * @see    sbncode/Metadata/RepositoryVersionReportUtils.h
 */

// library header
#include "sbncode/Metadata/RepositoryVersionReportUtils.h"

// SBN libraries
#include "sbncode/Metadata/RepositoryVersionReporter.h"

// framework libraries
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard library
#include <memory> // std::unique_ptr<>


// -----------------------------------------------------------------------------
bool sbn::addVersionFromRepository
  (sbn::OrderedPairList& versionList, std::string const& repoName)
{

  std::string const toolName = repoName + "RepositoryVersion";
    
  fhicl::ParameterSet config;
  config.put("tool_type", toolName);
  
  std::unique_ptr<sbn::RepositoryVersionReporter> reportTool;
  try {
    reportTool = art::make_tool<sbn::RepositoryVersionReporter>(config);
  }
  catch(cet::exception const& e) {
    mf::LogWarning{ "RepositoryVersionReporter" }
      << "No report from repository '" << repoName << "': tool not found.";
    mf::LogDebug{ "RepositoryVersionReporter" } << "Error:\n" << e;
    return false;
  }
  
  // blindly add everything
  versionList.items.insert
    (versionList.items.end(), reportTool->begin(), reportTool->end());
  
  return true;
  
} // sbn::addVersionFromRepository()


// -----------------------------------------------------------------------------
std::vector<std::string> sbn::addVersionFromRepositories(
  sbn::OrderedPairList& versionList, std::vector<std::string> const& repoNames
) {
  
  std::vector<std::string> missingRepos;
  for (std::string const& repoName: repoNames) {
    if (addVersionFromRepository(versionList, repoName)) continue;
    missingRepos.push_back(repoName);
  }
  
  return repoNames;
  
} // sbn::addVersionFromRepositories()


// -----------------------------------------------------------------------------
