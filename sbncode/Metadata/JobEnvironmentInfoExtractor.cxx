/**
 * @file   sbncode/Metadata/JobEnvironmentInfoExtractor.cxx
 * @brief  Algorithm extracting information from the job execution environment.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 16, 2025
 * @see    sbncode/Metadata/JobEnvironmentInfoExtractor.h
 */

// library header
#include "sbncode/Metadata/JobEnvironmentInfoExtractor.h"

// SBN libraries
#include "sbncode/Metadata/RepositoryVersionReportUtils.h"
#include "sbncode/Metadata/RepositoryVersionReporter.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard library
#include <memory>
#include <string>
#include <utility> // std::move()

// system libraries
#include <string.h> // strchr()
#include <unistd.h> // environ


// -----------------------------------------------------------------------------
std::vector<std::string> const
sbn::JobEnvironmentInfoExtractor::Config::DefaultDropVars = {
    R"|(.*%%$)|" // BASH function pattern
  , R"|(.*\(\)$)|" // another BASH function pattern
};



// -----------------------------------------------------------------------------
sbn::JobEnvironmentInfoExtractor::JobEnvironmentInfoExtractor
  (Config const& config /* = Config{} */)
  : fDropVars{ prepareRegEx(config.dropVars) }
  , fRepositoryReports{ config.repositories }
  , fLogCategory{ config.logCategory }
{
  
}


// -----------------------------------------------------------------------------
sbn::JobEnvironmentInfo sbn::JobEnvironmentInfoExtractor::extract
  (art::ModuleDescription const& moduleInfo)
{
  
  sbn::JobEnvironmentInfo info;
  
  info.processName = moduleInfo.processName();
  info.artVersion = moduleInfo.releaseVersion();
  
  info.sources = extractSourceVersions();
  
  info.variables = extractEnvironmentVariables();
  
  return info;
  
} // sbn::JobEnvironmentInfoExtractor::extract()


// -----------------------------------------------------------------------------
auto sbn::JobEnvironmentInfoExtractor::extractSourceVersions() const
  -> sbn::OrderedPairList
{
  
  // take advantage of the basic 
  sbn::OrderedPairList sourceVersions;
  addVersionFromRepositories(sourceVersions, fRepositoryReports);
  sourceVersions.finish();
  
  return sourceVersions;
  
} // sbn::JobEnvironmentInfoExtractor::extractSourceVersions()


// -----------------------------------------------------------------------------
auto sbn::JobEnvironmentInfoExtractor::extractEnvironmentVariables() const
  -> sbn::OrderedPairList
{
  
  sbn::OrderedPairList vars;
  
  char** itemPtr = environ;
  while (const char* item = *(itemPtr++)) {
    
    const char* sep = strchr(item, '=');
    
    std::string name, value;
    if (sep) {
      name = std::string{ item, sep++ };
      value = std::string{ sep };
    }
    else {
      mf::LogDebug{ fLogCategory }
        << "Environment element '" << item << "' is not in name=value form.";
      name = item;
    }
    
    // apply drop filters
    if (matchPatterns(name, fDropVars) < fDropVars.size())
      continue;
    
    vars.items.emplace_back(std::move(name), std::move(value));
    
  } // while
  
  vars.finish();
  
  return vars;
  
} // sbn::JobEnvironmentInfoExtractor::extractEnvironmentVariables()


// -----------------------------------------------------------------------------
std::size_t sbn::JobEnvironmentInfoExtractor::matchPatterns
  (std::string const& value, std::vector<std::regex> const& patterns)
{
  for (auto const& [ iPattern, pattern ]: util::enumerate(patterns))
    if (std::regex_match(value, pattern)) return iPattern;
  return patterns.size();
}


// -----------------------------------------------------------------------------
std::vector<std::regex> sbn::JobEnvironmentInfoExtractor::prepareRegEx
  (std::vector<std::string> const& patterns)
{
  std::vector<std::regex> regex;
  regex.reserve(patterns.size());
  for (std::string const& pattern: patterns) regex.emplace_back(pattern);
  return regex;
}


// -----------------------------------------------------------------------------
