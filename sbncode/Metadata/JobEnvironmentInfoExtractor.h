/**
 * @file   sbncode/Metadata/JobEnvironmentInfoExtractor.h
 * @brief  Algorithm extracting information from the job execution environment.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 16, 2025
 * @see    sbncode/Metadata/JobEnvironmentInfoExtractor.cxx
 */

#ifndef SBNCODE_METADATA_JOBENVIRONMENTINFOEXTRACTOR_H
#define SBNCODE_METADATA_JOBENVIRONMENTINFOEXTRACTOR_H


// local libraries
#include "sbnobj/Common/Metadata/JobEnvironmentInfo.h"

// C++ standard libraries
#include <regex>
#include <string>
#include <vector>


// -----------------------------------------------------------------------------
namespace art { class ModuleDescription; } // forward declarations

// -----------------------------------------------------------------------------
namespace sbn { class JobEnvironmentInfoExtractor; }
/**
 * @brief Extracts job execution environment information.
 * 
 * In its simplest form, the information is saved into a
 * `sbn::JobEnvironmentInfo` object:
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * sbn::JobEnvironmentInfoExtractor extractor;
 * sbn::JobEnvironmentInfo const info = extractor(moduleDescription());
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * The `moduleDescription()` is a _art_ module member function, and as such this
 * call is bound to the _art_ framework. However, if the context does not allow
 * for a `art::ModuleDescription` object, partial information can still be
 * retrieved which does not need that object.
 * 
 * 
 * Configuration
 * --------------
 * 
 * The configuration object, `Config`, includes the following configurable
 * elements:
 *  * `dropVars`: a list of regular expression (`std::regex`) patterns to be
 *    matched to the environment variable names; if the name of a variable
 *    matches any of the the full patterns specified here, its value is not
 *    included in the metadata. This is useful to exclude known irrelevant
 *    variables (for example, shell functions).
 *  * `repositories`: a list of names of repositories to query for their
 *    version. For each repository `<reponame>`, an attempt to load a _art_ tool
 *    named `<reponame>RepositoryVersion` with an empty configuration parameter
 *    set is attempted, and on success all the information obtained from the
 *    tool is integrated into the metadata.
 *  * `logCategory`: this algorithm uses messagefacility library to emit
 *    messages to console. This is the category these messages are assigned to.
 * 
 * 
 */
class sbn::JobEnvironmentInfoExtractor {
  
    public:
  
  /// Configuration record.
  struct Config {
    
    /// Default value of `dropVars`.
    static std::vector<std::string> const DefaultDropVars;
    
    /// Remove environment variables with names matching these patterns.
    std::vector<std::string> dropVars = DefaultDropVars;
    
    /// Look for the GIT version of these repositories.
    std::vector<std::string> repositories{ "sbncode" };
    
    std::string logCategory = "JobEnvironmentInfoExtractor";
    
    // needed because of compiler bugs
    // https://bugs.llvm.org/show_bug.cgi?id=36684
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96645
    Config() {}
    
  }; // Config
  
  
  /// Constructor with configuration.
  JobEnvironmentInfoExtractor(Config const& config = Config{});
  
  /// Extracts all the information from the environment.
  sbn::JobEnvironmentInfo extract
    (art::ModuleDescription const& moduleInfo); // not guaranteed to be const!
  
  /**
   * @brief Extracts all the information from the environment.
   * @param moduleInfo information about the _art_ module being executed
   * @return an object with all extracted metadata
   */
  sbn::JobEnvironmentInfo operator() (art::ModuleDescription const& moduleInfo)
    { return extract(moduleInfo); }
  
  // --- BEGIN ---  Partial metadata extraction  -------------------------------
  /**
   * @name Partial metadata extraction
   * 
   * The member functions in this group allow the extraction of part of the job
   * metadata. The function `extract()` includes the information from all of
   * them, plus some more.
   */
  /// @{
  
  /// Returns the environment variables.
  sbn::OrderedPairList extractEnvironmentVariables() const;
  
  /// Returns the version of the source repositories.
  sbn::OrderedPairList extractSourceVersions() const;
  
  /// @}
  // --- END -----  Partial metadata extraction  -------------------------------
  
    private:
  
  std::vector<std::regex> fDropVars; ///< Remove variables with matching names.
  
  ///< Names of the repositories to ask for a version report to.
  std::vector<std::string> fRepositoryReports;
  
  std::string fLogCategory;
  
  /// Converts all the patterns into regular expression objects.
  static std::vector<std::regex> prepareRegEx
    (std::vector<std::string> const& patterns);

  /// @brief Returns the index of the first pattern `value` matches.
  /// @return index of the first matched pattern, `patterns.size()` if none
  static std::size_t matchPatterns
    (std::string const& value, std::vector<std::regex> const& patterns);
  
}; // sbn::JobEnvironmentInfoExtractor


// -----------------------------------------------------------------------------

#endif // SBNCODE_METADATA_JOBENVIRONMENTINFOEXTRACTOR_H
