/**
 * @file   sbncode/Metadata/RepositoryVersionReportUtils.h
 * @brief  Interface for _art_ tool reporting the version of packages.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 18, 2025
 * @see    sbncode/Metadata/RepositoryVersionReportUtils.cxx
 */

#ifndef SBNCODE_METADATA_REPOSITORYVERSIONREPORTUTILS_H
#define SBNCODE_METADATA_REPOSITORYVERSIONREPORTUTILS_H


// SBN libraries
#include "sbnobj/Common/Metadata/OrderedPairList.h"

// C/C++ standard library
#include <string>
#include <vector>


// -----------------------------------------------------------------------------
namespace sbn {
  
  /**
   * @brief Adds versions from the specified repository.
   * @param[out] versionList version list to add to
   * @param repoName name of the repository to be queried
   * @return whether the repository plugin could be loaded
   * 
   * The plugin with name `repoName + "RepositoryVersion"` is loaded and its
   * report is appended to `versionList`.
   * It will be necessary to call `versionList.finish()` after all the additions
   * are done.
   * 
   * This utility is used in the
   * @ref SBNsourceMetadataSystem "repository version tracking system".
   */
  bool addVersionFromRepository
    (sbn::OrderedPairList& versionList, std::string const& repoName);
  
  
  /**
   * @brief Adds versions from the specified repositories.
   * @param[out] versionList version list to add to
   * @param repoNames name of the repositories to be queried
   * @return which repositories which could not be reached
   * @see addVersionFromRepository
   * 
   * This is equivalent to repeatedly calling `addVersionFromRepository()`.
   */
  std::vector<std::string> addVersionFromRepositories(
    sbn::OrderedPairList& versionList, std::vector<std::string> const& repoNames
    );
  
} // namespace sbn


// -----------------------------------------------------------------------------

#endif // SBNCODE_METADATA_REPOSITORYVERSIONREPORTUTILS_H
