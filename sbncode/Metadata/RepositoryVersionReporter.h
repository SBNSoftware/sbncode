/**
 * @file   sbncode/Metadata/RepositoryVersionReporter.h
 * @brief  Interface for _art_ tool reporting the version of packages.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 18, 2025
 * 
 */

#ifndef SBNCODE_METADATA_REPOSITORYVERSIONREPORTER_H
#define SBNCODE_METADATA_REPOSITORYVERSIONREPORTER_H


// SBN libraries
#include "sbnobj/Common/Metadata/OrderedPairList.h"

// C/C++ standard library
#include <iterator> // std::distance()


// -----------------------------------------------------------------------------
namespace sbn { class RepositoryVersionReporter; }
/**
 * @brief Interface for _art_ tool reporting the version of packages.
 * 
 * An object reporting package versions will return them as a range of pairs
 * package name/package version, each element being a string.
 * 
 * This interface is not at all generic, and the iterators that are returned
 * are exactly `const_iterator` ones defined in this interface, which are
 * iterators of the collection `vector_t` also defined here, which is exactly
 * a `std::vector`.
 * 
 * @note At this point there is no polymorphic interface: derived classes will
 *       fetch their data at construction time and keep it for anybody curious
 *       to learn it. In the future an actual virtual interface could become
 *       necessary, and at that point it would be clearer what it should
 *       contain.
 *
 */
struct sbn::RepositoryVersionReporter {
  
  /// Type of the name of a package (a `std::string`).
  using PackageName_t    = sbn::OrderedPairList::Key_t;
  
  /// Type of the version of a package (a `std::string`).
  using PackegaVersion_t = sbn::OrderedPairList::Value_t;
  
  // standard C++ container types
  using value_type = sbn::OrderedPairList::value_type;
  using size_type = sbn::OrderedPairList::size_type;
  using const_iterator = sbn::OrderedPairList::const_iterator;
  
  
  /// Collected packages and their versions.
  sbn::OrderedPairList packageVersions;
  
  
  /// Begin-iterator for the package-version pairs.
  const_iterator begin() const { return packageVersions.begin(); }
  
  /// End-iterator for the package-version pairs.
  const_iterator end() const { return packageVersions.end(); }
  
  /// Returns whether there is no entry in the report.
  bool empty() const noexcept { return begin() == end(); }
  
  /// Returns the number of items in the report.
  size_type size() const { return std::distance(begin(), end()); }
  
}; // sbn::RepositoryVersionReporter


// -----------------------------------------------------------------------------

#endif // SBNCODE_METADATA_REPOSITORYVERSIONREPORTER_H
