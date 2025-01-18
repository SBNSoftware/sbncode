/**
 * @file   sbncode/Metadata/GlobalMetadataRegistry.h
 * @brief  Centralized metadata tracker.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 18, 2025
 */


#ifndef SBNCODE_METADATA_GLOBALMETADATAREGISTRY_H
#define SBNCODE_METADATA_GLOBALMETADATAREGISTRY_H


// SBN libraries
#include "sbnobj/Common/Metadata/OrderedPairList.h"

// C++ standard libraries
#include <mutex>


// -----------------------------------------------------------------------------
namespace sbn { class GlobalMetadataRegistry; }
/**
 * @brief A simple class with singleton access storing data in concurrent way.
 * 
 * The registry is a simple key-value (string-string) sorted list (implemented
 * from `OrderedPairList`).
 * 
 * It provides a singleton instance (`instance()`) and it is designed to be
 * thread-safe (by blunt mutex'ing).
 * 
 */
class sbn::GlobalMetadataRegistry {
  
  static GlobalMetadataRegistry GlobalRegistry; ///< Singleton instance.
  
  mutable std::mutex fAccessMutex;
  
  sbn::OrderedPairList fRegistry;
  
  
  /// Grant read access.
  [[nodiscard]] std::lock_guard<std::mutex> readAccess() const;
  
  /// Grants (exclusive) write access.
  [[nodiscard]] std::lock_guard<std::mutex> writeAccess() const;
  
  
    public:
  
  // These types need to match the ones in sbn::OrderedPairList
  using Name_t = std::string; ///< Type of metadata name.
  using Value_t = std::string; ///< Type of metadata value.
  
  
  /// Adds a key-value pair to the registry
  /// @return this object
  GlobalMetadataRegistry& registerMetadata(Name_t name, Value_t value);
  
  
  // --- BEGIN ---  Queries  ---------------------------------------------------
  /// @name Queries
  /// @{
  
  /// Returns the complete registry of metadata. **This is not thread-safe!**
  sbn::OrderedPairList const& registry() const
    { return fRegistry; }
  
  /// Returns a pointer to the value of the requested metadata `name`,
  /// or `nullptr` if none.
  Value_t const* getPtr(Name_t const& name) const;
  
  /// Returns the value of requested metadata `name`, `std::nullopt` if none.
  std::optional<Value_t> get(Name_t const& name) const;
  
  /// Returns the value of the requested metadata `name`, or `defValue` if none.
  Value_t get(Name_t const& name, Value_t const& defValue) const;
  
  /// @}
  // --- END -----  Queries  ---------------------------------------------------
  
  
  /// Returns the global instance of this class (read/write).
  static GlobalMetadataRegistry& instance() { return GlobalRegistry; }
  
}; // GlobalMetadataRegistry


#endif // SBNCODE_METADATA_GLOBALMETADATAREGISTRY_H
