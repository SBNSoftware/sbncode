/**
 * @file   sbncode/Metadata/GlobalMetadataRegistry.cxx
 * @brief  Centralized metadata tracker (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 18, 2025
 * @file   sbncode/Metadata/GlobalMetadataRegistry.h
 */

// library header
#include "sbncode/Metadata/GlobalMetadataRegistry.h"

// -----------------------------------------------------------------------------
// global instance:
sbn::GlobalMetadataRegistry sbn::GlobalMetadataRegistry::GlobalRegistry;


// -----------------------------------------------------------------------------
std::lock_guard<std::mutex> sbn::GlobalMetadataRegistry::readAccess() const {
  return std::lock_guard{ fAccessMutex };
}


// -----------------------------------------------------------------------------
std::lock_guard<std::mutex> sbn::GlobalMetadataRegistry::writeAccess() const {
  return std::lock_guard{ fAccessMutex };
}


// -----------------------------------------------------------------------------
auto sbn::GlobalMetadataRegistry::registerMetadata
  (Name_t name, Value_t value) -> GlobalMetadataRegistry&
{
  
  auto const accessGrant = writeAccess();
  
  fRegistry.insert(std::move(name), std::move(value));
  
  return *this;
  
} // sbn::GlobalMetadataRegistry::registerMetadata()


// -----------------------------------------------------------------------------
auto sbn::GlobalMetadataRegistry::getPtr(Name_t const& name) const
  -> Value_t const*
{
  auto const accessGrant = readAccess();
  return registry().getPtr(name);
}

// -----------------------------------------------------------------------------
auto sbn::GlobalMetadataRegistry::get(Name_t const& name) const
  -> std::optional<Value_t>
{
  auto const accessGrant = readAccess();
  return registry().get(name);
}


// -----------------------------------------------------------------------------
auto sbn::GlobalMetadataRegistry::get
  (Name_t const& name, Value_t  const& defValue) const -> Value_t
{
  auto const accessGrant = readAccess();
  return registry().get(name, defValue); 
}


// -----------------------------------------------------------------------------
