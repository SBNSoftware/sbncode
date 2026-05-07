/**
 * @file   sbncode/CAFMaker/SRDefaults.h
 * @brief  Discovery of default values in Standard Record classes.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 6, 2026
 * @see    sbncode/CAFMaker/SRDefaults.cxx
 */


#ifndef SBNCODE_CAFMAKER_SRDEFAULTS_H
#define SBNCODE_CAFMAKER_SRDEFAULTS_H

// Standard Record objects
#include "sbnanaobj/StandardRecord/SRFlashMatch.h"
#include "sbnanaobj/StandardRecord/SRTPCPMTBarycenterMatch.h"

// C/C++ standard library
#include <tuple>
#include <type_traits> // std::decay_t

// -----------------------------------------------------------------------------
namespace caf { class SRDefaults; }
/**
 * @brief Collection of Standard Record objects with defaulted values.
 * 
 * Example of usage: to test if `caf::SRFlashMatch::time` is still the default
 * value:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * if (match.time == caf::SRDefaults::For<caf::SRFlashMatch>().time) { ... }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Note that this may be, and usually is, different from the
 * default-constructed value (e.g. `caf::SRFlashMatch{}.time`).
 * 
 * 
 * Extension directions
 * ---------------------
 * 
 * Only classes that offer `setDefault()` member function can be included in
 * this class. To add such a class:
 *  1. include the appropriate header above;
 *  2. add the name of the class as a `tuple` type in the definition of
 *     `SupportedTypes`.
 * 
 * In case a class is used that is not supported, a static error on `std::get()`
 * is typically emitted. In GCC 12.1.0:
 * ```
 * static assertion failed: the type T in std::get<T> must occur exactly once in the tuple
 * ```
 * 
 */
class caf::SRDefaults {
  
    public:
  
  /// A tuple with all the types included in this object.
  using SupportedTypes = std::tuple<
    caf::SRFlashMatch,
    caf::SRTPCPMTBarycenterMatch
    >;
  
  /// Single copy of defaulted objects.
  static SupportedTypes const defaultObjs;
  
  
  /**
   * @brief Returns an object of type `SRobj` with its default values set.
   * @tparam SRobj type of the Standard Record object to obtain
   * @return a static object of type `SRobj` with its default values set
   * 
   * A call of this method must explicitly include the template type `SRobj`,
   * e.g. `caf::SRDefaults::For<caf::SRTrack>()`.
   */
  template <typename SRobj>
  static SRobj const& For();
  
  /**
   * @brief Returns an object of type `SRobj` with its default values set.
   * @tparam SRobj type of the Standard Record object to obtain
   * @return a static object of type `SRobj` with its default values set
   * 
   * A call to this method should have as argument an existing object of the
   * `SRobj` type, that can be used to omit the explicit type in the call;
   * e.g., `caf::SRDefaults::For(track)`.
   */
  template <typename SRobj>
  static SRobj const& For(SRobj const&);
  
};


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
template <typename SRobj>
SRobj const& caf::SRDefaults::For() { return std::get<SRobj>(defaultObjs); }

template <typename SRobj>
SRobj const& caf::SRDefaults::For(SRobj const&)
  { return For<std::decay_t<SRobj>>(); }

// -----------------------------------------------------------------------------

#endif // SBNCODE_CAFMAKER_SRDEFAULTS_H
