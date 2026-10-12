/**
 * @file   sbncode/Utilities/AssnsUtils.h
 * @brief  Framework-level utilities dealing with associations.
 * @date   April 2, 2026
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * 
 * This library is currently header-only.
 */

#ifndef SBNCODE_UTILITIES_ASSNSUTILS_H
#define SBNCODE_UTILITIES_ASSNSUTILS_H

// framework libraries
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"

// C/C++ libraries
#include <cstddef> // std::size_t
#include <tuple> // std::get()
#include <type_traits> // std::is_same_v, ...


namespace sbn {
  /**
   * @brief Creates an association of a new product mirroring an existing one.
   * @tparam DiscardData (default: `false`) whether to omit copying the data
   * @tparam L type of left-side data in the association
   * @tparam R type of right-side data in the association
   * @tparam Data type of metadata in the association
   * @tparam LorR type whose pointer is being rebound
   * @param assns the existing association
   * @param reboundPtrMaker tool to create a _art_ pointer to the new product
   * @return the newly created association
   * 
   * The existing association `assns` is replicated replacing all the pointers
   * in the rebinding side (left or right) with their corresponding pointers
   * from another collection of the same type.
   * 
   * The application is for one-to-one processing of data product elements,
   * where the element resulting from the processing is still associated to the
   * same object as the processed element.
   * For example, let's have a reprocessing that modifies the signal on a
   * optical detector waveform (`raw::OpDetWaveform`), leaving the baseline
   * unchanged.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using Assns_t = art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>;
   * 
   * Assns_t const& oldAssns = event.getProduct<Assns_t>(assnsTag);
   * art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr{ event };
   * 
   * auto newAssns
   *   = std::make_unique<Assns_t>(sbn::RebindAssociatedProducts(oldAssns, makeWaveformPtr));
   * 
   * event.put(std::move(newAssns));
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The alias `Assns_t` was used only to make the example more readable.
   * 
   * As shown in the example, which side to rebind is determined by the type
   * of the `art::PtrMaker` object. If its operand type is neither of the
   * association types, a compilation error will occur.
   * 
   * Note that the invalid pointers are also rebound, but their rebound pointer
   * should be as invalid as the old one.
   * 
   * If metadata is present in the associations, it is _copied_ element by
   * element into the new one. It is possible to omit the metadata by calling
   * with `DiscardData` set to `true` (`RebindAssociatedProducts<true>(...)`),
   * in which case the returned type will always be an association without data.
   * 
   */
  template <
    bool DiscardData = false,
    typename L, typename R, typename Data,
    typename LorR
    >
  art::Assns<L, R, std::conditional_t<DiscardData, void, Data>>
  RebindAssociatedProducts(
    art::Assns<L, R, Data> const& assns,
    art::PtrMaker<LorR> const& reboundPtrMaker
    );

} // namespace sbn


// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
template <
  bool DiscardData /* = false */,
  typename L, typename R, typename Data,
  typename LorR
  >
art::Assns<L, R, std::conditional_t<DiscardData, void, Data>>
sbn::RebindAssociatedProducts(
  art::Assns<L, R, Data> const& assns,
  art::PtrMaker<LorR> const& reboundPtrMaker
) {
  static_assert(std::is_same_v<L, LorR> || std::is_same_v<R, LorR>,
    "The PtrMaker tool must operate on either of the association data types.");
  
  constexpr bool RebindLeft = std::is_same_v<L, LorR>;
  constexpr bool hasMetadata = !std::is_void_v<Data>;
  constexpr bool wantsMetadata = hasMetadata  && !DiscardData;
  
  using DestData_t = std::conditional_t<wantsMetadata, Data, void>;
  
  /*
   * How to deal with invalid pointers? (there should be none though!)
   *  a) preserve them (drawback: different product IDs are mixed in)
   *  b) replace them with default-constructed (which is also invalid)
   *  c) rebind it (drawback: rely on PtrMaker doing the right thing with them)
   * Here I chose (c).
   */
  auto rebind = [&reboundPtrMaker](auto& ptr)
    { /* if (ptr) */ ptr = reboundPtrMaker(ptr.key()); };
  
  art::Assns<L, R, DestData_t> newAssns;
  for (auto elements: assns) {
    
    art::Ptr<L> leftPtr  = std::move(std::get<0>(elements));
    art::Ptr<R> rightPtr = std::move(std::get<1>(elements));
    
    if constexpr(RebindLeft) { rebind(leftPtr);  }
    else                     { rebind(rightPtr); }
    
    if constexpr(wantsMetadata) {
      newAssns.addSingle
        (std::move(leftPtr), std::move(rightPtr), std::move(std::get<2>(elements)));
    }
    else {
      newAssns.addSingle(std::move(leftPtr), std::move(rightPtr));
    }
    
  } // for assn pairs
  
  return newAssns;
  
} // RebindAssociatedProducts()


// -----------------------------------------------------------------------------

#endif // SBNCODE_UTILITIES_ASSNSUTILS_H
