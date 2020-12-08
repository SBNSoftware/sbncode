#pragma once

// This file defines the basic Cut object. For specific cuts, and examples of
// how to implement your own, see Cuts.h

#include "CAFAnaCore/CAFAna/Core/Cut.h"

#include "StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  /// \brief Representation of a cut (selection) to be applied to a \ref
  /// caf::StandardRecord object
  ///
  /// A Cut consists of a function, taking a StandardRecord and returning a
  /// boolean indicating if that event passes the cut.
  ///
  /// Cut objects may be combined with the standard boolean operations && ||
  /// and !
  typedef _Cut<caf::SRSliceProxy> SliceCut;
  typedef _Cut<caf::SRSliceProxy> Cut;

  /// \brief Equivalent of \ref Cut acting on \ref caf::SRSpill. For use in
  /// spill-by-spill data quality cuts
  typedef _Cut<caf::SRSpillProxy> SpillCut;

  /// The simplest possible cut: pass everything, used as a default
  const Cut kNoCut([](const caf::SRSliceProxy*){return true;});

  /// The simplest possible cut: pass everything, used as a default
  const SpillCut kNoSpillCut([](const caf::SRSpillProxy*){return true;});
} // namespace
