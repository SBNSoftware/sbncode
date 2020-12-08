#pragma once

// This file defines the basic Var object. For specific variables, and examples
// of how to implement your own, see Vars.h

#include "CAFAnaCore/CAFAna/Core/Var.h"

#include "StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  /// \brief Representation of a variable to be retrieved from a \ref
  /// caf::StandardRecord object
  ///
  /// A Var consists of a function, taking a StandardRecord and returning the
  /// value of the variable (which may be some complicated function).
  typedef _Var<caf::SRSliceProxy> Var;

  /// \brief Equivalent of \ref Var acting on \ref caf::SRSpill
  typedef _Var<caf::SRSpillProxy> SpillVar;

  /// \brief For Vars where literally all you need is a single CAF variable
  ///
  /// eg Var myVar = SIMPLEVAR(my.var.str);
  /// NB lack of quotes quotes around my.var.str
#define SIMPLEVAR(CAFNAME) Var([](const caf::SRSliceProxy* sr){return sr->CAFNAME;})

#define SIMPLESPILLVAR(CAFNAME) SpillVar([](const caf::SRSpillProxy* sr){return sr->CAFNAME;})

  /// The simplest possible Var, always 1. Used as a default weight.
  const Var kUnweighted = Unweighted<caf::SRSliceProxy>();

  const SpillVar kSpillUnweighted = Unweighted<caf::SRSpillProxy>();
} // namespace
