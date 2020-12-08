#pragma once

#include "CAFAnaCore/CAFAna/Core/HistAxis.h"

#include "CAFAna/Core/Var.h" // TODO do we want our own FwdDeclare.h?

namespace ana
{
  typedef _HistAxis<Var> HistAxis;
  typedef _HistAxis<SpillVar> SpillHistAxis;
}
