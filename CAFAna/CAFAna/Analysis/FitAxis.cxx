#include "CAFAna/Analysis/FitAxis.h"

namespace ana
{
  // --------------------------------------------------------------------------
  FitAxis::FitAxis(const IFitVar* var_,
                   int nbins_, double min_, double max_,
                   bool islog_)
    : var(var_),
      nbins(nbins_),
      min(min_),
      max(max_),
      islog(islog_)
  {
  }
}
