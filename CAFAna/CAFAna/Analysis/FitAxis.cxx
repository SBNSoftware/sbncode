#include "CAFAna/Analysis/FitAxis.h"

#include <iostream>

namespace ana
{
  FitAxis::FitAxis(const IFitVar* var,
          int nbins, double xmin, double xmax,
          bool isLogScale)
    : fVar(var),
      fnbins(nbins),
      fxmin(xmin),
      fxmax(xmax),
      fIsLogScale(isLogScale)
  {
  }

}
