#pragma once

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/IFitVar.h"

namespace ana
{
  /// \brief Collect information describing the x-axis of an fit variable
  ///
  /// This will provide the support to creat e.g. log-scale fitting variables

  class FitAxis
  {
  public:

  FitAxis(const IFitVar* var,
            int nbins, double xmin, double xmax,
            bool isLogScale = false);
  
  ~FitAxis() {  };

  const IFitVar* fVar;
  int fnbins;
  double fxmin;
  double fxmax;
  bool fIsLogScale;

  };
}
