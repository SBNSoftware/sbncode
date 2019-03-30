#pragma once

namespace ana
{
  class IFitVar;

  /// \brief Collect information describing the axis of a fit variable
  ///
  /// This will provide the support to create e.g. log-scale fitting variables
  class FitAxis
  {
  public:
    FitAxis(const IFitVar* var,
            int nbins, double min, double max,
            bool islog = false);
    
    const IFitVar* var;
    int nbins;
    double min;
    double max;
    bool islog;
  };
}
