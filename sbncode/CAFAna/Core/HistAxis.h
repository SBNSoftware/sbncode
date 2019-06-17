#pragma once

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

#include <string>

namespace ana
{
  /// \brief Collect information describing the x-axis of an analysis histogram
  ///
  /// That is, what it should be labelled, what the binning should be, and what
  /// variable will be being filled into it.
  template<class T> class GenericHistAxis
  {
  public:
    GenericHistAxis(const std::string& label,
                    const Binning& bins,
                    const T& var);

    GenericHistAxis(const std::vector<std::string>& labels,
                    const std::vector<Binning>& bins,
                    const std::vector<T>& vars);

    GenericHistAxis(const std::string& labelX,
                    const Binning& binsX,
                    const T& varX,
                    const std::string& labelY,
                    const Binning& binsY,
                    const T& varY);

    /// Shortcut for simple binnings
    GenericHistAxis(const std::string& label,
                    int nx, double x0, double x1,
                    const T& var);

    /// Shortcut for simple binnings
    GenericHistAxis(const std::string& labelX,
                    int nx, double x0, double x1,
                    const T& varX,
                    const std::string& labelY,
                    int ny, double y0, double y1,
                    const T& varY);


    unsigned int NDimensions() const{return fLabels.size();}

    std::vector<std::string> GetLabels() const {return fLabels;}
    std::vector<Binning> GetBinnings() const {return fBins;}
    std::vector<T> GetVars() const {return fVars;}

    /// A variable "flattening" all the dimensions into one 1D value. Use
    /// sparingly.
    T GetMultiDVar() const;

  protected:
    std::vector<std::string> fLabels;
    std::vector<Binning> fBins;
    std::vector<T> fVars;
  };

  typedef GenericHistAxis<Var> HistAxis;
  typedef GenericHistAxis<SpillVar> SpillHistAxis;
  typedef GenericHistAxis<SpillTruthVar> SpillTruthHistAxis;
}
