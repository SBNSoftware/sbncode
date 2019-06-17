#include "CAFAna/Core/HistAxis.h"

#include <iostream>

namespace ana
{
  // explicitly instantiate the template for the types we know we have
  template class GenericHistAxis<Var>;
  template class GenericHistAxis<SpillVar>;
  template class GenericHistAxis<SpillTruthVar>;

  //----------------------------------------------------------------------
  template<class T> GenericHistAxis<T>::
  GenericHistAxis(const std::string& label,
                  const Binning& bins,
                  const T& var)
    : fLabels(1, label),
      fBins(1, bins),
      fVars(1, var)
  {
  }

  //----------------------------------------------------------------------
  template<class T> GenericHistAxis<T>::
  GenericHistAxis(const std::vector<std::string>& labels,
                  const std::vector<Binning>& bins,
                  const std::vector<T>& vars)
    : fLabels(labels), fBins(bins), fVars(vars)
  {
    assert(fLabels.size() == fBins.size());
    assert(fBins.size() == fVars.size());
  }

  //----------------------------------------------------------------------
  template<class T> GenericHistAxis<T>::
  GenericHistAxis(const std::string& labelX,
                  const Binning& binsX,
                  const T& varX,
                  const std::string& labelY,
                  const Binning& binsY,
                  const T& varY)
    : fLabels({labelX, labelY}),
      fBins({binsX, binsY}),
      fVars({varX, varY})
  {
  }

  //----------------------------------------------------------------------
  template<class T> GenericHistAxis<T>::
  GenericHistAxis(const std::string& label,
                  int nx, double x0, double x1,
                  const T& var)
    : GenericHistAxis(label, Binning::Simple(nx, x0, x1), var)
  {
  }

  //----------------------------------------------------------------------
  template<class T> GenericHistAxis<T>::
  GenericHistAxis(const std::string& labelX,
                  int nx, double x0, double x1,
                  const T& varX,
                  const std::string& labelY,
                  int ny, double y0, double y1,
                  const T& varY)
    : GenericHistAxis(labelX, Binning::Simple(nx, x0, x1), varX,
                      labelY, Binning::Simple(ny, y0, y1), varY)
  {
  }

  //----------------------------------------------------------------------
  template<class T> T GenericHistAxis<T>::GetMultiDVar() const
  {
    switch(fVars.size()){
    case 1:
      return fVars[0];
    case 2:
      return Var2D(fVars[0], fBins[0],
                   fVars[1], fBins[1]);
    case 3:
      return Var3D(fVars[0], fBins[0],
                   fVars[1], fBins[1],
                   fVars[2], fBins[2]);
    default:
      std::cout << "Error: HistAxis::GetMultiDVar() doesn't support "
                << fVars.size() << "-dimensional axes" << std::endl;
      abort();
    }
  }
}
