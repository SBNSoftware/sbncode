#pragma once

#include "CAFAna/Core/ISyst.h"

#include <vector>
#include <unordered_map>

namespace caf
{
  class PairProxy;
  template<class T> class VectorProxy;
}

namespace ana
{
  class SBNWeightSyst: public ISyst
  {
  public:
    SBNWeightSyst(const std::string& name);

    void Shift(double x, caf::SRProxy* sr, double& weight) const override;

  protected:
    mutable int fIdx;

    int GetIdx(const caf::VectorProxy<caf::PairProxy>& ws) const;


    struct Univs
    {
      int i0, i1;
      double w0, w1;
    };

    mutable std::unordered_map<double, Univs> fUnivs;

    Univs GetUnivs(double x) const;
  };

  const std::vector<const ISyst*>& GetSBNWeightSysts();
}
