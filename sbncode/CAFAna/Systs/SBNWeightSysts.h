#pragma once

#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Core/Var.h"

#include <vector>
#include <unordered_map>

namespace caf
{
  class PairProxy;
  template<class T> class VectorProxy;
}

namespace ana
{
  class UniverseWeight
  {
  public:
    UniverseWeight(const std::vector<std::string>& systs, int univIdx);
    UniverseWeight(const std::vector<const ISyst*>& systs, int univIdx);

    double operator()(const caf::SRProxy* sr) const;
  protected:
    int GetIdx(const caf::VectorProxy<caf::PairProxy>& ws, int isyst) const;

    std::vector<std::string> fNames;
    int fUnivIdx;
    mutable std::vector<int> fSystIdxs;
  };

  Var GetUniverseWeight(const std::string& syst, int univIdx)
  {
    return Var(UniverseWeight(std::vector<std::string>(1, syst), univIdx));
  }

  Var GetUniverseWeight(const std::vector<std::string>& systs, int univIdx)
  {
    return Var(UniverseWeight(systs, univIdx));
  }

  Var GetUniverseWeight(const std::vector<const ISyst*> systs, int univIdx)
  {
    return Var(UniverseWeight(systs, univIdx));
  }


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

  const std::vector<const ISyst*>& GetSBNGenieWeightSysts();
  const std::vector<const ISyst*>& GetSBNFluxWeightSysts();
  const std::vector<const ISyst*>& GetSBNWeightSysts(); // genie+flux
}
