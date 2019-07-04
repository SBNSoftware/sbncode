#include "CAFAna/Core/ISyst.h"

using namespace ana;

class ToyEnergyScaleSyst: public ISyst
{
 public:
 ToyEnergyScaleSyst() : ISyst("toyEScale", "Toy Energy Scale") {}
  void Shift(double sigma,
             caf::SRProxy* sr,
             double& weight) const override
  {
    const double scale = 1 + .05*sigma; // 5% E scale syst.                    

    sr->reco.reco_energy *= scale;
  }
};
const ToyEnergyScaleSyst& GetESyst()
{
  static const ToyEnergyScaleSyst eSyst;
  return eSyst;
}

class ToyNormSyst: public ISyst
{
 public:
 ToyNormSyst() : ISyst("toyNorm", "Toy Norm Scale") {}
  void Shift(double sigma,
             caf::SRProxy* sr,
             double& weight) const override
  {
    weight *= TMath::Max(0., 1+0.1*sigma);
  }
};
const ToyNormSyst& GetNSyst()
{
  static const ToyNormSyst nSyst;
  return nSyst;
}

// Make a vector with all (here only two) the systematics
std::vector<const ISyst*> allSysts{&GetESyst(), &GetNSyst()};
