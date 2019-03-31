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
    const double scale = 1 + .10*sigma; // 10% E scale syst.                                       
    sr->truth[0].neutrino.energy *= scale;
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
    if(sr->truth[0].neutrino.energy > 2) weight *= TMath::Max(0., 1+0.2*sigma);
    else weight *= TMath::Max(0., 1+0.1*sigma);
  }
};
const ToyNormSyst& GetNSyst()
{
  static const ToyNormSyst nSyst;
  return nSyst;
}
