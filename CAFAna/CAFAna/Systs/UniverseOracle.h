#pragma once

#include <map>
#include <string>
#include <vector>

namespace ana
{
  class UniverseOracle
  {
  public:
    UniverseOracle();

    bool SystExists(const std::string& name) const;
    std::vector<std::string> Systs() const;

    const std::vector<double>& ShiftsForSyst(const std::string& name) const;
    unsigned int ClosestIndex(const std::string& name,
                              double shift,
                              double* trueShift = 0) const;
  protected:
    std::map<std::string, std::vector<double>> fData;
  };
}
