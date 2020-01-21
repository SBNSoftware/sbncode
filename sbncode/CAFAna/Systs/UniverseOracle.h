#pragma once

#include "CAFAna/Core/SystShifts.h"

#include <map>
#include <string>
#include <vector>

namespace ana
{
  enum class ESide{
    kAbove, kBelow, kEither
  };

  class UniverseOracle
  {
  public:
    static UniverseOracle& Instance();

    bool SystExists(const std::string& name) const;
    std::vector<std::string> Systs() const;

    const std::vector<double>& ShiftsForSyst(const std::string& name) const;

    std::vector<SystShifts> ShiftsForSysts(const std::vector<const ISyst*>& systs, int nUniv) const;

    unsigned int ClosestIndex(const std::string& name,
                              double shift,
                              ESide side = ESide::kEither,
                              double* trueShift = 0) const;
  protected:
    UniverseOracle();

    std::map<std::string, std::vector<double>> fData;
  };
}
