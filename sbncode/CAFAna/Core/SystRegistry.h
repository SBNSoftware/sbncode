#pragma once

#include "CAFAna/Core/ISyst.h"

#include <map>
#include <string>

namespace ana
{
  class SystRegistry
  {
  public:
    static void Register(const ISyst* s);

    static void UnRegister(const ISyst* s);

    static const ISyst* ShortNameToSyst(const std::string& s,
                                        bool allowFail = false);

    static void Print();
  protected:
    static std::map<std::string, const ISyst*>& Map();
  };
}
