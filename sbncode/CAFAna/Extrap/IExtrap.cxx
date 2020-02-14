#include "CAFAna/Extrap/IExtrap.h"

#include "CAFAna/Core/LoadFromFile.h"

#include "TDirectory.h"
#include "TObjString.h"

#include <iostream>

// To implement LoadFrom()
#include "CAFAna/Extrap/TrivialExtrap.h"

namespace ana
{
  //----------------------------------------------------------------------
  // Definition to satisfy declaration in Core/LoadFromFile.h
  template<> std::unique_ptr<IExtrap> LoadFrom<IExtrap>(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);

    const TString tag = ptag->GetString();

    if(tag == "TrivialExtrap") return TrivialExtrap::LoadFrom(dir);

    std::cerr << "Unknown Extrap type '" << tag << "'" << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  void IExtrap::SaveTo(TDirectory* dir) const
  {
    assert(0 && "Not implemented");
  }
}
