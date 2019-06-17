#pragma once

#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include "TFile.h"

class TDirectory;

namespace osc{class IOscCalculator;}

namespace ana
{
  //----------------------------------------------------------------------
  // Most classes are happy to load themselves
  template<class T> std::unique_ptr<T> LoadFrom(TDirectory* dir)
  {
    return T::LoadFrom(dir);
  }

  //----------------------------------------------------------------------
  // But if you're trying to load a base class we need to figure out which
  // derived class is actually in the file and hand off to that. The
  // implementations of these are in the cxx files for the base classes in
  // question.
  class IDecomp;
  template<> std::unique_ptr<IDecomp> LoadFrom<IDecomp>(TDirectory* dir);
  class IExtrap;
  template<> std::unique_ptr<IExtrap> LoadFrom<IExtrap>(TDirectory* dir);
  class IPrediction;
  template<> std::unique_ptr<IPrediction> LoadFrom<IPrediction>(TDirectory* dir);
  class IExperiment;
  template<> std::unique_ptr<IExperiment> LoadFrom<IExperiment>(TDirectory* dir);
  class ModularExtrapComponent;
  template<> std::unique_ptr<ModularExtrapComponent>
    LoadFrom<ModularExtrapComponent>(TDirectory* dir);
  class IBkgdEstimator;
  template<> std::unique_ptr<IBkgdEstimator> LoadFrom<IBkgdEstimator>(TDirectory* dir);

  // This one is actually implemented in LoadFromFile.cxx to avoid polluting
  // OscLib with CAFAna conventions.
  template<> std::unique_ptr<osc::IOscCalculator> LoadFrom<osc::IOscCalculator>(TDirectory* dir);

  //----------------------------------------------------------------------
  // For symmetry
  template<class T> void SaveTo(const T& x, TDirectory* dir)
  {
    x.SaveTo(dir);
  }

  // Also in the cxx, to avoid having to put this logic into OscLib
  template<> void SaveTo(const osc::IOscCalculator& x, TDirectory* dir);

  //----------------------------------------------------------------------
  template<class T> std::unique_ptr<T> LoadFromFile(const std::string& fname,
                                                    const std::string& label)
  {
    TFile fin(fname.c_str());
    assert(!fin.IsZombie());
    TDirectory* dir = fin.GetDirectory(label.c_str());
    if(!dir){
      std::cerr << "Didn't find '" << label << "' in " << fname << std::endl;
      abort();
    }
    return LoadFrom<T>(dir);
  }

  //----------------------------------------------------------------------
  template<class T> void SaveToFile(const T& x,
                                    const std::string& fname,
                                    const std::string& label)
  {
    TFile fout(fname.c_str(), "RECREATE");
    x.SaveTo(fout.mkdir(label.c_str()));
  }
}
