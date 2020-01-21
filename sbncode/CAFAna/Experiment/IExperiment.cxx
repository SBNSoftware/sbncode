#include "CAFAna/Experiment/IExperiment.h"

#include "CAFAna/Core/Utilities.h"

#include "TFile.h"
#include "TObjString.h"

#include <cassert>
#include <iostream>

// To implement LoadFrom()
#include "CAFAna/Experiment/CountingExperiment.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/SolarConstraints.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/ReactorExperiment.h"

namespace ana
{
  //----------------------------------------------------------------------
  // Definition to satisfy declaration in Core/LoadFromFile.h
  template<> std::unique_ptr<IExperiment> LoadFrom<IExperiment>(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);

    const TString tag = ptag->GetString();

    if(tag == "CountingExperiment") return CountingExperiment::LoadFrom(dir);
    if(tag == "ReactorExperiment") return ReactorExperiment::LoadFrom(dir);
    if(tag == "SingleSampleExperiment") return SingleSampleExperiment::LoadFrom(dir);
    if(tag == "SolarConstraints") return SolarConstraints::LoadFrom(dir);
    if(tag == "MultiExperiment") return MultiExperiment::LoadFrom(dir);

    std::cerr << "Unknown Experiment type '" << tag << "'" << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  void IExperiment::SaveTo(TDirectory* dir) const
  {
    assert(0 && "Not implemented");
  }
}
