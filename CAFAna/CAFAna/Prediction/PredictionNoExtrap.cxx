#include "CAFAna/Prediction/PredictionNoExtrap.h"

#include "CAFAna/Extrap/IExtrap.h"
#include "CAFAna/Core/LoadFromFile.h"

#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Extrap/TrivialExtrap.h"

#include "TDirectory.h"
#include "TObjString.h"


namespace ana
{
  //----------------------------------------------------------------------
  PredictionNoExtrap::PredictionNoExtrap(SpectrumLoaderBase& loaderNonswap,
                                         SpectrumLoaderBase& loaderNue,
                                         SpectrumLoaderBase& loaderNuTau,
                                         const std::string& label,
                                         const Binning& bins,
                                         const Var& var,
                                         const Cut& cut,
                                         const SystShifts& shift,
                                         const Var& wei)
    : PredictionExtrap(new TrivialExtrap(loaderNonswap, loaderNue, loaderNuTau,
                                         label, bins, var, cut, shift, wei))
  {
  }

  //----------------------------------------------------------------------
  PredictionNoExtrap::PredictionNoExtrap(SpectrumLoaderBase& loaderNonswap,
                                         SpectrumLoaderBase& loaderNue,
                                         SpectrumLoaderBase& loaderNuTau,
					 const HistAxis& axis,
                                         const Cut& cut,
                                         const SystShifts& shift,
                                         const Var& wei)
    : PredictionExtrap(new TrivialExtrap(loaderNonswap, loaderNue, loaderNuTau,
                                         axis, cut, shift, wei))
  {
  }

  //----------------------------------------------------------------------
  PredictionNoExtrap::PredictionNoExtrap(PredictionExtrap* pred) : PredictionExtrap(pred->GetExtrap())
  {
  }

  //----------------------------------------------------------------------
  PredictionNoExtrap::PredictionNoExtrap(Loaders& loaders,
                                         const std::string& label,
                                         const Binning& bins,
                                         const Var& var,
                                         const Cut& cut,
                                         const SystShifts& shift,
                                         const Var& wei)
    : PredictionNoExtrap(loaders, HistAxis(label, bins, var), cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  PredictionNoExtrap::PredictionNoExtrap(Loaders& loaders,
                                         const HistAxis& axis,
                                         const Cut& cut,
                                         const SystShifts& shift,
                                         const Var& wei)
    : PredictionExtrap(new TrivialExtrap(loaders, axis, cut, shift, wei))
  {
  }

  //----------------------------------------------------------------------
  void PredictionNoExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("PredictionNoExtrap").Write("type");

    fExtrap->SaveTo(dir->mkdir("extrap"));

    tmp->cd();
  }


  //----------------------------------------------------------------------
  std::unique_ptr<PredictionNoExtrap> PredictionNoExtrap::LoadFrom(TDirectory* dir)
  {
    assert(dir->GetDirectory("extrap"));
    IExtrap* extrap = ana::LoadFrom<IExtrap>(dir->GetDirectory("extrap")).release();
    PredictionExtrap* pred = new PredictionExtrap(extrap);
    return std::unique_ptr<PredictionNoExtrap>(new PredictionNoExtrap(pred));
  }


  //----------------------------------------------------------------------
  PredictionNoExtrap::~PredictionNoExtrap()
  {
    // We created this in the constructor so it's our responsibility
    delete fExtrap;
  }
}
