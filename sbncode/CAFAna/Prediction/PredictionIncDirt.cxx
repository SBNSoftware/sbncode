#include "CAFAna/Prediction/PredictionIncDirt.h"

#include "CAFAna/Core/LoadFromFile.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  // --------------------------------------------------------------------------
  PredictionIncDirt::PredictionIncDirt(SpectrumLoaderBase& loaderNonswap,
                                       SpectrumLoaderBase& loaderNue,
                                       SpectrumLoaderBase& loaderNuTau,
                                       SpectrumLoaderBase& loaderIntrinsic,
                                       SpectrumLoaderBase& loaderDirt,
                                       const HistAxis& axis,
                                       const Cut& cut,
                                       const SystShifts& shift,
                                       const Var& wei)
    : fDet(loaderNonswap, loaderNue, loaderNuTau, loaderIntrinsic,
           axis, cut, shift, wei),
      fDirt(loaderDirt, kNullLoader, kNullLoader, kNullLoader,
            axis, cut, shift, wei)
  {
  }

  // --------------------------------------------------------------------------
  PredictionIncDirt::PredictionIncDirt(Loaders& loaders,
                                       SpectrumLoaderBase& loaderDirt,
                                       const HistAxis& axis,
                                       const Cut& cut,
                                       const SystShifts& shift,
                                       const Var& wei)
    : fDet(loaders, axis, cut, shift, wei),
      fDirt(loaderDirt, kNullLoader, kNullLoader, kNullLoader,
            axis, cut, shift, wei)
  {
  }

  // --------------------------------------------------------------------------
  PredictionIncDirt::~PredictionIncDirt()
  {
  }

  // --------------------------------------------------------------------------
  std::unique_ptr<PredictionIncDirt>
  PredictionIncDirt::LoadFrom(TDirectory* dir)
  {
    assert(dir->GetDirectory("det") && dir->GetDirectory("dirt"));

    return std::unique_ptr<PredictionIncDirt>(new PredictionIncDirt(ana::LoadFrom<PredictionNoExtrap>(dir->GetDirectory("det")),
                                                                    ana::LoadFrom<PredictionNoExtrap>(dir->GetDirectory("dirt"))));
  }

  // --------------------------------------------------------------------------
  void PredictionIncDirt::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("PredictionIncDirt").Write("type");

    fDet.SaveTo(dir->mkdir("det"));
    fDirt.SaveTo(dir->mkdir("dirt"));

    tmp->cd();
  }
}
