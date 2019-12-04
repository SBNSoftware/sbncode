#if 0
// Currently dead code

#pragma once

#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Extrap/IExtrap.h"
#include "CAFAna/Extrap/ModularExtrapComponent.h"
#include <memory>

class TDirectory;

namespace ana
{

  class Loaders;
  class SpectrumLoaderBase;
  class IDecomp;
  class OscillatableSpectrum;

  /// \brief Extrapolate each component using a separate ModularExtrapComponent
  ///
  /// Only extrapolates one sample.
  class ModularExtrap: public IExtrap
  {
    public:

      // /// Creates a nue-like extraploation.
      // /** numuDecomp and numuNDCut are for the signal component (mu->e).
      //     nueDecomp and nueNDCut are for the backgound components. */
      // static ModularExtrap Nue(
      //   Loaders& loaders,
      //   const IDecomp& nueDecomp,
      //   const IDecomp& numuDecomp,
      //   const HistAxis& axis,
      //   const HistAxis& axisNumuND,
      //   const Cut& fdcut,
      //   const Cut& nueNDcut,
      //   const Cut& numuNDcut,
      //   const SystShifts& shiftMC = kNoShift,
      //   const Var& weight = kUnweighted
      // );

      /// Creates a numu-like extrapolation.
      static ModularExtrap Numu(
        Loaders& loaders,
        const IDecomp& numuDecomp,
        const HistAxis& axis,
        const Cut& fdcut,
        const Cut& ndcut,
        const SystShifts& shiftMC = kNoShift,
        const Var& weight = kUnweighted
      );

      // /// Creates a nue-like extraploation with individual spectrum loaders.
      // /** numuDecomp and numuNDCut are for the signal component (mu->e).
      //     nueDecomp and nueNDCut are for the backgound components. */
      // static ModularExtrap Nue(
      //   SpectrumLoaderBase& nearMCLoader,
      //   SpectrumLoaderBase& farMCswapLoader,
      //   SpectrumLoaderBase& farMCnonswapLoader,
      //   SpectrumLoaderBase& farMCtauswapLoader,
      //   const IDecomp& nueDecomp,
      //   const IDecomp& numuDecomp,
      //   const HistAxis& axis,
      //   const HistAxis& axisNumuND,
      //   const Cut& fdcut,
      //   const Cut& nueNDcut,
      //   const Cut& numuNDcut,
      //   const SystShifts& shiftMC = kNoShift,
      //   const Var& weight = kUnweighted
      // );

      /// Creates a numu-like extrapolation with individual spectrum loaders.
      static ModularExtrap Numu(
        SpectrumLoaderBase& nearMCLoader,
        SpectrumLoaderBase& farMCnonswapLoader,
        SpectrumLoaderBase& farMCnueLoader,
        SpectrumLoaderBase& farMCnutauLoader,
        SpectrumLoaderBase& farMCncLoader,
        const IDecomp& numuDecomp,
        const HistAxis& axis,
        const Cut& fdcut,
        const Cut& ndcut,
        const SystShifts& shiftMC = kNoShift,
        const Var& weight = kUnweighted
      );

      // Prevent copying because we own objects on the free store.
      ModularExtrap(const ModularExtrap&) = delete;
      ModularExtrap& operator=(const ModularExtrap&) = delete;
      ModularExtrap(ModularExtrap&&) = default;
      ModularExtrap& operator=(ModularExtrap&&) = default;
      virtual ~ModularExtrap() = default;

      void SaveTo(TDirectory* dir) const override;
      void SavePlotsNue( TDirectory* dir, double potFD ) const;
      void SavePlotsNumu( TDirectory* dir, double potFD ) const;
      static std::unique_ptr<ModularExtrap> LoadFrom(TDirectory* dir);

      // Override abstract methods.
      OscillatableSpectrum NueSurvComponent()       override;
      OscillatableSpectrum AntiNueSurvComponent()   override;
      OscillatableSpectrum NumuSurvComponent()      override;
      OscillatableSpectrum AntiNumuSurvComponent()  override;
      OscillatableSpectrum NueAppComponent()        override;
      OscillatableSpectrum AntiNueAppComponent()    override;
      OscillatableSpectrum NumuAppComponent()       override;
      OscillatableSpectrum AntiNumuAppComponent()   override;
      OscillatableSpectrum TauFromMuComponent()     override;
      OscillatableSpectrum AntiTauFromMuComponent() override;
      OscillatableSpectrum TauFromEComponent()      override;
      OscillatableSpectrum AntiTauFromEComponent()  override;
      Spectrum             NCComponent()            override;

      std::vector<ModularExtrapComponent*> GetModExtrapComponents() const
      {
      return {
	fEEextrap.get(), fEEAntiextrap.get(),
	fMMextrap.get(), fMMAntiextrap.get(),
	fMEextrap.get(), fMEAntiextrap.get(),
	fEMextrap.get(), fEMAntiextrap.get(),
	fNCextrap.get(),
	fMTextrap.get(), fMTAntiextrap.get(),
	fETextrap.get(), fETAntiextrap.get()};
      }

    protected:

      /// Sets up all components to use FD MC--internal use only.
      /** Use a named constructor (the static function Nue() or Numu()) to
          create a ModularExtrap.  This function is protected. */
      ModularExtrap(
        SpectrumLoaderBase& farMCnonswapLoader,
        SpectrumLoaderBase& farMCnueLoader,
        SpectrumLoaderBase& farMCnutauLoader,
        SpectrumLoaderBase& farMCncLoader,
        const HistAxis& axis,
        const Cut& fdcut,
        const SystShifts& shiftMC,
        const Var& weight
      );

      std::unique_ptr<ModularExtrapComponent> fEEextrap;
      std::unique_ptr<ModularExtrapComponent> fEEAntiextrap;
      std::unique_ptr<ModularExtrapComponent> fMMextrap;
      std::unique_ptr<ModularExtrapComponent> fMMAntiextrap;
      std::unique_ptr<ModularExtrapComponent> fMEextrap;
      std::unique_ptr<ModularExtrapComponent> fMEAntiextrap;
      std::unique_ptr<ModularExtrapComponent> fEMextrap;
      std::unique_ptr<ModularExtrapComponent> fEMAntiextrap;
      std::unique_ptr<ModularExtrapComponent> fNCextrap;
      std::unique_ptr<ModularExtrapComponent> fMTextrap;
      std::unique_ptr<ModularExtrapComponent> fMTAntiextrap;
      std::unique_ptr<ModularExtrapComponent> fETextrap;
      std::unique_ptr<ModularExtrapComponent> fETAntiextrap;

    private:

      ModularExtrap(){};

  };

}

#endif
