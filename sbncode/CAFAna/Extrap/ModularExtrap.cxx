#if 0
// Currently dead code

#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Extrap/ModularExtrap.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "TObjString.h"
#include "TDirectory.h"

namespace ana
{

  //---------------------------------------------------------------------------
  /*
  ModularExtrap ModularExtrap::Nue(
    Loaders& loaders,
    const IDecomp& nueDecomp,
    const IDecomp& numuDecomp,
    const HistAxis& axis,
    const HistAxis& axisNumuND,
    const Cut& fdcut,
    const Cut& nueNDcut,
    const Cut& numuNDcut,
    const SystShifts& shiftMC,
    const Var& weight
  ){
    return ModularExtrap::Nue(
      loaders.GetLoader(
        caf::kNEARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap),
      loaders.GetLoader(
        caf::kFARDET,  Loaders::kMC, ana::kBeam, Loaders::kFluxSwap),
      loaders.GetLoader(
        caf::kFARDET,  Loaders::kMC, ana::kBeam, Loaders::kNonSwap),
      loaders.GetLoader(
        caf::kFARDET,  Loaders::kMC, ana::kBeam, Loaders::kTauSwap),
      nueDecomp,
      numuDecomp,
      axis,
      axisNumuND,
      fdcut,
      nueNDcut,
      numuNDcut,
      shiftMC,
      weight
    );
  }
  */
  //---------------------------------------------------------------------------

  ModularExtrap ModularExtrap::Numu(
    Loaders& loaders,
    const IDecomp& numuDecomp,
    const HistAxis& axis,
    const Cut& fdcut,
    const Cut& ndcut,
    const SystShifts& shiftMC,
    const Var& weight
  ){
    return ModularExtrap::Numu(
      loaders.GetLoader(
        caf::kNEARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap),
      loaders.GetLoader(
        caf::kFARDET,  Loaders::kMC, ana::kBeam, Loaders::kNonSwap),
      loaders.GetLoader(
        caf::kFARDET,  Loaders::kMC, ana::kBeam, Loaders::kNue),
      loaders.GetLoader(
        caf::kFARDET,  Loaders::kMC, ana::kBeam, Loaders::kTau),
      loaders.GetLoader(
        caf::kFARDET,  Loaders::kMC, ana::kBeam, Loaders::kNC),
      numuDecomp,
      axis,
      fdcut,
      ndcut,
      shiftMC,
      weight
    );
  }

  //---------------------------------------------------------------------------
  /*
  ModularExtrap ModularExtrap::Nue(
    SpectrumLoaderBase& nearMC,
    SpectrumLoaderBase& farMCswap,
    SpectrumLoaderBase& farMCnonswap,
    SpectrumLoaderBase& farMCtauswap,
    const IDecomp& nueDecomp,
    const IDecomp& numuDecomp,
    const HistAxis& axis,
    const HistAxis& axisNumuND,
    const Cut& fdcut,
    const Cut& nueNDcut,
    const Cut& numuNDcut,
    const SystShifts& shiftMC,
    const Var& weight
  ){

    ModularExtrap extrap(
      farMCswap,
      farMCnonswap,
      farMCtauswap,
      axis,
      fdcut,
      shiftMC,
      weight
    );

    // mu -> mu  ----
    extrap.fMMextrap = std::unique_ptr<ModularExtrapComponent>(
      new RecoReweight(
        nearMC, axis, fdcut, shiftMC, weight,
        "mu -> mu", "#nu_{#mu} #rightarrow #nu_{#mu}",
        nueNDcut, nueDecomp,                         // nue selection in ND
        DecompResult::numu, kIsNumuCC && !kIsAntiNu, // numu truth in ND
        farMCnonswap, kIsNumuCC && !kIsAntiNu        // mu->mu in FD
      )
    );

    // mu -> e ----
    extrap.fMEextrap = std::unique_ptr<ModularExtrapComponent>(
      new TruthReweight(
        nearMC, axis, axisNumuND, fdcut, shiftMC, weight,
        "mu -> e", "#nu_{#mu} #rightarrow #nu_{e}",
        numuNDcut, numuDecomp,                       // numu selection in ND
        DecompResult::numu, kIsNumuCC && !kIsAntiNu, // numu truth in ND
        farMCswap, kIsSig && !kIsAntiNu              // mu->e in FD
      )
    );
    extrap.fMEAntiextrap = std::unique_ptr<ModularExtrapComponent>(
      new TruthReweight(
        nearMC, axis, axisNumuND, fdcut, shiftMC, weight,
        "mubar -> ebar", "#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}",
        numuNDcut, numuDecomp,                         // numu selection in ND
        DecompResult::numubar, kIsNumuCC && kIsAntiNu, // numubar truth in ND
        farMCswap, kIsSig && kIsAntiNu                 // mubar->ebar in FD
      )
    );

    // NC -> NC ----
    extrap.fNCextrap = std::unique_ptr<ModularExtrapComponent>(
      new RecoReweight(
        nearMC, axis, fdcut, shiftMC, weight,
        "NC -> NC", "NC #rightarrow NC",
        nueNDcut, nueDecomp,           // nue selection in ND
        DecompResult::NC, kIsNC,       // NC truth in ND
        farMCswap, kIsNC,              // NC->NC in FD
        farMCnonswap, farMCtauswap     // extra NC stats
      )
    );

    // e -> e ----
    extrap.fEEextrap = std::unique_ptr<ModularExtrapComponent>(
      new RecoReweight(
        nearMC, axis, fdcut, shiftMC, weight,
        "e -> e", "#nu_{e} #rightarrow #nu_{e}",
        nueNDcut, nueDecomp,                         // nue selection in ND
        DecompResult::nue, kIsBeamNue && !kIsAntiNu, // nue truth in ND
        farMCnonswap, kIsBeamNue && !kIsAntiNu       // e->e in FD
      )
    );

    return extrap;

  }
  */
  //---------------------------------------------------------------------------

  ModularExtrap ModularExtrap::Numu(
    SpectrumLoaderBase& nearMC,
    SpectrumLoaderBase& farMCnonswap,
    SpectrumLoaderBase& farMCnue,
    SpectrumLoaderBase& farMCnutau,
    SpectrumLoaderBase& farMCnc,
    const IDecomp& numuDecomp,
    const HistAxis& axis,
    const Cut& fdcut,
    const Cut& numuNDcut,
    const SystShifts& shiftMC,
    const Var& weight
  ){

    ModularExtrap extrap(
      farMCnonswap,
      farMCnue,
      farMCnutau,
      farMCnc,
      axis,
      fdcut,
      shiftMC,
      weight
    );

    // mu -> mu ----
    extrap.fMMextrap = std::unique_ptr<ModularExtrapComponent>(
      new TruthReweight(
        nearMC, axis, axis, fdcut, shiftMC, weight,
        "mu -> mu", "#nu_{#mu} #rightarrow #nu_{#mu}",
        numuNDcut, numuDecomp,                       // numu selection in ND
        DecompResult::numu, kIsNumuCC && !kIsAntiNu, // numu truth in ND
        farMCnonswap, kIsNumuCC && !kIsAntiNu        // mu->mu in FD
      )
    );
    extrap.fMMAntiextrap = std::unique_ptr<ModularExtrapComponent>(
      new TruthReweight(
        nearMC, axis, axis, fdcut, shiftMC, weight,
        "mubar -> mubar", "#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{#mu}",
        numuNDcut, numuDecomp,                         // numu selection in ND
        DecompResult::numubar, kIsNumuCC && kIsAntiNu, // numubar truth in ND
        farMCnonswap, kIsNumuCC && kIsAntiNu           // mubar->mubar in FD
      )
    );

    return extrap;
  }

  //---------------------------------------------------------------------------

  void ModularExtrap::SaveTo(TDirectory* dir) const
  {

    TDirectory* tmp = gDirectory;
    dir->cd();
    TObjString("ModularExtrap").Write("type");

    fEEextrap->SaveTo(dir->mkdir("EEextrap"));
    fEMextrap->SaveTo(dir->mkdir("EMextrap"));
    fMEextrap->SaveTo(dir->mkdir("MEextrap"));
    fMMextrap->SaveTo(dir->mkdir("MMextrap"));
    fEEAntiextrap->SaveTo(dir->mkdir("EEAntiextrap"));
    fEMAntiextrap->SaveTo(dir->mkdir("EMAntiextrap"));
    fMEAntiextrap->SaveTo(dir->mkdir("MEAntiextrap"));
    fMMAntiextrap->SaveTo(dir->mkdir("MMAntiextrap"));
    fMTextrap->SaveTo(dir->mkdir("MTextrap"));
    fETextrap->SaveTo(dir->mkdir("ETextrap"));
    fMTAntiextrap->SaveTo(dir->mkdir("MTAntiextrap"));
    fETAntiextrap->SaveTo(dir->mkdir("ETAntiextrap"));
    fNCextrap->SaveTo(dir->mkdir("NCextrap"));

    tmp->cd();

  }

  //---------------------------------------------------------------------------

  void  ModularExtrap::SavePlotsNue( TDirectory* dir, double potFD ) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();
    fMEextrap->SavePlots( dir->mkdir("MEextrap"), potFD );
    fMEAntiextrap->SavePlots( dir->mkdir("MEAntiextrap"), potFD );
    fEEextrap->SavePlots( dir->mkdir("EEextrap"), potFD );
    fMMextrap->SavePlots( dir->mkdir("MMextrap"), potFD );
    fNCextrap->SavePlots( dir->mkdir("NCextrap"), potFD );
    tmp->cd();
  }

  //---------------------------------------------------------------------------

  void  ModularExtrap::SavePlotsNumu( TDirectory* dir, double potFD ) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();
    fMMextrap->SavePlots( dir->mkdir("MMextrap"), potFD );
    fMMAntiextrap->SavePlots( dir->mkdir("MMAntiextrap"), potFD );
    tmp->cd();
  }

  //---------------------------------------------------------------------------

  std::unique_ptr<ModularExtrap> ModularExtrap::LoadFrom(TDirectory* dir)
  {

    std::unique_ptr<ModularExtrap> ret(new ModularExtrap);

    assert(dir->GetDirectory("EEextrap"));
    ret->fEEextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("EEextrap") );
 
    assert(dir->GetDirectory("EMextrap"));
    ret->fEMextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("EMextrap") );
 
    assert(dir->GetDirectory("MEextrap"));
    ret->fMEextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("MEextrap") );
 
    assert(dir->GetDirectory("MMextrap"));
    ret->fMMextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("MMextrap") );
 
    assert(dir->GetDirectory("EEAntiextrap"));
    ret->fEEAntiextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("EEAntiextrap") );
 
    assert(dir->GetDirectory("EMAntiextrap"));
    ret->fEMAntiextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("EMAntiextrap") );
 
    assert(dir->GetDirectory("MEAntiextrap"));
    ret->fMEAntiextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("MEAntiextrap") );
 
    assert(dir->GetDirectory("MMAntiextrap"));
    ret->fMMAntiextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("MMAntiextrap") );
 
    assert(dir->GetDirectory("NCextrap"));
    ret->fNCextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("NCextrap") );

    assert(dir->GetDirectory("MTextrap"));
    ret->fMTextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("MTextrap") );
 
    assert(dir->GetDirectory("ETextrap"));
    ret->fETextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("ETextrap") );
 
    assert(dir->GetDirectory("MTAntiextrap"));
    ret->fMTAntiextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("MTAntiextrap") );
 
    assert(dir->GetDirectory("ETAntiextrap"));
    ret->fETAntiextrap = ana::LoadFrom<ModularExtrapComponent>(
      dir->GetDirectory("ETAntiextrap") );
 
  return ret;

  }
 
  //---------------------------------------------------------------------------

  ModularExtrap::ModularExtrap(
    SpectrumLoaderBase& farMCnonswap,
    SpectrumLoaderBase& farMCnue,
    SpectrumLoaderBase& farMCnutau,
    SpectrumLoaderBase& farMCnc,
    const HistAxis& axis,
    const Cut& fdcut,
    const SystShifts& shiftMC,
    const Var& weight
  ) :

      // e -> e ----
      fEEextrap( new NoReweight(
        farMCnonswap, axis, fdcut, shiftMC, weight, kIsBeamNue && !kIsAntiNu)),
      fEEAntiextrap( new NoReweight(
        farMCnonswap, axis, fdcut, shiftMC, weight, kIsBeamNue && kIsAntiNu)),

      // mu -> mu  ----
      fMMextrap( new NoReweight(
        farMCnonswap, axis, fdcut, shiftMC, weight, kIsNumuCC && !kIsAntiNu)),
      fMMAntiextrap( new NoReweight(
        farMCnonswap, axis, fdcut, shiftMC, weight, kIsNumuCC && kIsAntiNu)),

      // mu -> e ----
      fMEextrap( new NoReweight(
        farMCnue, axis, fdcut, shiftMC, weight, kIsSig && !kIsAntiNu)),
      fMEAntiextrap( new NoReweight(
        farMCnue, axis, fdcut, shiftMC, weight, kIsSig && kIsAntiNu)),

      // e -> mu ----
      fEMextrap( new NoReweight(
        farMCnutau, axis, fdcut, shiftMC, weight, kIsNumuApp && !kIsAntiNu)),
      fEMAntiextrap( new NoReweight(
        farMCnutau, axis, fdcut, shiftMC, weight, kIsNumuApp && kIsAntiNu)),

      // NC -> NC ----
      fNCextrap( new NoReweight(
        farMCnc, axis, fdcut, shiftMC, weight, kIsNC)),

      // mu -> tau ----
      fMTextrap( new NoReweight(
        farMCnutau, axis, fdcut, shiftMC, weight, kIsTauFromMu && !kIsAntiNu)),
      fMTAntiextrap( new NoReweight(
        farMCnutau, axis, fdcut, shiftMC, weight, kIsTauFromMu && kIsAntiNu)),

      // e -> tau ----
      fETextrap( new NoReweight(
        farMCnue, axis, fdcut, shiftMC, weight, kIsTauFromE && !kIsAntiNu)),
      fETAntiextrap( new NoReweight(
        farMCnue, axis, fdcut, shiftMC, weight, kIsTauFromE && kIsAntiNu))

  {}

  //---------------------------------------------------------------------------

  OscillatableSpectrum ModularExtrap::NueSurvComponent()
    {return fEEextrap->Return();}

  OscillatableSpectrum ModularExtrap::AntiNueSurvComponent()
    {return fEEAntiextrap->Return();}

  OscillatableSpectrum ModularExtrap::NumuSurvComponent()
    {return fMMextrap->Return();}

  OscillatableSpectrum ModularExtrap::AntiNumuSurvComponent()
    {return fMMAntiextrap->Return();}

  OscillatableSpectrum ModularExtrap::NueAppComponent()
    {return fMEextrap->Return();}

  OscillatableSpectrum ModularExtrap::AntiNueAppComponent()
    {return fMEAntiextrap->Return();}

  OscillatableSpectrum ModularExtrap::NumuAppComponent()
    {return fEMextrap->Return();}

  OscillatableSpectrum ModularExtrap::AntiNumuAppComponent()
    {return fEMAntiextrap->Return();}

  Spectrum ModularExtrap::NCComponent()
    {return fNCextrap->Return().Unoscillated();}

  OscillatableSpectrum ModularExtrap::TauFromMuComponent()
    {return fMTextrap->Return();}

  OscillatableSpectrum ModularExtrap::AntiTauFromMuComponent()
    {return fMTAntiextrap->Return();}

  OscillatableSpectrum ModularExtrap::TauFromEComponent()
    {return fETextrap->Return();}

  OscillatableSpectrum ModularExtrap::AntiTauFromEComponent()
    {return fETAntiextrap->Return();}

}

#endif
