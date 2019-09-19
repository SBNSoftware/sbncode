#if 0
// Currently dead code

#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Extrap/ModularExtrapComponent.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"
#include "CAFAna/Decomp/IDecomp.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Utilities.h"

#include "TH1.h"
#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  //---------------------------------------------------------------------------
  // Definition to satisfy declaration in Core/LoadFromFile.h
  template<> std::unique_ptr<ModularExtrapComponent>
    LoadFrom<ModularExtrapComponent>(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);

    const TString tag = ptag->GetString();

    if(tag == "NoReweight")
      return NoReweight::LoadFrom(dir);
    if(tag == "TruthReweight")
      return TruthReweight::LoadFrom(dir);
    if(tag == "RecoReweight")
      return RecoReweight::LoadFrom(dir);

    std::cerr << "Unknown Component Extrapolation type '"
              << tag << "'" << std::endl;
    abort();
  }

  //---------------------------------------------------------------------------

  ModularExtrapComponent::DivByZeroCounter::~DivByZeroCounter()
  {
    if (fQuiet)
      return;

    if (this->fBins.size() == 0)
      return;

    std::cerr << "\nWARNING: There were attempted divisions by empty bins (for which a fallback to no reweighting was used)"
              << "\n         during extrapolation in the following bins (check your MC stats):"
              << std::endl;
    for (const auto & tuple : this->fBins)
    {
      std::cerr << "    Channel: "  << std::get<0>(tuple)
                << "    Location: "   << std::get<1>(tuple)
                << "    Bin Center: " << std::get<2>(tuple)
                << "\n" << std::endl;
    }

  }

  //---------------------------------------------------------------------------

  const OscillatableSpectrum& ModularExtrapComponent::Return() const
  {
    if (!fEvaluated)
    {
      fCache = Eval();
      fEvaluated = true;
    }
    return fCache;
  }

  //---------------------------------------------------------------------------

  Spectrum ModularExtrapComponent::GetDecompResult(
    const IDecomp& decomp,
    const DecompResult dr
  ){
    switch (dr){
      case DecompResult::nue     : return decomp.NueComponent();
      case DecompResult::numu    : return decomp.NumuComponent();
      case DecompResult::nuebar  : return decomp.AntiNueComponent();
      case DecompResult::numubar : return decomp.AntiNumuComponent();
      case DecompResult::NC      : return decomp.NCComponent();
    }
    assert( 0 && "Bad DecompResult" );
  }

  //---------------------------------------------------------------------------

  std::string ModularExtrapComponent::DRToString(
    const DecompResult dr
  ){
    switch (dr){
      case DecompResult::nue     : return "Nue";
      case DecompResult::numu    : return "Numu";
      case DecompResult::nuebar  : return "NueBar";
      case DecompResult::numubar : return "NumuBar";
      case DecompResult::NC      : return "NC";
    }
    assert( 0 && "Bad DecompResult" );
  }

  //---------------------------------------------------------------------------

  DecompResult ModularExtrapComponent::StringToDR(
    const std::string str
  ){
    if      (str=="Nue")     return DecompResult::nue;
    else if (str=="Numu")    return DecompResult::numu;
    else if (str=="NueBar")  return DecompResult::nuebar;
    else if (str=="NumuBar") return DecompResult::numubar;
    else if (str=="NC")      return DecompResult::NC;
    else assert( 0 && "Bad DecompResult String" );
  }

  //---------------------------------------------------------------------------

  Ratio ModularExtrapComponent::FormSmartRatio(
    const Spectrum& num,
    const Spectrum& denom,
    const std::string component,
    const std::string location,
    const Spectrum& mult
  ){

    DontAddDirectory guard;

    std::unique_ptr<TH1D> numh(num.ToTH1(1e20));
    std::unique_ptr<TH1D> denomh(denom.ToTH1(1e20));
    std::unique_ptr<TH1D> multh(mult.ToTH1(1e20));

    std::unique_ptr<TH1D> ratioh((TH1D*)numh->Clone( UniqueName().c_str() ));
    ratioh->Reset();

    int nbins = numh->GetNbinsX();
    assert( (nbins == denomh->GetNbinsX()) && "Bin Mismatch" );
    assert( (nbins == multh->GetNbinsX()) && "Bin Mismatch" );

    static DivByZeroCounter counter;
    for (int bin(0); bin<=nbins+1; ++bin)
    {
      if ( denomh->GetBinContent(bin) != 0 ){
        ratioh->SetBinContent( bin,   numh->GetBinContent(bin)
                                   / denomh->GetBinContent(bin) );
      } else {
        ratioh->SetBinContent( bin, 1. );
        if (    numh->GetBinContent(bin) != 0
             || multh->GetBinContent(bin) != 0 )
          counter.fBins.insert(std::make_tuple(component, location, denomh->GetBinCenter(bin)));
      }
    }

    return Ratio( ratioh.get() );
  }

  bool ModularExtrapComponent::fQuiet = false;

  //---------------------------------------------------------------------------

  NoReweight::NoReweight(
    SpectrumLoaderBase& loader,
    const HistAxis& axis,
    const Cut& fdcut,
    const SystShifts& shiftMC,
    const Var& weight,
    const Cut& flavors )
    : fRecoFD( loader, axis, fdcut && flavors, shiftMC, weight )
  {
  }

  OscillatableSpectrum NoReweight::Eval() const
    {return fRecoFD;}

  void NoReweight::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();
    TObjString("NoReweight").Write("type");
    fRecoFD.SaveTo(dir->mkdir("RecoFD"));
    tmp->cd();
  }

  std::unique_ptr<NoReweight>
    NoReweight::LoadFrom(TDirectory* dir)
  {
    assert(dir->GetDirectory("RecoFD"));
    return std::unique_ptr<NoReweight>(new NoReweight(
      *(OscillatableSpectrum::LoadFrom(dir->GetDirectory("RecoFD")).release())
    ));
  }

  //---------------------------------------------------------------------------

  TruthReweight::TruthReweight(
    SpectrumLoaderBase& ndloader,
    const HistAxis& axisFD,
    const HistAxis& axisND,
    const Cut& fdcut,
    const SystShifts& shiftMC,
    const Var& weight,
    std::string label,
    std::string latex,
    const Cut& ndcut,
    const IDecomp& decomposition,
    const DecompResult dr,
    const Cut& ndflavor,
    SpectrumLoaderBase& fdloader,
    const Cut& fdflavors
  )
    : fRecoToTrueND(ndloader, axisND, ndcut && ndflavor,  shiftMC, weight),
      fTrueToRecoFD(fdloader, axisFD, fdcut && fdflavors, shiftMC, weight),
      fDecomp(decomposition),
      fDecompRes(dr),
      fLabel(label),
      fLatex(latex)
  {}

  OscillatableSpectrum TruthReweight::Eval() const
  {

    //Copy to local variables because reweighting is in-place
    OscillatableSpectrum recoToTrueND(fRecoToTrueND);
    OscillatableSpectrum trueToRecoFD(fTrueToRecoFD);

    //Get ND data from Decomp
    Spectrum decompresult(GetDecompResult(fDecomp,fDecompRes));

    //Compute Data/MC Ratio in reco energy bins to get divide-by-zero warnings
    FormSmartRatio(
      decompresult, fRecoToTrueND.Unoscillated(),
      fLabel, "MC ND Reco",
      fRecoToTrueND.Unoscillated() );

    //ND Reco->True
    recoToTrueND.ReweightToRecoSpectrum( decompresult );

    //Compute Data/MC Ratio in true energy bins
    Ratio dataMCtrue = FormSmartRatio(
      recoToTrueND.TrueEnergy(), fRecoToTrueND.TrueEnergy(),
      fLabel, "MC ND Truth",
      fTrueToRecoFD.TrueEnergy() );

    // Multiply by Data/MC Ratio and add in FD truth information
    trueToRecoFD.ReweightToTrueSpectrum(   fTrueToRecoFD.TrueEnergy()
                                         * dataMCtrue );

    return trueToRecoFD;

  }

  void TruthReweight::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();
    TObjString("TruthReweight").Write("type");
    fRecoToTrueND.SaveTo(dir->mkdir("RecoToTrueND"));
    fTrueToRecoFD.SaveTo(dir->mkdir("TrueToRecoFD"));
    fDecomp.SaveTo(dir->mkdir("Decomp"));
    TObjString(DRToString(fDecompRes).c_str()).Write("DecompRes");
    TObjString(fLabel.c_str()).Write("Label");
    TObjString(fLatex.c_str()).Write("Latex");
    tmp->cd();
  }

  std::unique_ptr<TruthReweight>
    TruthReweight::LoadFrom(TDirectory* dir)
  {

    assert(dir->GetDirectory("RecoToTrueND"));
    assert(dir->GetDirectory("TrueToRecoFD"));
    assert(dir->GetDirectory("Decomp"));
    TObjString* dr = (TObjString*)dir->Get("DecompRes");
    assert(dr);
    TObjString* label = (TObjString*)dir->Get("Label");
    TObjString* latex = (TObjString*)dir->Get("Latex");
    assert(label);
    assert(latex);

    return std::unique_ptr<TruthReweight>(new TruthReweight(
      *(OscillatableSpectrum::LoadFrom(dir->GetDirectory("RecoToTrueND"))),
      *(OscillatableSpectrum::LoadFrom(dir->GetDirectory("TrueToRecoFD"))),
      *(ana::LoadFrom<IDecomp>(
        dir->GetDirectory("Decomp")).release()), //leaks!
      StringToDR(dr->GetString().Data()),
      label->GetString().Data(),
      latex->GetString().Data()
    ));

  }

  //---------------------------------------------------------------------------

  RecoReweight::RecoReweight(
    SpectrumLoaderBase& ndloader,
    const HistAxis& axis,
    const Cut& fdcut,
    const SystShifts& shiftMC,
    const Var& weight,
    std::string label,
    std::string latex,
    const Cut& ndcut,
    const IDecomp& decomposition,
    const DecompResult dr,
    const Cut& ndflavor,
    SpectrumLoaderBase& fdloader,
    const Cut& fdflavors
  )
    : fRecoND(ndloader, axis, ndcut && ndflavor, shiftMC, weight),
      fTrueToRecoFD(fdloader, axis, fdcut && fdflavors, shiftMC, weight),
      fDecomp(&decomposition),
      fDecompRes(dr),
      fLabel(label),
      fLatex(latex)
  {
  }

  OscillatableSpectrum RecoReweight::Eval() const
  {

    //Copy to local variable because reweighting is in-place
    OscillatableSpectrum result(fTrueToRecoFD);

    //Get ND data from Decomp
    Spectrum decompresult(GetDecompResult(*fDecomp,fDecompRes));

    //Compute Data/MC Ratio
    Ratio dataMC = FormSmartRatio(
      decompresult, fRecoND,
      fLabel, "MC ND Reco",
      fTrueToRecoFD.Unoscillated() );

    // Multiply by Data/MC Ratio and add in FD truth information
    result.ReweightToRecoSpectrum( fTrueToRecoFD.Unoscillated() * dataMC );

    return result;

  }

  void RecoReweight::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();
    TObjString("RecoReweight").Write("type");
    fRecoND.SaveTo(dir->mkdir("RecoND"));
    fTrueToRecoFD.SaveTo(dir->mkdir("TrueToRecoFD"));
    fDecomp->SaveTo(dir->mkdir("Decomp"));
    TObjString(DRToString(fDecompRes).c_str()).Write("DecompRes");
    TObjString(fLabel.c_str()).Write("Label");
    TObjString(fLatex.c_str()).Write("Latex");
    tmp->cd();
  }

  std::unique_ptr<RecoReweight>
    RecoReweight::LoadFrom(TDirectory* dir)
  {

    assert(dir->GetDirectory("RecoND"));
    assert(dir->GetDirectory("TrueToRecoFD"));
    assert(dir->GetDirectory("Decomp"));
    TObjString* dr = (TObjString*)dir->Get("DecompRes");
    assert(dr);
    TObjString* label = (TObjString*)dir->Get("Label");
    TObjString* latex = (TObjString*)dir->Get("Latex");
    assert(label);
    assert(latex);

    return std::unique_ptr<RecoReweight>(new RecoReweight(
      *(Spectrum::LoadFrom(dir->GetDirectory("RecoND"))),
      *(OscillatableSpectrum::LoadFrom(dir->GetDirectory("TrueToRecoFD"))),
      *(ana::LoadFrom<IDecomp>(
        dir->GetDirectory("Decomp")).release()), //leaks!
      StringToDR(dr->GetString().Data()),
      label->GetString().Data(),
      latex->GetString().Data()
    ));

  }

}

#endif
