#include "CAFAna/Prediction/PredictionLinFit.h"

#include "CAFAna/Core/HistCache.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/SystRegistry.h"

#include "OscLib/IOscCalculator.h"

#include "TDecompLU.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TObjString.h"

namespace ana
{
  // --------------------------------------------------------------------------
  PredictionLinFit::PredictionLinFit(const std::vector<const ISyst*>& systs,
                                     const IPrediction* pnom,
                                     const std::vector<std::pair<SystShifts, const IPrediction*>>& univs)
    : fSysts(systs), fNom(pnom), fUnivs(univs)
  {
    osc::NoOscillations calc; // TODO pass this in?

    const Spectrum snom = pnom->Predict(&calc);
    TH1D* hnom = snom.ToTH1(1); // only need this to count the bins :(
    const int Nbins = hnom->GetNbinsX();
    HistCache::Delete(hnom);

    std::vector<std::vector<double>> coords(univs.size());

    for(unsigned int univIdx = 0; univIdx < univs.size(); ++univIdx){
      coords[univIdx].reserve(systs.size());
      for(const ISyst* s: systs) coords[univIdx].push_back(univs[univIdx].first.GetShift(s));
    }

    TMatrixD M(systs.size(), systs.size());
    for(unsigned int univIdx = 0; univIdx < univs.size(); ++univIdx){
      for(unsigned int i = 0; i < systs.size(); ++i){
        for(unsigned int j = 0; j < systs.size(); ++j){
          M(i, j) += coords[univIdx][i] * coords[univIdx][j];
        }
      }
    }

    TDecompLU decomp(M);

    std::vector<TVectorD> vs(Nbins+2, TVectorD(systs.size()));

    for(unsigned int univIdx = 0; univIdx < univs.size(); ++univIdx){
      const Ratio r(univs[univIdx].second->Predict(&calc), snom);
      TH1D* hr = r.ToTH1();

      for(int binIdx = 0; binIdx < Nbins+2; ++binIdx){
        const double d = hr->GetBinContent(binIdx);

        for(unsigned int i = 0; i < systs.size(); ++i){
          vs[binIdx](i) += coords[univIdx][i] * log(d);
        }
      }

      HistCache::Delete(hr);
    }

    for(TVectorD& v: vs){
      decomp.Solve(v);
      fCoeffs.push_back(v);
    }
  }

  // --------------------------------------------------------------------------
  PredictionLinFit::~PredictionLinFit()
  {
  }

  // --------------------------------------------------------------------------
  Spectrum PredictionLinFit::Predict(osc::IOscCalculator* calc) const
  {
    return fNom->Predict(calc);
  }

  // --------------------------------------------------------------------------
  Spectrum PredictionLinFit::PredictSyst(osc::IOscCalculator* calc,
                                         const SystShifts& syst) const
  {
    return GetRatio(syst) * Predict(calc);
  }

  // --------------------------------------------------------------------------
  Spectrum PredictionLinFit::PredictComponent(osc::IOscCalculator* calc,
                                              Flavors::Flavors_t flav,
                                              Current::Current_t curr,
                                              Sign::Sign_t sign) const
  {
    return fNom->PredictComponent(calc, flav, curr, sign);
  }

  // --------------------------------------------------------------------------
  Spectrum PredictionLinFit::PredictComponentSyst(osc::IOscCalculator* calc,
                                                  const SystShifts& syst,
                                                  Flavors::Flavors_t flav,
                                                  Current::Current_t curr,
                                                  Sign::Sign_t sign) const
  {
    return GetRatio(syst) * PredictComponent(calc, flav, curr, sign);
  }

  // --------------------------------------------------------------------------
  Ratio PredictionLinFit::GetRatio(const SystShifts& shift) const
  {
    // To get correct binning
    TH1D* hret = fNom->PredictUnoscillated().ToTH1(1);
    hret->Reset();

    for(unsigned int binIdx = 0; binIdx < fCoeffs.size(); ++binIdx){
      double factor = 0;
      // TODO consider looping through the set entries in 'shift', rather than
      // fSysts
      for(unsigned int systIdx = 0; systIdx < fSysts.size(); ++systIdx){
        const double x = shift.GetShift(fSysts[systIdx]);
        factor += fCoeffs[binIdx][systIdx] * x;
      }
      hret->SetBinContent(binIdx, exp(factor));
    }

    const Ratio ret(hret);
    HistCache::Delete(hret);
    return ret;
  }

  //----------------------------------------------------------------------
  void PredictionLinFit::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();
    TObjString("PredictionLinFit").Write("type");

    fNom->SaveTo(dir->mkdir("nom"));

    for(unsigned int univIdx = 0; univIdx < fUnivs.size(); ++univIdx){
      TDirectory* ud = dir->mkdir(TString::Format("univ_%d", univIdx).Data());
      fUnivs[univIdx].first.SaveTo(ud->mkdir("shift"));
      fUnivs[univIdx].second->SaveTo(ud->mkdir("pred"));
    } // end for it

    if(!fSysts.empty()){
      TH1F hSystNames("syst_names", ";Syst names", fSysts.size(), 0, fSysts.size());
      int binIdx = 1;
      for(auto it: fSysts){
        hSystNames.GetXaxis()->SetBinLabel(binIdx++, it->ShortName().c_str());
      }
      hSystNames.Write("syst_names");
    }

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<PredictionLinFit> PredictionLinFit::LoadFrom(TDirectory* dir)
  {
    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "PredictionLinFit");

    std::unique_ptr<IPrediction> nom = ana::LoadFrom<IPrediction>(dir->GetDirectory("nom"));

    std::vector<const ISyst*> systs;
    TH1* hSystNames = (TH1*)dir->Get("syst_names");
    if(hSystNames){
      for(int systIdx = 0; systIdx < hSystNames->GetNbinsX(); ++systIdx){
        systs.push_back(SystRegistry::ShortNameToSyst(hSystNames->GetXaxis()->GetBinLabel(systIdx+1)));
      }
    }

    std::vector<std::pair<SystShifts, const IPrediction*>> univs;

    for(int univIdx = 0; ; ++univIdx){
      TDirectory* ud = dir->GetDirectory(TString::Format("univ_%d", univIdx).Data());
      if(!ud) break; // out of universes

      univs.emplace_back(*ana::LoadFrom<SystShifts>(ud->GetDirectory("shift")),
                         ana::LoadFrom<IPrediction>(ud->GetDirectory("pred")).release());
    }

    // TODO think about memory management
    return std::make_unique<PredictionLinFit>(systs, nom.release(), univs);
  }

}
