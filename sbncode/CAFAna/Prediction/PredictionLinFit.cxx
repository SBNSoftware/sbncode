#include "CAFAna/Prediction/PredictionLinFit.h"

#include "CAFAna/Core/HistCache.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/MathUtil.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/SystRegistry.h"

#include "CAFAna/Prediction/PredictionGenerator.h"

#include "CAFAna/Prediction/IncrementalCholeskyDecomp.h"

#include "CAFAna/Systs/SBNWeightSysts.h"
#include "CAFAna/Systs/UniverseOracle.h"

#include "OscLib/IOscCalculator.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH2.h"
#include "TObjString.h"

namespace ana
{
  // --------------------------------------------------------------------------
  PredictionLinFit::PredictionLinFit(const std::vector<const ISyst*>& systs,
                                     const IPrediction* pnom,
                                     const std::vector<std::pair<SystShifts, const IPrediction*>>& univs)
    : fSysts(systs), fNom(pnom), fUnivs(univs)
  {
  }

  // --------------------------------------------------------------------------
  PredictionLinFit::PredictionLinFit(const std::vector<const ISyst*>& systs,
                                     const IPredictionGenerator& predGen,
                                     Loaders& loaders,
                                     int nUniv)
    : fSysts(systs)
  {
    fNom = predGen.Generate(loaders, SystShifts::Nominal()).release();

    // We really want to apply a weight - but PredictionGenerator is set up to
    // take a SystShifts, so...
    class DummyUnivWeightSyst: public ISyst
    {
    public:
      DummyUnivWeightSyst(const std::vector<const ISyst*> systs, int univIdx)
        : ISyst(UniqueName(), UniqueName()),
          fVar(GetUniverseWeight(systs, univIdx))
      {
      }

      void Shift(double sigma, caf::SRProxy* sr, double& weight) const override
      {
        weight *= fVar(sr);
      }
    protected:
      Var fVar;
    };

    for(int univIdx = 0; univIdx < nUniv; ++univIdx){
      // This one is an accurate label for what the shift is
      const SystShifts s1 = UniverseOracle::Instance().ShiftsForSysts(systs, nUniv)[univIdx];
      // But this is the sane way to compute its effect
      const SystShifts s2(new DummyUnivWeightSyst(systs, univIdx), 1);
      fUnivs.emplace_back(s1, predGen.Generate(loaders, s2).release());
    }
  }

  // --------------------------------------------------------------------------
  PredictionLinFit::~PredictionLinFit()
  {
  }

  // --------------------------------------------------------------------------
  void PredictionLinFit::InitFits() const
  {
    if(!fCoeffs.empty()) return;

    std::vector<std::vector<double>> coords;
    coords.reserve(fUnivs.size());
    for(const auto& it: fUnivs) coords.push_back(GetCoords(it.first));

    const unsigned int N = coords[0].size();

    std::vector<std::vector<double>> M(N, std::vector<double>(N));

    for(unsigned int univIdx = 0; univIdx < fUnivs.size(); ++univIdx){
      for(unsigned int i = 0; i < N; ++i){
        for(unsigned int j = 0; j < N; ++j){
          M[i][j] += coords[univIdx][i] * coords[univIdx][j];
        }
      }
    }

    osc::NoOscillations calc; // TODO class member?

    const Spectrum snom = fNom->Predict(&calc);

    TH1D* hnom = snom.ToTH1(1); // only need this to count the bins :(
    const int Nbins = hnom->GetNbinsX();
    HistCache::Delete(hnom);

    // The data (ratios) that we're trying to fit
    std::vector<std::vector<double>> ds(Nbins+2, std::vector<double>(fUnivs.size()));

    for(unsigned int univIdx = 0; univIdx < fUnivs.size(); ++univIdx){
      const Ratio r(fUnivs[univIdx].second->Predict(&calc), snom);
      TH1D* hr = r.ToTH1();

      for(int binIdx = 0; binIdx < Nbins+2; ++binIdx){
        ds[binIdx][univIdx] = log(hr->GetBinContent(binIdx));
      }

      HistCache::Delete(hr);
    }

    for(int binIdx = 0; binIdx < Nbins+2; ++binIdx){
      fCoeffs.push_back(InitFitsBin(M, ds[binIdx], coords));
    }
  }

  // --------------------------------------------------------------------------
  std::vector<double> PredictionLinFit::
  InitFitsBin(const std::vector<std::vector<double>>& M,
              const std::vector<double>& ds,
              const std::vector<std::vector<double>>& coords) const
  {
    const unsigned int N = M.size();

    std::vector<double> v(N);

    for(unsigned int univIdx = 0; univIdx < fUnivs.size(); ++univIdx){
      for(unsigned int i = 0; i < N; ++i){
        v[i] += coords[univIdx][i] * ds[univIdx];
      }
    }

    std::cout << "-----" << std::endl;
    double best_mse = std::numeric_limits<double>::infinity();
    std::vector<bool> already(N);

    IncrementalCholeskyDecomp icd;

    std::vector<double> vsub;

    std::vector<int> vars;

    std::vector<std::vector<double>> coords_sub(fUnivs.size());

    for(unsigned int matSize = 1; matSize <= 50/*std::min(N, fUnivs.size())*/; ++matSize){
      vsub.push_back(0);
      vars.push_back(-1);
      icd.Extend();
      for(std::vector<double>& c: coords_sub) c.push_back(0);

      int bestVarIdx = -1;
      for(unsigned int varIdx = 0; varIdx < N; ++varIdx){
        if(already[varIdx]) continue;

        // Provisionally add this var to the list
        vars.back() = varIdx;

        icd.SetLastRow(M[varIdx], vars);

        vsub.back() = v[varIdx];

        const std::vector<double> a = icd.Solve(vsub);

        // Mean squared error
        double mse = 0;
        for(unsigned int univIdx = 0; univIdx < fUnivs.size(); ++univIdx){
          const double target = ds[univIdx];

          std::vector<double>& c = coords_sub[univIdx];
          c.back() = coords[univIdx][varIdx];

          double estimate = 0;
          for(unsigned int i = 0; i < matSize; ++i){
            estimate += a[i] * c[i];
          }

          mse += util::sqr(estimate-target);
        }

        mse /= fUnivs.size();

        //        std::cout << "    " << fSysts[varIdx]->ShortName() << " " << sqrt(mse) << std::endl;

        if(mse < best_mse){
          best_mse = mse;
          bestVarIdx = varIdx;
        }
      } // end for varIdx

      if(bestVarIdx == -1){
        std::cout << "No var improved. Done" << std::endl;
        vars.pop_back();
        break;
      }

      std::cout << "  " << matSize << ": best var " << bestVarIdx/*fSysts[bestVarIdx]->ShortName()*/ << " " << sqrt(best_mse) << std::endl;
      vars.back() = bestVarIdx;

      // And set the correct row back into the matrix and vector
      vsub.back() = v[bestVarIdx];

      for(unsigned int univIdx = 0; univIdx < coords.size(); ++univIdx){
        coords_sub[univIdx].back() = coords[univIdx][bestVarIdx];
      }

      icd.SetLastRow(M[bestVarIdx], vars);

      already[bestVarIdx] = true;
    } // end for matSize

    std::cout << "VARS:";
    for(int v: vars) std::cout << " " << v;
    std::cout << " MSE: " << sqrt(best_mse) << std::endl;
    std::cout << std::endl;

    // Solve one last time
    vsub = icd.Solve(vsub);

    // And package up
    std::vector<double> ret(N);
    for(unsigned int i = 0; i < vars.size(); ++i) ret[vars[i]] = vsub[i];
    return ret;
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
  std::vector<double> PredictionLinFit::GetCoords(const SystShifts& shift) const
  {
    const unsigned int N = fSysts.size();

    std::vector<double> coords;
    coords.reserve(N + (N*(N+1))/2);

    // Add all the linear terms
    for(const ISyst* s: fSysts) coords.push_back(shift.GetShift(s));

    // Now add all the quadratic and cross terms
    for(unsigned int i = 0; i < N; ++i){
      for(unsigned int j = i; j < N; ++j){
        coords.push_back(coords[i] * coords[j]);
      }
    }

    return coords;
  }

  // --------------------------------------------------------------------------
  Ratio PredictionLinFit::GetRatio(const SystShifts& shift) const
  {
    if(fCoeffs.empty()) InitFits();

    // To get correct binning
    TH1D* hret = fNom->PredictUnoscillated().ToTH1(1);
    hret->Reset();

    const std::vector<double> coords = GetCoords(shift);
    const unsigned int N = coords.size();

    for(unsigned int binIdx = 0; binIdx < fCoeffs.size(); ++binIdx){
      double factor = 0;
      for(unsigned int i = 0; i < N; ++i) factor += fCoeffs[binIdx][i] * coords[i];

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
  void PredictionLinFit::DebugPlot(const ISyst* syst,
                                   osc::IOscCalculator* calc,
                                   Flavors::Flavors_t flav,
                                   Current::Current_t curr,
                                   Sign::Sign_t sign) const
  {
    DontAddDirectory guard;

    if(fCoeffs.empty()) InitFits();

    std::unique_ptr<TH1> nom(fNom->PredictComponent(calc, flav, curr, sign).ToTH1(18e20));
    const int nbins = nom->GetNbinsX();

    TGraph* curves[nbins];

    for(int i = 0; i <= 80; ++i){
      const double x = .1*i-4;
      const SystShifts ss(syst, x);
      std::unique_ptr<TH1> h(PredictComponentSyst(calc, ss, flav, curr, sign).ToTH1(18e20));

      for(int bin = 0; bin < nbins; ++bin){
        if(i == 0){
          curves[bin] = new TGraph;
        }

        const double ratio = h->GetBinContent(bin+1)/nom->GetBinContent(bin+1);

        if(!std::isnan(ratio)) curves[bin]->SetPoint(curves[bin]->GetN(), x, ratio);
        else curves[bin]->SetPoint(curves[bin]->GetN(), x, 1);
      } // end for bin
    } // end for i (x)

    int nx = int(sqrt(nbins));
    int ny = int(sqrt(nbins));
    if(nx*ny < nbins) ++nx;
    if(nx*ny < nbins) ++ny;

    TCanvas* c = new TCanvas;
    c->Divide(nx, ny);

    for(int bin = 0; bin < nbins; ++bin){
      c->cd(bin+1);
      (new TH2F(UniqueName().c_str(),
                TString::Format("%s %g < %s < %g;Shift;Ratio",
                                syst->ShortName().c_str(),
                                nom->GetXaxis()->GetBinLowEdge(bin+1),
                                nom->GetXaxis()->GetTitle(),
                                nom->GetXaxis()->GetBinUpEdge(bin+1)),
                100, -4, +4, 100, .5, 1.5))->Draw();
      curves[bin]->Draw("l same");
    } // end for bin

    c->cd(0);
  }

  //----------------------------------------------------------------------
  void PredictionLinFit::DebugPlots(osc::IOscCalculator* calc,
				    const std::string& savePattern,
				    Flavors::Flavors_t flav,
				    Current::Current_t curr,
				    Sign::Sign_t sign) const
  {
    const bool save = !savePattern.empty();
    const bool multiFile = savePattern.find("%s") != std::string::npos;

    if(save && !multiFile){
      new TCanvas;
      gPad->Print((savePattern+"[").c_str());
    }

    for(const ISyst* s: fSysts){
      DebugPlot(s, calc, flav, curr, sign);

      if(save){
        if(multiFile){
          gPad->Print(TString::Format(savePattern.c_str(), s->ShortName().c_str()).Data());
        }
        else{
          gPad->Print(TString::Format(savePattern.c_str(), s->ShortName()).Data());
        }
      }
    } // end for s

    if(save && !multiFile){
      gPad->Print((savePattern+"]").c_str());
    }
  }

  //----------------------------------------------------------------------
  void PredictionLinFit::DebugPlotColz(const ISyst* syst,
                                       osc::IOscCalculator* calc,
                                       Flavors::Flavors_t flav,
                                       Current::Current_t curr,
                                       Sign::Sign_t sign) const
  {
    InitFits();

    std::unique_ptr<TH1> nom(fNom->PredictComponent(calc, flav, curr, sign).ToTH1(18e20));
    const int nbins = nom->GetNbinsX();

    TH2* h2 = new TH2F("", (syst->LatexName()+";;#sigma").c_str(),
                       nbins, nom->GetXaxis()->GetXmin(), nom->GetXaxis()->GetXmax(),
                       80, -4, +4);
    h2->GetXaxis()->SetTitle(nom->GetXaxis()->GetTitle());

    for(int i = 1; i <= 80; ++i){
      const double y = h2->GetYaxis()->GetBinCenter(i);
      const SystShifts ss(syst, y);
      std::unique_ptr<TH1> h(PredictComponentSyst(calc, ss, flav, curr, sign).ToTH1(18e20));

      for(int bin = 0; bin < nbins; ++bin){
        const double ratio = h->GetBinContent(bin+1)/nom->GetBinContent(bin+1);

        if(!isnan(ratio) && !isinf(ratio))
          h2->Fill(h2->GetXaxis()->GetBinCenter(bin), y, ratio);
      } // end for bin
    } // end for i (x)

    h2->Draw("colz");
    h2->SetMinimum(0.5);
    h2->SetMaximum(1.5);
  }

  //----------------------------------------------------------------------
  void PredictionLinFit::DebugPlotsColz(osc::IOscCalculator* calc,
                                        const std::string& savePattern,
                                        Flavors::Flavors_t flav,
                                        Current::Current_t curr,
                                        Sign::Sign_t sign) const
  {
    const bool save = !savePattern.empty();
    const bool multiFile = savePattern.find("%s") != std::string::npos;

    if(save && !multiFile){
      new TCanvas;
      gPad->Print((savePattern+"[").c_str());
    }

    for(const ISyst* s: fSysts){
      new TCanvas;
      DebugPlotColz(s, calc, flav, curr, sign);

      if(save){
        if(multiFile){
          gPad->Print(TString::Format(savePattern.c_str(), s->ShortName().c_str()).Data());
        }
        else{
          gPad->Print(savePattern.c_str());
        }
      }
    } // end for it

    if(save && !multiFile){
      gPad->Print((savePattern+"]").c_str());
    }
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
