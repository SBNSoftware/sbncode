#include "CAFAna/Experiment/RatioExperiment.h"

#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/HistCache.h"

#include "CAFAna/Prediction/IPrediction.h"

#include "TDirectory.h"
#include "TObjString.h"

#include "TH1.h"

#include <cassert>

namespace ana
{
  //----------------------------------------------------------------------
  double RatioExperiment::ChiSq(osc::IOscCalculatorAdjustable* osc,
                                const SystShifts& syst) const
  {
    osc->SetL(kBaselineSBND);
    Spectrum predND = fPredND->PredictSyst(osc, syst);
    osc->SetL(kBaselineIcarus);
    Spectrum predFD = fPredFD->PredictSyst(osc, syst);
    
    predFD *= Ratio(fSpectND, predND);

    TH1D* hpred = predFD.ToTH1(fSpectFD.POT());
    TH1D* hdata = fSpectFD.ToTH1(fSpectFD.POT());

    const double ret = LogLikelihood(hpred, hdata);

    HistCache::Delete(hpred);
    HistCache::Delete(hdata);

    return ret;
  }

  //----------------------------------------------------------------------
  // void RatioExperiment::
  // Derivative(osc::IOscCalculator* calc,
  //            const SystShifts& shift,
  //            std::unordered_map<const ISyst*, double>& dchi) const
  // {
  //   // Each one should sum into the total so far
  //   for(unsigned int n = 0; n < fExpts.size(); ++n){
  //     if(fSystCorrelations[n].empty()){
  //       // If there are no adjustments needed it's easy
  //       fExpts[n]->Derivative(calc, shift, dchi);
  //       if(dchi.empty()) return;
  //       continue;
  //     }

  //     auto dchiLocal = dchi;
  //     for(auto& it: dchiLocal) it.second = 0;
  //     SystShifts localShifts = shift;
  //     std::unordered_map<const ISyst*, const ISyst*> reverseMap;
  //     for(auto it: fSystCorrelations[n]){
  //       // We're mapping prim -> sec
  //       const ISyst* prim = it.first;
  //       const ISyst* sec = it.second;
  //       if(dchi.count(prim)){
  //         dchiLocal.erase(dchiLocal.find(prim));
  //         if(sec){
  //           dchiLocal.emplace(sec, 0);
  //           reverseMap[sec] = prim;
  //         }
  //       }
  //       if(sec) localShifts.SetShift(sec, shift.GetShift(prim));
  //       // We've either translated or discarded prim, so drop it here.
  //       localShifts.SetShift(prim, 0);
  //     }

  //     fExpts[n]->Derivative(calc, localShifts, dchiLocal);

  //     // One of the components doesn't support derivatives, don't bother asking
  //     // the rest.
  //     if(dchiLocal.empty()){
  //       dchi.clear();
  //       return;
  //     }

  //     // And translate back
  //     for(auto it: dchiLocal){
  //       if(reverseMap.count(it.first) > 0){
  //         dchi[reverseMap[it.first]] += it.second;
  //       }
  //       else{
  //         dchi[it.first] += it.second;
  //       }
  //     }
  //   }
  // }

  // //----------------------------------------------------------------------
  // void RatioExperiment::SaveTo(TDirectory* dir) const
  // {
  //   bool hasCorr = false;
  //   for(auto it: fSystCorrelations) if(!it.empty()) hasCorr = true;

  //   if(hasCorr){
  //     std::cerr << "Warning in RatioExperiment: systematic correlations are set and will not be serialized by this call to SaveTo(). You will have to re-set them once you load the experiment back in." << std::endl;
  //   }

  //   TDirectory* tmp = dir;

  //   dir->cd();
  //   TObjString("RatioExperiment").Write("type");

  //   for(unsigned int i = 0; i < fExpts.size(); ++i){
  //     fExpts[i]->SaveTo(dir->mkdir(TString::Format("expt%d", i)));
  //   }

  //   tmp->cd();
  // }

  // //----------------------------------------------------------------------
  // std::unique_ptr<RatioExperiment> RatioExperiment::LoadFrom(TDirectory* dir)
  // {
  //   TObjString* ptag = (TObjString*)dir->Get("type");
  //   assert(ptag);
  //   assert(ptag->GetString() == "RatioExperiment");

  //   std::vector<const IExperiment*> expts;

  //   for(int i = 0; ; ++i){
  //     TDirectory* subdir = dir->GetDirectory(TString::Format("expt%d", i));
  //     if(!subdir) break;

  //     expts.push_back(ana::LoadFrom<IExperiment>(subdir).release());
  //   }

  //   assert(!expts.empty());

  //   return std::unique_ptr<RatioExperiment>(new RatioExperiment(expts));
  // }
}
