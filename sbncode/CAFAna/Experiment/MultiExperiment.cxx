#include "CAFAna/Experiment/MultiExperiment.h"

#include "CAFAna/Core/LoadFromFile.h"

#include "TDirectory.h"
#include "TObjString.h"

#include <cassert>

namespace ana
{
  //----------------------------------------------------------------------
  double MultiExperiment::ChiSq(osc::IOscCalcAdjustable* osc,
                                const SystShifts& syst) const
  {
    double ret = 0;
    for(unsigned int n = 0; n < fExpts.size(); ++n){
      // Don't waste time fiddling with the systematics if for sure there's
      // nothing to do.
      if(fSystCorrelations[n].empty()){
        ret += fExpts[n]->ChiSq(osc, syst);
      }
      else{
        // Make a local copy we're going to rewrite into the terms this
        // sub-experiment will accept.
        SystShifts localShifts = syst;
        for(auto it: fSystCorrelations[n]){
          // We're mapping prim -> sec
          const ISyst* prim = it.first;
          const ISyst* sec = it.second;
          if(syst.GetShift(prim) != 0){
            // sec can be unset, which means there's no representation needed
            // of prim in the sub-experiment.
            if(sec) localShifts.SetShift(sec, syst.GetShift(prim));
            // We've either translated or discarded prim, so drop it here.
            localShifts.SetShift(prim, 0);
          }
        }
        ret += fExpts[n]->ChiSq(osc, localShifts);
      }
    }
    return ret;
  }

  //----------------------------------------------------------------------
  void MultiExperiment::
  SetSystCorrelations(int idx,
                      const std::vector<std::pair<const ISyst*, const ISyst*>>& corrs)
  {
    // Sanity-check the mapping
    std::map<const ISyst*, const ISyst*> already;
    for(auto it: corrs){
      assert(it.first != it.second);

      // Don't worry if second element is null pointer
      if (!it.second) continue;
      if(already.find(it.second) == already.end()){
        already[it.second] = it.first;
      }
      else{
        std::cout << "MultiExperiment::SetSystCorrelations(): Warning!\n"
                  << "In experiment " << idx << " both "
                  << already[it.second]->ShortName() << " and "
                  << it.first->ShortName()
                  << " are configured to map to " << it.second->ShortName()
                  << ". That's probably not what you want." << std::endl;
      }
    }

    // Apply it
    fSystCorrelations[idx] = corrs;
  }

  //----------------------------------------------------------------------
  void MultiExperiment::SaveTo(TDirectory* dir) const
  {
    bool hasCorr = false;
    for(auto it: fSystCorrelations) if(!it.empty()) hasCorr = true;

    if(hasCorr){
      std::cerr << "Warning in MultiExperiment: systematic correlations are set and will not be serialized by this call to SaveTo(). You will have to re-set them once you load the experiment back in." << std::endl;
    }

    TDirectory* tmp = dir;

    dir->cd();
    TObjString("MultiExperiment").Write("type");

    for(unsigned int i = 0; i < fExpts.size(); ++i){
      fExpts[i]->SaveTo(dir->mkdir(TString::Format("expt%d", i)));
    }

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<MultiExperiment> MultiExperiment::LoadFrom(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    assert(ptag->GetString() == "MultiExperiment");

    std::vector<const IExperiment*> expts;

    for(int i = 0; ; ++i){
      TDirectory* subdir = dir->GetDirectory(TString::Format("expt%d", i));
      if(!subdir) break;

      expts.push_back(ana::LoadFrom<IExperiment>(subdir).release());
    }

    assert(!expts.empty());

    return std::unique_ptr<MultiExperiment>(new MultiExperiment(expts));
  }
}
