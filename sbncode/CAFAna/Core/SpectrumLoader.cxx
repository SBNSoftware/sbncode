#include "CAFAna/Core/SpectrumLoader.h"

#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/ReweightableSpectrum.h"
#include "CAFAna/Core/SAMProjectSource.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Utilities.h"

#include "CAFAna/Core/GenieWeightList.h"

#include "StandardRecord/StandardRecord.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "TTreeFormula.h"

namespace ana
{
  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::string& wildcard, DataSource src, int max)
    : SpectrumLoaderBase(wildcard, src), max_entries(max)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::vector<std::string>& fnames,
                                 DataSource src, int max)
    : SpectrumLoaderBase(fnames, src), max_entries(max)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(DataSource src)
    : SpectrumLoaderBase(src)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader SpectrumLoader::FromSAMProject(const std::string& proj,
                                                DataSource src,
                                                int fileLimit)
  {
    SpectrumLoader ret;
    ret.fSource = src;
    ret.fWildcard = "project "+proj;
    ret.fFileSource = std::unique_ptr<IFileSource>(new SAMProjectSource(proj, fileLimit));
    return ret;
  }

  //----------------------------------------------------------------------
  SpectrumLoader::~SpectrumLoader()
  {
  }

  struct CompareByID
  {
    bool operator()(const Cut& a, const Cut& b)
    {
      return a.ID() < b.ID();
    }
  };

  //----------------------------------------------------------------------
  void SpectrumLoader::Go()
  {
    if(fGone){
      std::cerr << "Error: can only call Go() once on a SpectrumLoader" << std::endl;
      abort();
    }
    fGone = true;

    // Find all the unique cuts
    std::set<Cut, CompareByID> cuts;
    for(auto& shiftdef: fHistDefs)
      for(auto& cutdef: shiftdef.second)
        cuts.insert(cutdef.first);
    for(const Cut& cut: cuts) fAllCuts.push_back(cut);

    fLivetimeByCut.resize(fAllCuts.size());
    fPOTByCut.resize(fAllCuts.size());


    const int Nfiles = NFiles();

    Progress* prog = 0;

    caf::SRBranchRegistry::clear();

    int fileIdx = -1;
    while(TFile* f = GetNextFile()){
      ++fileIdx;

      if(Nfiles >= 0 && !prog) prog = new Progress(TString::Format("Filling %lu spectra from %d files matching '%s'", fHistDefs.TotalSize(), Nfiles, fWildcard.c_str()).Data());

      HandleFile(f, Nfiles == 1 ? prog : 0);

      if(Nfiles > 1 && prog) prog->SetProgress((fileIdx+1.)/Nfiles);
    } // end for fileIdx

    StoreExposures();

    if(prog){
      prog->Done();
      delete prog;
    }

    ReportExposures();

    fHistDefs.RemoveLoader(this);
    fHistDefs.Clear();
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleFile(TFile* f, Progress* prog)
  {
    assert(!f->IsZombie());

    // Test for flat (has extra trees) or nested cases.
    TDirectory* dir = f->Get("sbnana.reco") ? f : 0;
    TTree* tr = (TTree*)f->Get("sbnana");
    assert(tr);

    long n;
    caf::SRProxy sr(dir, tr, dir ? "sbnana" : "events", n, 0);

    //    FloatingExceptionOnNaN fpnan;

    long Nentries = tr->GetEntries();
    if (max_entries != 0 && max_entries < Nentries) Nentries = max_entries;

    for(n = 0; n < Nentries; ++n){
      if(!dir) tr->LoadTree(n); // nested mode

      HandleRecord(&sr);

      if(prog) prog->SetProgress(double(n)/Nentries);
    } // end for n
  }

  //----------------------------------------------------------------------
  /// Helper for \ref HandleRecord
  template<class T, class U> class CutVarCache
  {
  public:
    CutVarCache() : fVals(U::MaxID()+1), fValsSet(U::MaxID()+1, false) {}

    inline T Get(const U& var, const caf::SRProxy* sr)
    {
      const unsigned int id = var.ID();

      if(fValsSet[id]){
        return fVals[id];
      }
      else{
        const T val = var(sr);
        fVals[id] = val;
        fValsSet[id] = true;
        return val;
      }
    }

  protected:
    // Seems to be faster to do this than [unordered_]map
    std::vector<T> fVals;
    std::vector<bool> fValsSet;
  };

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleRecord(caf::SRProxy* sr)
  {
    // Some shifts only adjust the weight, so they're effectively nominal, but
    // aren't grouped with the other nominal histograms. Keep track of the
    // results for nominals in these caches to speed those systs up.
    CutVarCache<bool, Cut> nomCutCache;
    CutVarCache<double, Var> nomWeiCache;
    CutVarCache<double, Var> nomVarCache;

    for(auto& shiftdef: fHistDefs){
      const SystShifts& shift = shiftdef.first;

      // Need to provide a clean slate for each new set of systematic shifts to
      // work from. Copying the whole StandardRecord is pretty expensive, so
      // modify it in place and revert it afterwards.

      bool shifted = false;

      double systWeight = 1;
      // Can special-case nominal to not pay cost of Shift()
      if(!shift.IsNominal()){
        shift.Shift(sr, systWeight);
        // If there were only weighting systs applied then the cached nominal
        // values are still valid.
        shifted = caf::SRProxySystController::AnyShifted();
      }

      for(auto& cutdef: shiftdef.second){
        const Cut& cut = cutdef.first;

        const bool pass = shifted ? cut(sr) : nomCutCache.Get(cut, sr);
        // Cut failed, skip all the histograms that depended on it
        if(!pass) continue;

        for(auto& weidef: cutdef.second){
          const Var& weivar = weidef.first;

          double wei = shifted ? weivar(sr) : nomWeiCache.Get(weivar, sr);

          wei *= systWeight;
          if(wei == 0) continue;

          for(auto& vardef: weidef.second){
            if(vardef.first.IsMulti()){
              for(double val: vardef.first.GetMultiVar()(sr)){
                for(Spectrum* s: vardef.second.spects)
                  s->Fill(val, wei);
              }
              continue;
            }

            const Var& var = vardef.first.GetVar();

            const double val = shifted ? var(sr) : nomVarCache.Get(var, sr);

            if(std::isnan(val) || std::isinf(val)){
              std::cerr << "Warning: Bad value: " << val
                        << " returned from a Var. The input variable(s) could "
                        << "be NaN in the CAF, or perhaps your "
                        << "Var code computed 0/0?";
              std::cout << " Not filling into this histogram for this slice." << std::endl;
              continue;
            }

            for(Spectrum* s: vardef.second.spects) s->Fill(val, wei);

            for(ReweightableSpectrum* rw: vardef.second.rwSpects){
              const double yval = rw->ReweightVar()(sr);

              if(std::isnan(yval) || std::isinf(yval)){
                std::cerr << "Warning: Bad value: " << yval
                          << " for reweighting Var";
                std::cout << ". Not filling into histogram." << std::endl;
                continue;
              }

              // TODO: ignoring events with no true neutrino etc
              if(yval != 0) rw->fHist->Fill(val, yval, wei);
            } // end for rw
          } // end for vardef
        } // end for weidef
      } // end for cutdef

       
      // Return StandardRecord to its unshifted form ready for the next
      // histogram.
      caf::SRProxySystController::ResetSysts();
    } // end for shiftdef
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::ReportExposures()
  {
    // The POT member variables we use here were filled as part of
    // SpectrumLoaderBase::GetNextFile() as we looped through the input files.

    // Let's just assume no-one is using the Cut::POT() function yet, so this
    // printout remains relevant...

    std::cout << fPOT << " POT" << std::endl;
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::AccumulateExposures(const caf::SRSpill* spill)
  {
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::StoreExposures()
  {
    for(auto& shiftdef: fHistDefs){
      for(auto& cutdef: shiftdef.second){
        for(auto& weidef: cutdef.second){
          for(auto& vardef: weidef.second){
            for(Spectrum* s: vardef.second.spects) s->fPOT += fPOT;
            for(ReweightableSpectrum* rw: vardef.second.rwSpects) rw->fPOT += fPOT;
          }
        }
      }
    }
  }
} // namespace
