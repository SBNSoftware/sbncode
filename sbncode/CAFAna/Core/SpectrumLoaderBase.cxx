#include "CAFAna/Core/SpectrumLoaderBase.h"

#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/ReweightableSpectrum.h"
#include "CAFAna/Core/SAMQuerySource.h"
#include "CAFAna/Core/SAMProjectSource.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Utilities.h"
#include "CAFAna/Core/WildcardSource.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include "ifdh.h"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include "boost/algorithm/string.hpp"

namespace ana
{

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectList::~SpectList()
  {
    // We don't seem to be able to do this - maybe SpectList gets copied?
    // Let's just bite the bullet and leak a few pointers. This all works
    // completely differently in NOvA CAFAna with refactored SpectrumLoader.
    //
    // Clean up the memory allocated for the pointers themselves
    //    for(Spectrum** s: spects) delete *s;
    //    for(auto rv: rwSpects) delete *rv.first;
  }

  // Work around ReweightableSpectrum's friend requirements
  struct ReweightableSpectrumSink
  {
    static void AddLoader(ReweightableSpectrum** rw){(*rw)->AddLoader(rw);}
    static void RemoveLoader(ReweightableSpectrum** rw){(*rw)->RemoveLoader(rw);}
  };

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::RemoveLoader(SpectrumLoaderBase* l)
  {
    for(Spectrum** s: spects) if(*s) (*s)->RemoveLoader(s);
    for(auto rv: rwSpects) if(*rv.first) ReweightableSpectrumSink::RemoveLoader(rv.first);
  }

  //----------------------------------------------------------------------
  size_t SpectrumLoaderBase::SpectList::TotalSize() const
  {
    return spects.size() + rwSpects.size();
  }

  //----------------------------------------------------------------------
  template<class T, class U> U& SpectrumLoaderBase::IDMap<T, U>::
  operator[](const T& key)
  {
    for(auto& it: fElems){
      if(it.first.ID() == key.ID()) return it.second;
    }
    fElems.push_back(std::make_pair(key, U()));
    return fElems.back().second;
  }

  //----------------------------------------------------------------------
  template<class T, class U> void SpectrumLoaderBase::IDMap<T, U>::
  RemoveLoader(SpectrumLoaderBase* l)
  {
    for(auto& it: fElems) it.second.RemoveLoader(l);
  }

  //----------------------------------------------------------------------
  template<class T, class U> void SpectrumLoaderBase::IDMap<T, U>::Clear()
  {
    fElems.clear();
  }

  //----------------------------------------------------------------------
  template<class T, class U> size_t SpectrumLoaderBase::IDMap<T, U>::
  TotalSize()
  {
    size_t ret = 0;
    for(auto& it: fElems) ret += it.second.TotalSize();
    return ret;
  }

  // Start of SpectrumLoaderBase proper

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase()
    : fGone(false), fPOT(0)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(const std::string& wildcard)
    : SpectrumLoaderBase()
  {
    fWildcard = wildcard;
    fFileSource = std::unique_ptr<IFileSource>(WildcardOrSAMQuery(wildcard));
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(const std::vector<std::string>& fnames) : SpectrumLoaderBase()
  {
    fWildcard = "file list";
    fFileSource = std::unique_ptr<IFileSource>(new FileListSource(fnames));

    assert(!fnames.empty());
    std::cout << "Loading from " << fnames.size() << " files" << std::endl;
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::~SpectrumLoaderBase()
  {
    fHistDefs.RemoveLoader(this);
    fHistDefs.Clear();
  }

  //----------------------------------------------------------------------
  IFileSource* SpectrumLoaderBase::
  WildcardOrSAMQuery(const std::string& str) const
  {
    int stride = -1;
    int offset = -1;
    if(getenv("CAFANA_STRIDE")){
      stride = atoi(getenv("CAFANA_STRIDE"));
      if(stride > 1 && getenv("CAFANA_OFFSET")){
        offset = atoi(getenv("CAFANA_OFFSET"));
      }
    }

    // stat() blows up on strings with spaces
    if(str.find(' ') == std::string::npos){
      WildcardSource* ret = new WildcardSource(str, stride, offset);
      if(ret->NFiles() > 0) return ret;
      delete ret;
    }

    // Maybe this the name of a SAM project?
    {
      IFDHSilent silent; // the usual case is for this to fail
      ifdh_ns::ifdh i;
      // findProject always gives back an address just by gluing bits together.
      // (a tad annoying, because it _does_ go and look for the project, and
      // even would print its 'didn't-find-this-project' error out to stderr if
      // not for the IFDHSilent, but without scraping its stderr there's no way
      // to know whether the project is there or not -- you still get the URL.)
      // however, the WebAPI call looking for the /status will return a 404 if
      // the project doesn't exist.  (suggested by Robert I. in
      // INC000000925362)
      try{
        ifdh_util_ns::WebAPI webapi(i.findProject(str, Experiment()) +  "/status");
        return new SAMProjectSource(str);
      }
      catch(ifdh_util_ns::WebAPIException &e){
        ;
      }
    }

    // Maybe this is a SAM dataset or query?
    return new SAMQuerySource(str, stride, offset);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const Var& var,
                                       const SpillCut& spillcut,
                                       const Cut& cut,
                                       const SystShifts& shift,
                                       const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    Spectrum** ps = new Spectrum*;
    *ps = &spect;
    fHistDefs[spillcut][shift][cut][wei][var].spects.push_back(ps);

    // Remember we have a Go() pending
    spect.AddLoader(ps);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const MultiVar& var,
                                       const SpillCut& spillcut,
                                       const Cut& cut,
                                       const SystShifts& shift,
                                       const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    Spectrum** ps = new Spectrum*;
    *ps = &spect;
    fHistDefs[spillcut][shift][cut][wei][var].spects.push_back(ps);

    // Remember we have a Go() pending
    spect.AddLoader(ps);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const SpillVar& var,
                                       const SpillCut& cut,
                                       const SpillVar& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    Spectrum** ps = new Spectrum*;
    *ps = &spect;
    fSpillHistDefs[cut][wei][var].spects.push_back(ps);

    // Remember we have a Go() pending
    spect.AddLoader(ps);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const SpillMultiVar& var,
                                       const SpillCut& cut,
                                       const SpillVar& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    Spectrum** ps = new Spectrum*;
    *ps = &spect;
    fSpillHistDefs[cut][wei][var].spects.push_back(ps);

    // Remember we have a Go() pending
    spect.AddLoader(ps);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddReweightableSpectrum(ReweightableSpectrum& spect,
                                                   const Var& xvar,
                                                   const Var& yvar,
                                                   const Cut& cut,
                                                   const SystShifts& shift,
                                                   const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    ReweightableSpectrum** prw = new ReweightableSpectrum*;
    *prw = &spect;
    fHistDefs[kNoSpillCut][shift][cut][wei][xvar].rwSpects.emplace_back(prw, yvar);

    // Remember we have a Go() pending
    ReweightableSpectrumSink::AddLoader(prw);
  }

  //----------------------------------------------------------------------
  int SpectrumLoaderBase::NFiles() const
  {
    return fFileSource->NFiles();
  }

  //----------------------------------------------------------------------
  TFile* SpectrumLoaderBase::GetNextFile()
  {
    TFile* f = fFileSource->GetNextFile();
    if(!f) return 0; // out of files

    TH1* hPOT = (TH1*)f->Get("TotalPOT");
    assert(hPOT);
    fPOT += hPOT->Integral(0, -1);
    /*
    TTree* trPot = new TTree();
    if (f->GetListOfKeys()->Contains("sbnsubrun"))
      trPot = (TTree*)f->Get("sbnsubrun");
    assert(trPot);

    long n;
    // TODO should be totgoodpot?
    caf::Proxy<double> pot(0, trPot, "totpot", n, 0);

    for(n = 0; n < trPot->GetEntries(); n++){
      trPot->LoadTree(n);

      fPOT += pot;
    }
    */
    return f;
  }

  //----------------------------------------------------------------------
  void NullLoader::Go()
  {
  }

  //----------------------------------------------------------------------
  NullLoader::~NullLoader()
  {
  }

  // Apparently the existence of fSpillDefs isn't enough and I need to spell
  // this out to make sure the function bodies are generated.
  template struct SpectrumLoaderBase::IDMap<SpillCut, SpectrumLoaderBase::IDMap<SystShifts, SpectrumLoaderBase::IDMap<Cut, SpectrumLoaderBase::IDMap<Var, SpectrumLoaderBase::IDMap<SpectrumLoaderBase::VarOrMultiVar, SpectrumLoaderBase::SpectList>>>>>;

  template struct SpectrumLoaderBase::IDMap<SpillCut, SpectrumLoaderBase::IDMap<SpillVar, SpectrumLoaderBase::IDMap<SpectrumLoaderBase::SpillVarOrMultiVar, SpectrumLoaderBase::SpectList>>>;
} // namespace
