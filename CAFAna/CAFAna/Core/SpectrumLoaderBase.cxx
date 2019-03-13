#include "CAFAna/Core/SpectrumLoaderBase.h"

#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/ReweightableSpectrum.h"
#include "CAFAna/Core/SAMQuerySource.h"
#include "CAFAna/Core/SAMProjectSource.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Utilities.h"
#include "CAFAna/Core/WildcardSource.h"

#include "StandardRecord/StandardRecord.h"

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
  // Apparently the existence of fSpillDefs isn't enough and I need to spell
  // this out to make sure the function bodies are generated.
  template class SpectrumLoaderBase::IDMap<SystShifts, SpectrumLoaderBase::IDMap<Cut, SpectrumLoaderBase::IDMap<Var, SpectrumLoaderBase::IDMap<SpectrumLoaderBase::VarOrMultiVar, SpectrumLoaderBase::SpectList>>>>;

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::Erase(Spectrum* s)
  {
    auto it = std::find(spects.begin(), spects.end(), s);
    if(it != spects.end()) spects.erase(it);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::Erase(ReweightableSpectrum* rs)
  {
    auto it = std::find(rwSpects.begin(), rwSpects.end(), rs);
    if(it != rwSpects.end()) rwSpects.erase(it);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::RemoveLoader(SpectrumLoaderBase* l)
  {
    for(Spectrum* s: spects) s->RemoveLoader(l);
    for(ReweightableSpectrum* rs: rwSpects) rs->RemoveLoader(l);
  }

  //----------------------------------------------------------------------
  size_t SpectrumLoaderBase::SpectList::TotalSize() const
  {
    return spects.size() + rwSpects.size();
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::GetSpectra(std::vector<Spectrum*>& ss)
  {
    ss.insert(ss.end(), spects.begin(), spects.end());
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::
  GetReweightableSpectra(std::vector<ReweightableSpectrum*>& ss)
  {
    ss.insert(ss.end(), rwSpects.begin(), rwSpects.end());
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
  template<class T, class U> template<class V> void SpectrumLoaderBase::IDMap<T, U>::Erase(const V& v)
  {
    for(auto& it: fElems) it.second.Erase(v);
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

  //----------------------------------------------------------------------
  template<class T, class U> void SpectrumLoaderBase::IDMap<T, U>::
  GetSpectra(std::vector<Spectrum*>& ss)
  {
    for(auto& it: fElems) it.second.GetSpectra(ss);
  }

  //----------------------------------------------------------------------
  template<class T, class U> void SpectrumLoaderBase::IDMap<T, U>::
  GetReweightableSpectra(std::vector<ReweightableSpectrum*>& ss)
  {
    for(auto& it: fElems) it.second.GetReweightableSpectra(ss);
  }

  // Start of SpectrumLoaderBase proper

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(DataSource src)
    : fSource(src), fGone(false), fPOT(0)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(const std::string& wildcard,
                                         DataSource src)
    : SpectrumLoaderBase(src)
  {
    fWildcard = wildcard;
    fFileSource = std::unique_ptr<IFileSource>(WildcardOrSAMQuery(wildcard));
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(const std::vector<std::string>& fnames,
                                         DataSource src)
    : SpectrumLoaderBase(src)
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
      else {
	std::cout << "Warning: " << str << " does not exist!" << std::endl;
	abort();
      }
      delete ret;
    }

    // Maybe this the name of a SAM project?
    ifdh_ns::ifdh i;
    const std::string info = i.dumpProject(i.findProject(str, "nova"));
    // findProject always gives back an address just by gluing bits
    // together. But dumpProject will give an empty result for a nonexistent
    // project.
    if(!info.empty()){
      if(stride > 1)
        std::cout << "Warning: --stride has no effect on SAM projects"
                  << std::endl;

      return new SAMProjectSource(str);
    }

    // Maybe this is a SAM dataset or query?
    return new SAMQuerySource(str, stride, offset);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const Var& var,
                                       const Cut& cut,
                                       const SystShifts& shift,
                                       const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fHistDefs[shift][cut][wei][var].spects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const MultiVar& var,
                                       const Cut& cut,
                                       const SystShifts& shift,
                                       const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fHistDefs[shift][cut][wei][var].spects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::RemoveSpectrum(Spectrum* spect)
  {
    fHistDefs.Erase(spect);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddReweightableSpectrum(ReweightableSpectrum& spect,
                                                   const Var& var,
                                                   const Cut& cut,
                                                   const SystShifts& shift,
                                                   const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fHistDefs[shift][cut][wei][var].rwSpects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::
  RemoveReweightableSpectrum(ReweightableSpectrum* spect)
  {
    fHistDefs.Erase(spect);
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

    TTree* trPot = new TTree();
    if (f->GetListOfKeys()->Contains("sbnsubrun"))
      trPot = (TTree*)f->Get("sbnsubrun");
    assert(trPot);

    double pot;
    trPot->SetBranchAddress("totpot", &pot); // should this be totgoodpot?

    for(int n = 0; n < trPot->GetEntries(); n++){
      trPot->GetEntry(n);

      fPOT += pot;
    }

    // This won't work for old files, where we would need to set the POT
    // to an arbitrary number (e.g. fPOT = 1.e20)

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

} // namespace

