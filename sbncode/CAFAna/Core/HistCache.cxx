#include "CAFAna/Core/HistCache.h"

#include "CAFAna/Core/Utilities.h"

#include "TH2.h"

#include <iostream>

namespace ana
{
  std::multimap<int, std::unique_ptr<TH1D>> HistCache::fgMap;
  std::multimap<std::pair<int, int>, std::unique_ptr<TH2D>> HistCache::fgMap2D;

  int HistCache::fgOut = 0;
  int HistCache::fgIn = 0;

  long HistCache::fgEstMemUsage = 0;
  long HistCache::fgMemHandedOut = 0;

  //---------------------------------------------------------------------
  TH1D* HistCache::New(const std::string& title, const Binning& bins)
  {
    ++fgOut;
    fgMemHandedOut += 16*bins.NBins();

    // Look in the cache
    auto it = fgMap.find(bins.ID());
    if(it != fgMap.end()){
      TH1D* ret = it->second.release();
      fgMap.erase(it);
      ret->Reset();
      ret->SetTitle(title.c_str());

      fgEstMemUsage -= 16*bins.NBins();

      return ret;
    }

    // If not, create a new one directly
    return MakeTH1D(UniqueName().c_str(), title.c_str(), bins);
  }

  //---------------------------------------------------------------------
  TH1D* HistCache::New(const std::string& title, const TAxis* bins)
  {
    return New(title, Binning::FromTAxis(bins));
  }

  //---------------------------------------------------------------------
  TH2D* HistCache::NewTH2D(const std::string& title, const Binning& bins)
  {
    return NewTH2D(title, bins, kTrueEnergyBins);
  }

  //---------------------------------------------------------------------
  TH2D* HistCache::NewTH2D(const std::string& title, const TAxis* bins)
  {
    return NewTH2D(title, Binning::FromTAxis(bins));
  }

  //---------------------------------------------------------------------
  TH2D* HistCache::NewTH2D(const std::string& title, const Binning& xbins, const Binning& ybins)
  {
    ++fgOut;
    fgMemHandedOut += 16*xbins.NBins()*ybins.NBins();

    std::pair<int, int> IDs (xbins.ID(), ybins.ID());
    auto it = fgMap2D.find(IDs);
    if(it != fgMap2D.end()){
      TH2D* ret = it->second.release();
      fgMap2D.erase(it);
      ret->Reset();
      ret->SetTitle(title.c_str());
      fgEstMemUsage -= 16*xbins.NBins()*ybins.NBins();
      return ret;
    }

    return MakeTH2D(UniqueName().c_str(), title.c_str(), xbins, ybins);
  }

  //---------------------------------------------------------------------
  TH2D* HistCache::NewTH2D(const std::string& title, const TAxis* xbins, const TAxis* ybins)
  {
    return NewTH2D(title, Binning::FromTAxis(xbins), Binning::FromTAxis(ybins));
  }

  //---------------------------------------------------------------------
  TH1D* HistCache::Copy(const TH1D* h)
  {
    TH1D* ret = New(h->GetTitle(), h->GetXaxis());
    *ret = *h;
    return ret;
  }

  //---------------------------------------------------------------------
  TH2D* HistCache::Copy(const TH2D* h)
  {
    TH2D* ret = NewTH2D(h->GetTitle(), h->GetXaxis(), h->GetYaxis());
    *ret = *h;
    return ret;
  }

  //---------------------------------------------------------------------
  void HistCache::Delete(TH1D*& h)
  {
    if(!h) return;

    ++fgIn;
    fgMemHandedOut -= 16*h->GetNbinsX();

    fgMap.emplace(Binning::FromTAxis(h->GetXaxis()).ID(),
                  std::unique_ptr<TH1D>(h));

    fgEstMemUsage += 16*h->GetNbinsX();
    CheckMemoryUse();

    h = 0;
  }

  //---------------------------------------------------------------------
  void HistCache::Delete(TH2D*& h)
  {
    if(!h) return;

    ++fgIn;
    fgMemHandedOut -= 16*h->GetNbinsX()*h->GetNbinsY();

    fgMap2D.emplace(std::pair<int, int>(Binning::FromTAxis(h->GetXaxis()).ID(), Binning::FromTAxis(h->GetYaxis()).ID()), 
                    std::unique_ptr<TH2D>(h));

    fgEstMemUsage += 16*h->GetNbinsX()*h->GetNbinsY();
    CheckMemoryUse();

    h = 0;
  }

  //---------------------------------------------------------------------
  void HistCache::CheckMemoryUse()
  {
    if(fgEstMemUsage > 500l*1024*1024){
      std::cerr << "Warning! HistCache memory usage exceeds 500MB. "
                << "That probably means histograms are being returned "
                << "to the cache that weren't originally handed out by it. "
                << std::endl;
      PrintStats();
      std::cerr << "Now clearing cache. This could take a long time..."
                << std::endl;
      ClearCache();
      std::cerr << "Done clearing cache" << std::endl;
    }
  }

  //---------------------------------------------------------------------
  void HistCache::ClearCache()
  {
    fgMap.clear();
    fgMap2D.clear();
    fgEstMemUsage = 0;
    fgOut = 0;
    fgIn = 0;
  }

  //---------------------------------------------------------------------
  void HistCache::PrintStats()
  {
    // Count number of unique keys
    std::set<int> keys;
    for(auto& it: fgMap) keys.insert(it.first);

    std::cout << "Gave out " << fgOut << " histograms, got back "
              << fgIn << " of them (" << fgOut-fgIn << " still out, totalling "
              << fgMemHandedOut << " bytes), in "
              << keys.size() << " different shapes." << std::endl
              << "Holding " << fgMap.size()+fgMap2D.size()
              << " histograms for an estimated memory usage of "
              << fgEstMemUsage << " bytes." << std::endl;
  }
}
