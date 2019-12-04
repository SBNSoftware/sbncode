#pragma once

#include <tuple>
#include <map>
#include <string>

#include "CAFAna/Core/Binning.h"

class TAxis;
class TH1D;
class TH2D;

namespace ana
{
  /// \brief Helper for \ref Spectrum
  ///
  /// ROOT's handling of allocations, and especially deletions, can be very
  /// slow. It keeps everything in a big map that it then has to lookups
  /// in. This class provides a simple cache of histograms, recycling an old
  /// histogram of the same binning instead of creating a new one.
  ///
  /// Allocate new histograms with \ref New, and return them to the cache with
  /// \ref Delete.
  class HistCache
  {
  public:
    static TH1D* New(const std::string& title, const Binning& bins);
    static TH1D* New(const std::string& title, const TAxis* bins);
    static TH1D* Copy(const TH1D* h);

    /// WARNING - the first two TH2D functions are kept for backward compatibility,
    /// they take kTrueEnergyBins on the y-axis.
    static TH2D* NewTH2D(const std::string& title, const Binning& bins);
    static TH2D* NewTH2D(const std::string& title, const TAxis* bins);
    /// The following two TH2D functions can have flexible y-axis. 
    static TH2D* NewTH2D(const std::string& title, const Binning& xbins, const Binning& ybins);
    static TH2D* NewTH2D(const std::string& title, const TAxis* xbins, const TAxis* ybins);
    static TH2D* Copy(const TH2D* h);

    static void Delete(TH1D*& h);
    static void Delete(TH2D*& h);
    static void PrintStats();
    static void ClearCache();
  protected:
    static void CheckMemoryUse();

    // Key to the maps is Binning::ID()
    static std::multimap<int, std::unique_ptr<TH1D>> fgMap;
    static std::multimap<std::pair<int, int>, std::unique_ptr<TH2D>> fgMap2D;

    static int fgOut, fgIn;

    static long fgEstMemUsage;
    static long fgMemHandedOut;
  };
}
