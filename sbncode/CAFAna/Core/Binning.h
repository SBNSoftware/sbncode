#pragma once

#include <cassert>
#include <map>
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>

class TAxis;
class TDirectory;

namespace ana
{
  /// \brief Represent the binning of a Spectrum's x-axis
  ///
  /// May be "Simple" (equally spaced) or "Custom" (arbitrary binning)
  class Binning
  {
  public:
    static Binning Simple(int n, double lo, double hi,
                          const std::vector<std::string>& labels = {});
    static Binning LogUniform(int n, double lo, double hi);
    static Binning Custom(const std::vector<double>& edges);
    static Binning FromTAxis(const TAxis* ax);

    int NBins() const {return fNBins;}
    double Min() const {return fMin;}
    double Max() const {return fMax;}
    int FindBin(float x) const;
    bool IsSimple() const {return fIsSimple;}
    const std::vector<double>& Edges() const
    {
      return fEdges;
    }

    const std::vector<std::string>& Labels() const {return fLabels;}

    void SaveTo(TDirectory* dir) const;
    static std::unique_ptr<Binning> LoadFrom(TDirectory* dir);

    int ID() const {return fID;}
    static int MaxID() {return fgNextID-1;}

    bool operator==(const Binning& rhs) const;
    bool operator<(const Binning& rhs) const;
  protected:
    Binning();

    static Binning SimpleHelper(int n, double lo, double hi,
                                const std::vector<std::string>& labels = {});

    static Binning CustomHelper(const std::vector<double>& edges);

    std::vector<double> fEdges;
    std::vector<std::string> fLabels;
    int fNBins;
    double fMin, fMax;
    bool fIsSimple;

    int fID;
    /// The next ID that hasn't yet been assigned
    static int fgNextID;

    static std::map<Binning, int>& IDMap();
  };

  /// Default true-energy bin edges
  Binning TrueEnergyBins();
  /// Default true-energy bin edges
  const Binning kTrueEnergyBins = TrueEnergyBins();
}
