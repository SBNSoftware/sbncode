#pragma once

#include "CAFAna/Core/Spectrum.h"

#include <string>

class TDirectory;
class TH2;
class TH2D;

namespace ana
{
  /// %Spectrum with the value of a second variable, allowing for reweighting
  class ReweightableSpectrum
  {
  public:
    friend class SpectrumLoaderBase;
    friend class SpectrumLoader;
    friend class NullLoader;
    friend class MRCCLoader;

    ReweightableSpectrum(SpectrumLoaderBase& loader,
                         const HistAxis& recoAxis,
                         const HistAxis& trueAxis,
                         const Cut& cut,
                         const SystShifts& shift = kNoShift,
                         const Var& wei = kUnweighted);

    ReweightableSpectrum(const Var& rwVar,
                         const std::string& xlabel, const std::string& ylabel,
                         double pot,
                         int nbinsx, double xmin, double xmax,
                         int nbinsy, double ymin, double ymax);

    ReweightableSpectrum(const Var& rwVar,
                         TH2* h,
                         const std::vector<std::string>& labels,
                         const std::vector<Binning>& bins,
                         double pot, double livetime);

    ReweightableSpectrum(const Var& rwVar,
                         std::unique_ptr<TH2D> h,
                         const std::vector<std::string>& labels,
                         const std::vector<Binning>& bins,
                         double pot, double livetime);

    virtual ~ReweightableSpectrum();

    ReweightableSpectrum(const ReweightableSpectrum& rhs);
    ReweightableSpectrum& operator=(const ReweightableSpectrum& rhs);

    /// \brief The variable that will be used to fill the y-axis
    ///
    /// By convention, return zero if the information can't be obtained, and
    /// this event will be skipped.
    const Var& ReweightVar() const {return fRWVar;}

    void Fill(double x, double y, double w = 1);

    TH2D* ToTH2(double pot) const;

    Spectrum UnWeighted() const;

    Spectrum WeightingVariable() const;

    Spectrum WeightedBy(const TH1* weights) const;

    /// Rescale bins so that \ref WeightingVariable will return \a target
    void ReweightToTrueSpectrum(const Spectrum& target);
    /// Recale bins so that \ref Unweighted will return \a target
    void ReweightToRecoSpectrum(const Spectrum& target);

    void Clear();

    /// Function to save a ReweightableSpectrum to file
    /// the fRWVar member is not written to file, so when
    /// the spectrum is loaded back from file, ReweightVar
    /// should not be accessed, but reweighting still works
    void SaveTo(TDirectory* dir) const;

    static std::unique_ptr<ReweightableSpectrum> LoadFrom(TDirectory* dir);

    unsigned int NDimensions() const{return fLabels.size();}
    std::vector<std::string> GetLabels() const {return fLabels;}
    std::vector<Binning> GetBinnings() const {return fBins;}

  protected:
    // Derived classes can be trusted take care of their own construction
    ReweightableSpectrum(const std::vector<std::string>& labels,
                         const std::vector<Binning>& bins,
                         const Var& rwVar)
      : fRWVar(rwVar),
        fHist(0), fPOT(0), fLivetime(0),
        fLabels(labels), fBins(bins)
    {
    }

    ReweightableSpectrum(const std::string& label,
                         const Binning& bins,
                         const Var& rwVar)
      : fRWVar(rwVar),
        fHist(0), fPOT(0), fLivetime(0),
        fLabels(1, label), fBins(1, bins)
    {
    }

    /// Constructor needed by LoadFrom. Since there's no good
    /// way to store a Var, ReweightVar will return nonsense
    /// for ReweightableSpectrum that are loaded from a file
    ReweightableSpectrum(TH2* h,
                         const std::vector<std::string>& labels,
                         const std::vector<Binning>& bins,
                         double pot, double livetime)
      : ReweightableSpectrum(kUnweighted, h, labels, bins, pot, livetime)
    {
    }

    void RemoveLoader(SpectrumLoaderBase*);
    void AddLoader(SpectrumLoaderBase*);

    Var fRWVar; ///< What goes on the y axis?

    TH2D* fHist;
    double fPOT;
    double fLivetime;

    std::vector<std::string> fLabels;
    std::vector<Binning> fBins;

    std::string fTrueLabel;

    /// This count is maintained by SpectrumLoader, as a sanity check
    std::set<SpectrumLoaderBase*> fLoaderCount;
  };
}
