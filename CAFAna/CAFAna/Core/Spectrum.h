#pragma once

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"
#include "CAFAna/Core/Utilities.h"

#include "TAttLine.h"
#include "THnSparse.h"

#include <memory>
#include <string>
#include <vector>

class TDirectory;
class TH1;
class TH2;
class TH3;
class TH1D;

/// Oscillation analysis framework, runs over CAF files outside of ART
namespace ana
{
  class Ratio;

  /// Representation of a spectrum in any variable, with associated POT
  class Spectrum
  {
  public:
    friend class SpectrumLoaderBase;
    friend class SpectrumLoader;
    friend class NullLoader;

    enum ESparse{kDense, kSparse};

    Spectrum(const std::string& label, const Binning& bins,
             SpectrumLoaderBase& loader,
             const Var& var,
             const Cut& cut,
             const SystShifts& shift = kNoShift,
             const Var& wei = kUnweighted);

    /// The only \ref MultiVar variant available
    Spectrum(const std::string& label, const Binning& bins,
             SpectrumLoaderBase& loader,
             const MultiVar& var,
             const Cut& cut,
             const SystShifts& shift = kNoShift,
             const Var& wei = kUnweighted);

    Spectrum(SpectrumLoaderBase& loader,
             const HistAxis& axis,
             const Cut& cut,
             const SystShifts& shift = kNoShift,
             const Var& wei = kUnweighted);

    Spectrum(const std::string& label, const Binning& bins, ESparse sparse = kDense);
    Spectrum(const std::string& label, double pot, double livetime, const Binning& bins);

    /// Copies \a h
    Spectrum(TH1* h,
             const std::vector<std::string>& labels,
             const std::vector<Binning>& bins,
             double pot, double livetime);

    /// Takes possession of \a h
    Spectrum(std::unique_ptr<TH1D> h,
             const std::vector<std::string>& labels,
             const std::vector<Binning>& bins,
             double pot, double livetime);

    /// 2D Spectrum of two Vars
    Spectrum(const std::string& label, SpectrumLoaderBase& loader,
             const Binning& binsx, const Var& varx,
             const Binning& binsy, const Var& vary,
             const Cut& cut,
             const SystShifts& shift = kNoShift,
             const Var& wei = kUnweighted);

    /// 2D Spectrum taking 2 HistAxis
    Spectrum(SpectrumLoaderBase& loader,
             const HistAxis& xAxis,
             const HistAxis& yAxis,
             const Cut& cut,
             const SystShifts& shift = kNoShift,
             const Var& wei = kUnweighted);

    Spectrum(const std::string& xLabel,
	     const std::string& yLabel,
	     SpectrumLoaderBase& loader,
	     const Binning& binsx, const Var& varx,
	     const Binning& binsy, const Var& vary,
	     const Cut& cut,
	     const SystShifts& shift = kNoShift,
	     const Var& wei = kUnweighted);

    /// 3D Spectrum of three Vars
    Spectrum(const std::string& label, SpectrumLoaderBase& loader,
             const Binning& binsx, const Var& varx,
             const Binning& binsy, const Var& vary,
             const Binning& binsz, const Var& varz,
             const Cut& cut,
             const SystShifts& shift = kNoShift,
             const Var& wei = kUnweighted,
	     ESparse sparse = kDense);

    Spectrum(const std::string& xLabel,
	     const std::string& yLabel,
	     const std::string& zLabel,
	     SpectrumLoaderBase& loader,
             const Binning& binsx, const Var& varx,
             const Binning& binsy, const Var& vary,
             const Binning& binsz, const Var& varz,
             const Cut& cut,
             const SystShifts& shift = kNoShift,
             const Var& wei = kUnweighted,
	     ESparse sparse = kDense);

    /// 3D Spectrum taking 3 HistAxis
    Spectrum(SpectrumLoaderBase& loader,
             const HistAxis& xAxis,
             const HistAxis& yAxis,
	     const HistAxis& zAxis,
             const Cut& cut,
             const SystShifts& shift = kNoShift,
             const Var& wei = kUnweighted,
	     ESparse sparse = kDense);


    virtual ~Spectrum();

    Spectrum(const Spectrum& rhs);
    Spectrum(Spectrum&& rhs);
    Spectrum& operator=(const Spectrum& rhs);
    Spectrum& operator=(Spectrum&& rhs);

    void Fill(double x, double w = 1);

    /// \brief Histogram made from this Spectrum, scaled to some exposure
    ///
    /// \param exposure POT or livetime (seconds)
    /// \param col Histogram color (default black)
    /// \param style Histogram line style (default solid)
    /// \param expotype How to interpret exposure (kPOT (default) or kLivetime)
    TH1D* ToTH1(double exposure,
                Color_t col = kBlack,
                Style_t style = kSolid,
                EExposureType expotype = kPOT,
                EBinType bintype = kBinContent) const;

    /// \brief Histogram made from this Spectrum, scaled to some exposure
    ///
    /// \param exposure POT or livetime (seconds)
    /// \param expotype How to interpret exposure (kPOT (default) or kLivetime)
    TH1D* ToTH1(double exposure,
                EExposureType expotype,
                EBinType bintype = kBinContent) const
    {
      return ToTH1(exposure, kBlack, kSolid, expotype, bintype);
    }

    /// Spectrum must be 2D to obtain TH2
    TH2*  ToTH2     (double exposure, EExposureType expotype = kPOT,
		     EBinType bintype = kBinContent) const;
    /// Spectrum must be 2D to obtain TH2. Normalized to X axis.
    TH2*  ToTH2NormX(double exposure, EExposureType expotype = kPOT) const;

    /// Spectrum must be 3D to obtain TH3
    TH3*  ToTH3     (double exposure, EExposureType expotype = kPOT,
		     EBinType bintype = kBinContent) const;

    /// Function decides what is the appropriate projection based on fBins, and does that
    TH1*  ToTHX     (double exposure, bool force1D = false, EExposureType expotype = kPOT) const;

    /// \brief Return total number of events scaled to \a pot
    ///
    /// \param exposure POT/livetime to scale to
    /// \param err      The statistical error on this total (optional)
    /// \param expotype What the first parameter represents
    double Integral(double exposure, double* err = 0,
		    EExposureType expotype = kPOT) const;

    /// \brief Return mean of 1D histogram
    double Mean() const;

    /// \brief Mock data is \ref FakeData with Poisson fluctuations applied
    ///
    /// Use for low-budget MDCs, or just getting a sense of the expected scale
    /// of statistical variation
    Spectrum MockData(double pot, bool makethrow=false, int seed=0) const;
    /// \brief Fake data is a MC spectrum scaled to the POT expected in the data
    ///
    /// Use for sensitivity plots and testing fit convergence
    Spectrum FakeData(double pot) const;

    double POT() const {return fPOT;}

    /// Seconds. For informational purposes only. No calculations use this.
    double Livetime() const {return fLivetime;}

    /// DO NOT USE UNLESS YOU ARE 110% CERTAIN THERE ISN'T A BETTER WAY!
    void OverridePOT(double newpot) {fPOT = newpot;}

    /// DO NOT USE UNLESS YOU ARE 110% CERTAIN THERE ISN'T A BETTER WAY!
    void OverrideLivetime(double newlive) {fLivetime = newlive;}

    void Clear();

    /// Multiply this spectrum by a constant c
    void Scale(double c);

    // Arithmetic operators are as if these are unlike samples, each a
    // contribution to one total, not seperate sources of stats for the same
    // sample.
    Spectrum& operator+=(const Spectrum& rhs);
    Spectrum operator+(const Spectrum& rhs) const;

    Spectrum& operator-=(const Spectrum& rhs);
    Spectrum operator-(const Spectrum& rhs) const;

    Spectrum& operator*=(const Ratio& rhs);
    Spectrum operator*(const Ratio& rhs) const;

    Spectrum& operator/=(const Ratio& rhs);
    Spectrum operator/(const Ratio& rhs) const;

    void SaveTo(TDirectory* dir) const;
    static std::unique_ptr<Spectrum> LoadFrom(TDirectory* dir);

    unsigned int NDimensions() const{return fLabels.size();}
    std::vector<std::string> GetLabels() const {return fLabels;}
    std::vector<Binning> GetBinnings() const {return fBins;}

  protected:
    Spectrum(const std::vector<std::string>& labels,
             const std::vector<Binning>& bins,
             ESparse sparse = kDense);

    void ConstructHistogram(ESparse sparse = kDense);

    void RemoveLoader(SpectrumLoaderBase*);
    void AddLoader(SpectrumLoaderBase*);

    /// Helper for operator+= and operator-=
    Spectrum& PlusEqualsHelper(const Spectrum& rhs, int sign);

    TH1D* fHist;
    THnSparseD* fHistSparse;
    double fPOT;
    double fLivetime;

    /// This count is maintained by SpectrumLoader, as a sanity check
    std::set<SpectrumLoaderBase*> fLoaderCount;

    std::vector<std::string> fLabels;
    std::vector<Binning> fBins;
  };

  // Commutative
  inline Spectrum operator*(const Ratio& lhs, const Spectrum& rhs){return rhs*lhs;}
  inline Spectrum operator/(const Ratio& lhs, const Spectrum& rhs){return rhs/lhs;}
}
