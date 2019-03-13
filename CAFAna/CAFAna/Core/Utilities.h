#pragma once

#include <fenv.h>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <memory>

// these are templated types.
// can't forward-declare them here
// because compiler errors result
// when the templates are introduced
#include "TMatrixD.h"
#include "TVectorD.h"

class TArrayD;
class TDirectory;
class TH1;
class TH2;
class TH3;
class TF1;
class TH1D;
class TH2F;
class TH2D;
class TH3D;
class TVector3;

namespace ana
{
  class Binning;
  class Spectrum;
  class Ratio;

  enum EBinType
  {
    kBinContent, ///< Regular histogram
    kBinDensity  ///< Divide bin contents by bin widths
  };


  /// For use as an argument to \ref Spectrum::ToTH1
  enum EExposureType{
    kPOT,
    kLivetime
  };


  /// Return a different string each time, for creating histograms
  std::string UniqueName();

  /// \brief Prevent histograms being added to the current directory
  ///
  /// Upon going out of scope, restores the previous setting
  class DontAddDirectory
  {
  public:
    DontAddDirectory();
    ~DontAddDirectory();
  protected:
    bool fBackup;
  };

  /// \brief Alter floating-point exception flag
  ///
  /// Upon going out of scope, restores the previous setting
  class FloatingExceptionOnNaN
  {
  public:
    FloatingExceptionOnNaN(bool enable = true);
    ~FloatingExceptionOnNaN();
  protected:
    fexcept_t fBackup;
  };

  /** \brief Compute bin-to-bin covariance matrix from a collection of sets of bin contents.

      \param binSets   Collection of sets of bins from which covariances should be calculated
      \param firstBin  The first bin that should be considered (inclusive)
      \param lastBin   The last bin that should be considered (inclusive).  -1 means "last in set"

      \returns  unique_ptr to TMatrixD containing computed covariance matrix unless binSets.size() < 2,
                in which case the unique_ptr's target is nullptr.

      Note TH1D is a child class of TArrayD -- so you can pass a vector
      of TH1D* to this method.
  **/
  std::unique_ptr<TMatrixD> CalcCovMx(const std::vector<TArrayD*> & binSets, int firstBin=0, int lastBin=-1);

  class LLPerBinFracSystErr
  {
  public:
    static void SetError(double e) {fgErr = e;}
    static double GetError() {return fgErr;}
  protected:
    static double fgErr;
  };

  /** \brief The log-likelihood formula from the PDG.

      \param exp The expected spectrum
      \param obs The corresponding observed spectrum

      \returns The log-likelihood formula from the PDG
      \f[ \chi^2=2\sum_i^{\rm bins}\left(e_i-o_i+o_i\ln\left({o_i\over e_i}\right)\right) \f]

      Includes underflow bin and an option for
      overflow bin (off by default) and handles
      zero observed or expected events correctly.
  **/
  double LogLikelihood(const TH1* exp, const TH1* obs, bool useOverflow = false);

  /** \brief The log-likelihood formula for a single bin

      \param exp Expected count
      \param obs Observed count

      \returns The log-likelihood formula from the PDG
      \f[ \chi^2=2\left(e-o+o\ln\left({o\over e}\right)\right) \f]

      Handles zero observed or expected events correctly.
  **/
  double LogLikelihood(double exp, double obs);

  double LogLikelihoodDerivative(double e, double o, double dedx);

  double LogLikelihoodDerivative(const TH1D* eh, const TH1D* oh,
                                 const std::vector<double>& dedx);

  /**  \brief Chi-squared calculation using a covariance matrix.

       \param exp   Expected bin counts
       \param obs   Observed bin counts
       \param covmxinv Inverse of covariance matrix.  Must have same dimensions as exp and obs

       \returns The chi^2 calculated according to the formula from the PDG:
       \f[ \chi^2 = (\vec{o}-\vec{e})^{T} V^{-1} (\vec{o} - \vec{e}) \]

       Note that this implicitly assumes Gaussian statistics for the bin counts!
  **/
  double Chi2CovMx(const TVectorD* exp, const TVectorD* obs, const TMatrixD* covmxinv);

  /// Chi-squared calculation using covariance matrix (calls the TVectorD version internally).
  double Chi2CovMx(const TH1* exp, const TH1* obs, const TMatrixD* covmxinv);

  /// \brief Internal helper for \ref Surface and \ref FCSurface
  ///
  /// Creates a histogram having bins \em centred at the min and max
  /// coordinates
  TH2F* ExpandedHistogram(const std::string& title,
                          int nbinsx, double xmin, double xmax,
                          int nbinsy, double ymin, double ymax);

  /// \brief Invert a symmetric matrix with possibly empty rows/columns.
  ///
  /// Invert a symmetric matrix that may have empty rows/columns,
  /// which (strictly speaking) make it impossible to invert the matrix.
  /// (This often arises when computing covariance matrices for predictions
  ///  which have empty bins; the covariance is 0 for the entire row/column
  ///  in that case.)
  /// Since those rows/cols are not useful, we can sidestep the problem
  /// by removing them (and the corresponding columns)
  /// from the matrix, inverting that, then re-inserting
  /// the null rows/columns.
  std::unique_ptr<TMatrixD> SymmMxInverse(const TMatrixD& mx);

  /// Utility function to avoid need to switch on bins.IsSimple()
  TH1D* MakeTH1D(const char* name, const char* title, const Binning& bins);
  /// Utility function to avoid 4-way combinatorial explosion on the bin types
  TH2D* MakeTH2D(const char* name, const char* title,
                 const Binning& binsx,
                 const Binning& binsy);

  /// \brief For use with \ref Var2D
  ///
  /// Re-expand a histogram flattened by \ref Var2D into a 2D histogram for
  /// plotting purposes. The binning scheme must match that used in the
  /// original Var.
  TH2* ToTH2(const Spectrum& s, double exposure, ana::EExposureType expotype,
             const Binning& binsx, const Binning& binsy,
	     ana::EBinType bintype = ana::EBinType::kBinContent);

  /// Same as ToTH2, but with 3 dimensions
  TH3* ToTH3(const Spectrum& s, double exposure, ana::EExposureType expotype,
             const Binning& binsx, const Binning& binsy, const Binning& binsz,
	     ana::EBinType bintype = ana::EBinType::kBinContent);

  /// \brief For use with \ref Var2D
  ///
  /// Re-expand a flatenned histogram into a 2D histogram for
  /// plotting purposes. The binning scheme must match that used in the
  /// original Var.
  TH2* ToTH2(const Ratio& r, const Binning& binsx, const Binning& binsy);

  /// Same as ToTH2, but with 3 dimensions
  TH3* ToTH3(const Ratio& r, const Binning& binsx,
	     const Binning& binsy, const Binning& binsz);

  /// Helper for ana::ToTH2
  TH2* ToTH2Helper(std::unique_ptr<TH1> h1,
		   const Binning& binsx,
		   const Binning& binsy,
		   ana::EBinType bintype = ana::EBinType::kBinContent);

  /// Helper for ana::ToTH3
  TH3* ToTH3Helper(std::unique_ptr<TH1> h1,
		   const Binning& binsx,
		   const Binning& binsy,
		   const Binning& binsz,
		   ana::EBinType bintype = ana::EBinType::kBinContent);

  /// Find files matching a UNIX glob, plus expand environment variables
  std::vector<std::string> Wildcard(const std::string& wildcardString);

  /// This is $SRT_PRIVATE_CONTEXT if a CAFAna directory exists there,
  /// otherwise $SRT_PUBLIC_CONTEXT
  std::string FindCAFAnaDir();

  /// Read list of input files from a text file, one per line
  std::vector<std::string> LoadFileList(const std::string& listfile);

  /// \brief Extract map of metadata parameters from a CAF file
  ///
  /// \param dir The "meta" directory from the CAF file
  /// \return    A map from metadata field name to metadata value
  std::map<std::string, std::string> GetCAFMetadata(TDirectory* dir);

  /// \brief \a base += \a add
  ///
  /// \param base The original source of strings, will be altered
  /// \param add  Strings to add to \a base if missing
  /// \param mask Fields for which there was a conflict, will be altered
  void CombineMetadata(std::map<std::string, std::string>& base,
                       const std::map<std::string, std::string>& add,
                       std::set<std::string>& mask);

  /// \brief Write map of metadata parameters into a CAF file
  ///
  /// \param dir  The "meta" directory of the CAF file
  /// \param meta Map from metadata field name to metadata value
  void WriteCAFMetadata(TDirectory* dir,
                        const std::map<std::string, std::string>& meta);

  /// Is this a grid (condor) job?
  bool RunningOnGrid();

  bool AlmostEqual(double a, double b);

  std::string pnfs2xrootd(std::string loc, bool unauth = false);

  // Calling this function will return a Fourier series, fit to the input
  // histogram.  Assumes x-axis covers one period
  class FitToFourier
  {
  public:
    FitToFourier(TH1* h, double xlo, double xhi, int NOsc);
    ~FitToFourier();
    TF1* Fit() const;
    double operator()(double *x, double *par) const;
  private:

    const TH1*   fHist; // Histogram to fit
    const double fxlo;  // Lower bound
    const double fxhi;  // Upper bound - assumed to be 1 osc from the low end
    const int    fNOsc; // Highest harmonic to include

  };

  void EnsurePositiveDefinite(TH2* mat);

  /// Returns a masking histogram based on axis limits
  TH1* GetMaskHist(const Spectrum& s,
		   double xmin=0, double xmax=-1,
		   double ymin=0, double ymax=-1);
}
