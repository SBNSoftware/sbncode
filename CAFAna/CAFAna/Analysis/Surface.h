#pragma once

#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Core/SystShifts.h"

#include "CAFAna/Analysis/Fit.h"

#include "Rtypes.h"
#include "TAttMarker.h"

#include <iostream>
#include <map>
#include <memory>

class TGraph;
class TH2;
class TH2F;

namespace osc{class IOscCalculatorAdjustable;}

namespace ana
{
  class IExperiment;
  class IFitVar;

  /// Log-likelihood scan across two parameters
  class Surface
  {
  public:
    friend class NumuSurface;
    friend class NueSurface;

    /// \param expt The experiment object to draw \f$ \chi^2 \f$ values from
    /// \param calc Values for oscillation parameters to be held fixed
    /// \param xvar Oscillation parameter to place on the x axis
    /// \param nbinsx Number of bins along x axis
    /// \param xmin Minimum value of x axis
    /// \param xmax Maximum value of x axis
    /// \param nbinsy Number of bins along y axis
    /// \param ymin Minimum value of y axis
    /// \param ymax Maximum value of y axis
    /// \param profVars Oscillation parameters to profile over
    /// \param profSysts Systematic parameters to profile over
    /// \param seedPts Try all combinations of these params as seeds
    /// \param systSeedPts Try all of these systematic combinations as seeds
    /// \param parallel Use all the cores on this machine? Be careful...
    Surface(const IExperiment* expt,
            osc::IOscCalculatorAdjustable* calc,
            const IFitVar* xvar, int nbinsx, double xmin, double xmax,
            const IFitVar* yvar, int nbinsy, double ymin, double ymax,
            const std::vector<const IFitVar*>& profVars = {},
            const std::vector<const ISyst*>& profSysts = {},
            const std::map<const IFitVar*, std::vector<double>>& seedPts = {},
            const std::vector<SystShifts>& systSeedPts = {},
            bool parallel = false,
            Fitter::Precision prec = Fitter::kNormal);

    /// Draw the surface itself
    void Draw() const;
    /// Draw the best fit point
    void DrawBestFit(Color_t color, Int_t marker=kFullCircle) const;
    double MinChi() const {return fMinChi;}
    double GetMinX() const {return fMinX;} 
    double GetMinY() const {return fMinY;}

    /// \param fc Surface to compare against for this significance level
    /// \param style Line style for TAttLine, solid, dotted, dashed etc
    /// \param color Line color for TAttLine
    /// \param minchi \f$\chi^2\f$ of best fit to compare against.
    ///               Default: best fit from this surface.
    void DrawContour(TH2* fc, Style_t style, Color_t color,
                     double minchi = -1);
    TH2* ToTH2(double minchi = -1) const;
    void SetTitle(const char* str);

    /// Maps of the values taken on by the profiled parameters
    std::vector<TH2*> GetProfiledHists() {return fProfHists;}
    /// Deprecated. Retained for backwards compatibility.
    std::vector<TH2*> GetMarginalizedHists() {return fProfHists;}

    /// For expert use, custom painting of contours
    std::vector<TGraph*> GetGraphs(TH2* fc, double minchi = -1);

    void SaveTo(TDirectory * dir)const;
    static std::unique_ptr< Surface > LoadFrom(TDirectory * dir);

  protected:
    Surface(){};

    void EnsureAxes() const;

    void FillSurface(const std::string& progTitle,
                     const IExperiment* expt,
                     osc::IOscCalculatorAdjustable* calc,
                     const IFitVar* xvar, const IFitVar* yvar,
                     const std::vector<const IFitVar*>& profVars,
                     const std::vector<const ISyst*>& profSysts,
                     const std::map<const IFitVar*, std::vector<double>>& seedPts,
                     const std::vector<SystShifts>& systSeedPts);

    void FillSurfacePoint(const IExperiment* expt,
                          osc::IOscCalculatorAdjustable* calc,
                          const IFitVar* xvar, double x,
                          const IFitVar* yvar, double y,
                          const std::vector<const IFitVar*>& profVars,
                          const std::vector<const ISyst*>& profSysts,
                          const std::map<const IFitVar*, std::vector<double>>& seedPts,
                          const std::vector<SystShifts>& systSeedPts);

    bool fParallel;
    Fitter::Precision fPrec;

    double fMinChi;
    double fMinX, fMinY; // Best fit point
    TH2F* fHist;
    std::vector<TH2*> fProfHists;
    std::vector<double> fSeedValues;
  };

  /// Up-value surface for 68% confidence in 2D in gaussian approximation
  TH2* Gaussian68Percent2D(const Surface& s);
  /// Up-value surface for 90% confidence in 2D in gaussian approximation
  TH2* Gaussian90Percent2D(const Surface& s);
  /// Up-value surface for 95% confidence in 2D in gaussian approximation
  TH2* Gaussian95Percent2D(const Surface& s);
  /// Up-value surface for 2 sigma confidence in 2D in gaussian approximation
  TH2* Gaussian2Sigma2D   (const Surface& s);
  /// Up-value surface for 99% confidence in 2D in gaussian approximation
  TH2* Gaussian99Percent2D(const Surface& s);
  /// Up-value surface for 3 sigma confidence in 2D in gaussian approximation
  TH2* Gaussian3Sigma2D   (const Surface& s);

  // First approximation of the correct up-values to use for ss2th13 vs delta

  /// Up-value surface for 68% confidence in 2D in gaussian approximation
  TH2* Gaussian68Percent1D(const Surface& s);
  /// Up-value surface for 90% confidence in 2D in gaussian approximation
  TH2* Gaussian90Percent1D(const Surface& s);
  /// Up-value surface for 95% confidence in 2D in gaussian approximation
  TH2* Gaussian95Percent1D(const Surface& s);
  /// Up-value surface for 2 sigma confidence in 2D in gaussian approximation
  TH2* Gaussian2Sigma1D   (const Surface& s);
  /// Up-value surface for 99% confidence in 2D in gaussian approximation
  TH2* Gaussian99Percent1D(const Surface& s);
  /// Up-value surface for 3 sigma confidence in 2D in gaussian approximation
  TH2* Gaussian3Sigma1D   (const Surface& s);
}
