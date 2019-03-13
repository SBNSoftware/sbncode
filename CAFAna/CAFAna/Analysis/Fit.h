#pragma once

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Core/SystShifts.h"

#include "Math/Minimizer.h"

#include <memory>

class TGraph;

namespace osc{class IOscCalculatorAdjustable;}

namespace ana
{
  class IExperiment;
  class IFitVar;

  /// \brief Figure-of-merit with no systematics, for binned data
  ///
  /// \param obs The observed data
  /// \param unosc A spectrum representing the null hypothesis
  /// \param pot POT to scale to. Leave at default to adopt POT from \a obs
  ///
  /// \returns Sum in quadrature over the bins
  /// \f[ \sqrt{\sum_i^{\rm bins}\left({s_i\over\sqrt{s_i+b_i}}\right)^2} \f]
  double SimpleFOM(const Spectrum& obs, const Spectrum& unosc, double pot = 0);

  namespace{SystShifts junkShifts;}

  /// Perform MINUIT fits in one or two dimensions
  class Fitter: public ROOT::Math::IGradientFunctionMultiDim
  {
  public:
    enum Verbosity{kQuiet, kVerbose};

    enum Precision{
      // You must select one of these. The first three codes match the settings
      // used by migrad. The fourth is a custom minimizer.
      kFast = 0,
      kNormal = 1,
      kCareful = 2,
      kGradDesc = 3,
      // Allow bitmask operations to extract these first four options
      kAlgoMask = 3,
      // You may optionally specify this (eg kNormal | kIncludeSimplex) to
      // improve the chances of escaping from invalid minima
      kIncludeSimplex = 4,
      // You may optionally specify these to improve the final error estimates
      kIncludeHesse = 8,
      kIncludeMinos = 16
    };
    void SetPrecision(Precision prec);

    Fitter(const IExperiment* expt,
           std::vector<const IFitVar*> vars,
           std::vector<const ISyst*> systs = {},
           Precision prec = kNormal);

    /// \param[out] seed Seed parameter and output best-fit point
    /// \param[out] bestSysts Best systematics result returned here
    /// \param      seedPts Set each var to each of the values. Try all
    ///                     combinations. Beware of combinatorical explosion...
    /// \param      systSeedPts If non-empty, try fit starting at each of these
    /// \param      verb If quiet, no printout
    /// \return     -2x the log-likelihood of the best-fit point
    double Fit(osc::IOscCalculatorAdjustable* seed,
               SystShifts& bestSysts = junkShifts,
               const std::map<const IFitVar*, std::vector<double>>& seedPts = {},
               const std::vector<SystShifts>& systSeedPts = {},
               Verbosity verb = kVerbose) const;

    /// Variant with no seedPts
    double Fit(osc::IOscCalculatorAdjustable* seed,
               SystShifts& systSeed,
               Verbosity verb) const
    {
      return Fit(seed, systSeed, {}, std::vector<SystShifts>(1, systSeed), verb);
    }

    /// Variant with no seedPts and no systematics result returned
    double Fit(osc::IOscCalculatorAdjustable* seed,
               Verbosity verb) const
    {
      return Fit(seed, junkShifts, {}, {}, verb);
    }

    /// Variant with no oscillations - useful for ND fits
    double Fit(SystShifts& systSeed, Verbosity verb = kVerbose) const
    {
      return Fit(0, systSeed, {}, std::vector<SystShifts>(1, systSeed), verb);
    }

    /// Return the fit covariance
    TMatrixDSym* GetCovariance(){ return this->fCovar;}

    /// covariance matrix status (0 = not valid, 1 approximate, 2, full but made pos def, 3 accurate and not pos def)
    int GetCovarianceStatus() {return this->fCovarStatus;}

    /// Return the fit names
    std::vector<std::string> GetParamNames(){ return this->fParamNames;}

    /// Return the prefit values
    std::vector<double> GetPreFitValues(){ return this->fPreFitValues;}

    /// Return the prefit errors
    std::vector<double> GetPreFitErrors(){ return this->fPreFitErrors;}

    /// Return the postfit values
    std::vector<double> GetPostFitValues(){ return this->fPostFitValues;}

    /// Return the postfit errors
    std::vector<double> GetPostFitErrors(){ return this->fPostFitErrors;}

    /// Return the minos errors
    std::vector<std::pair<double,double>> GetMinosErrors(){ return this->fMinosErrors;}

    /// Return number of function calls
    int GetNFCN(){return this->fNEval;}

    /// Return edm form the fit
    double GetEDM(){return this->fEdm;}

    /// Say whether the fit was good
    bool GetIsValid(){return this->fIsValid;}

    SystShifts GetSystShifts() const {return fShifts;}

    /// Evaluate the log-likelihood, as required by MINUT interface
    virtual double DoEval(const double* pars) const override;

    // Part of the fitter interface
    virtual unsigned int NDim() const override {return fVars.size()+fSysts.size();}

    virtual void Gradient(const double* x, double* grad) const override;

    virtual double DoDerivative(const double* x, unsigned int icoord) const override
    {
      std::cout << "Fitter::DoDerivative() not implemented" << std::endl;
      abort();
    }

    Fitter* Clone() const override
    {
      std::cout << "Fitter::Clone() not implemented" << std::endl;
      abort();
    }


    // TODO unused
    bool CheckGradient() const{return (fPrec & Fitter::kAlgoMask) != kFast;}
  protected:
    struct SeedPt
    {
      std::map<const IFitVar*, double> fitvars;
      SystShifts shift;
    };
    std::vector<SeedPt> ExpandSeeds(const std::map<const IFitVar*,
                                                   std::vector<double>>& seedPts,
                                    std::vector<SystShifts> systSeedPts) const;

    /// Helper for \ref FitHelper
    std::unique_ptr<ROOT::Math::Minimizer>
    FitHelperSeeded(osc::IOscCalculatorAdjustable* seed,
                    SystShifts& systSeed,
                    Verbosity verb) const;

    /// Helper for \ref Fit
    double FitHelper(osc::IOscCalculatorAdjustable* seed,
                     SystShifts& bestSysts,
                     const std::map<const IFitVar*, std::vector<double>>& seedPts,
                     std::vector<SystShifts> systSeedPts,
                     Verbosity verb) const;

    /// Updates mutable fCalc and fShifts
    void DecodePars(const double* pars) const;

    /// Intended to be called only once (from constructor) to initialize
    /// fSupportsDerivatives
    bool SupportsDerivatives() const;

    const IExperiment* fExpt;
    std::vector<const IFitVar*> fVars;
    std::vector<const ISyst*> fSysts;
    Precision fPrec = kNormal;
    mutable osc::IOscCalculatorAdjustable* fCalc;
    mutable SystShifts fShifts;

    bool fSupportsDerivatives;

    mutable int fNEval = 0;
    mutable int fNEvalGrad = 0;
    mutable int fNEvalFiniteDiff = 0;

    // Some information for post-fit evaluation if necessary
    mutable double fEdm = -1;
    mutable bool fIsValid = false;
    mutable TMatrixDSym* fCovar;
    mutable bool fCovarStatus;
    mutable std::vector<std::string> fParamNames;
    mutable std::vector<double> fPreFitValues;
    mutable std::vector<double> fPreFitErrors;
    mutable std::vector<double> fPostFitValues;
    mutable std::vector<double> fPostFitErrors;
    mutable std::vector<std::pair<double, double> > fMinosErrors;
    mutable std::vector<std::pair<double, double> > fTempMinosErrors; // Bit of a hack
  };

  // Modern C++ thinks that enum | enum == int. Make things work like we expect
  // for this bitmask.
  inline Fitter::Precision operator|(Fitter::Precision a, Fitter::Precision b)
  {
    return Fitter::Precision(int(a) | int(b));
  }

  // Default values for Profile()
  static std::map<const IFitVar*, TGraph*> empty_vars_map;
  static std::map<const ISyst*,   TGraph*> empty_syst_map;

  /// \brief \f$\chi^2\f$ scan in one variable, profiling over all others
  ///
  /// \param expt   The experiment to retrieve chisq values from
  /// \param calc   Initial values of all oscillation parameters
  /// \param v      Scan over this variable
  /// \param nbinsx Binning
  /// \param minx   Binning
  /// \param maxx   Binning
  /// \param minchi Set non-default to force a chisq value to evaluate delta
  ///               chisqs against. Useful for comparing two profiles. If set
  ///               to zero it will not zero-adjust and the axis will be
  ///               labelled without "delta"
  /// \param profVars  Profile over these variables
  /// \param profSysts Profile over these systematics
  /// \param seedPts   Set each var to each of the values. Try all
  ///                  combinations. Beware of combinatorical explosion...
  /// \param      systSeedPts If non-empty, try fit starting at each of these
  /// \param[out] profVarsMap Pass empty map. Returns best values of each var.
  /// \param[out] systsMap    Pass empty map. Returns best values of each syst.
  ///
  /// \return The best fit delta chisq as a function of \a a
  TH1* Profile(const IExperiment* expt,
	       osc::IOscCalculatorAdjustable* calc,
               const IFitVar* v,
	       int nbinsx, double minx, double maxx,
	       double minchi = -1,
               const std::vector<const IFitVar*>& profVars = {},
               const std::vector<const ISyst*>& profSysts = {},
               const std::map<const IFitVar*, std::vector<double>>& seedPts = {},
               const std::vector<SystShifts>& systsSeedPts = {},
               std::map<const IFitVar*, TGraph*>& profVarsMap = empty_vars_map,
               std::map<const ISyst*, TGraph*>& systsMap = empty_syst_map);

  /// Forward to \ref Profile but sqrt the result for a crude significance
  TH1* SqrtProfile(const IExperiment* expt,
                   osc::IOscCalculatorAdjustable* calc,
                   const IFitVar* v,
                   int nbinsx, double minx, double maxx,
                   double minchi = -1,
                   std::vector<const IFitVar*> profVars = {},
                   std::vector<const ISyst*> profSysts = {},
                   const std::map<const IFitVar*, std::vector<double>>& seedPts = {},
                   const std::vector<SystShifts>& systsSeedPts = {},
                   std::map<const IFitVar*, TGraph*>& profVarsMap = empty_vars_map,
                   std::map<const ISyst*, TGraph*>& systsMap = empty_syst_map);

  /// \f$\chi^2\f$ scan in one variable, holding all others constant
  TH1* Slice(const IExperiment* expt,
             osc::IOscCalculatorAdjustable* calc, const IFitVar* v,
             int nbinsx, double minx, double maxx,
             double minchi = -1);

  /// Forward to \ref Slice but sqrt the result for a crude significance
  TH1* SqrtSlice(const IExperiment* expt,
                 osc::IOscCalculatorAdjustable* calc, const IFitVar* v,
                 int nbinsx, double minx, double maxx, double minchi = -1);

  /// \brief Find the minimum in one variable as a function of another
  ///
  /// \param transpose plot \a scanVar on the y axis
  TGraph* FindValley(const IExperiment* expt,
		     osc::IOscCalculatorAdjustable* calc,
		     const IFitVar& scanVar,
		     const IFitVar& fitVar,
		     int nbinsx, double xmin, double xmax,
		     const std::vector<const IFitVar*>& profVars = {},
		     const std::vector<const ISyst*>& profSysts = {},
                     const std::map<const IFitVar*, std::vector<double>>& seedPts = {},
                     const std::vector<SystShifts>& systsSeedPts = {},
		     bool transpose = false);

  /// \brief Intended for use on the output of \ref Profile
  ///
  /// Returns a list of all the x-coordinates at which the curve described by
  /// \a h crosses \a critVal. eg using critVal=1 will find the 1sigma lower
  /// and upper bounds.
  std::vector<double> FindCurveCrossings(TH1* h, double critVal);
}
