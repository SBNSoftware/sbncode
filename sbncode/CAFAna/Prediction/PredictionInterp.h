#pragma once

#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Prediction/PredictionGenerator.h"

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/SystShifts.h"

#include <map>
#include <memory>
#include <unordered_map>

#include "TMD5.h"

class TH1;

namespace ana
{
  class Loaders;

  /// Implements systematic errors by interpolation between shifted templates
  class PredictionInterp: public IPrediction
  {
  public:
    enum EMode_t{
      kCombineSigns, kSplitBySign
    };

    /// \param systs What systematics we should be capable of interpolating
    /// \param osc The oscillation point to expand around
    /// \param predGen Construct an IPrediction from the following
    ///                information.
    /// \param loaders The loaders to be passed on to the underlying prediction
    /// \param shiftMC Underlying shift. Use with care. Mostly for
    ///                PredictionNumuFAHadE. Should not contain any of of
    ///                \a systs
    /// \param mode    In FHC the wrong-sign has bad stats and probably the
    ///                fits can't be split out reasonably. For RHC it's
    ///                important not to conflate them.
    PredictionInterp(std::vector<const ISyst*> systs,
                     osc::IOscCalculator* osc,
                     const IPredictionGenerator& predGen,
                     Loaders& loaders,
                     const SystShifts& shiftMC = kNoShift,
                     EMode_t mode = kCombineSigns);

    virtual ~PredictionInterp();



    virtual Spectrum Predict(osc::IOscCalculator* calc) const override;
    virtual Spectrum PredictSyst(osc::IOscCalculator* calc,
                                 const SystShifts& shift) const override;

    virtual Spectrum PredictComponent(osc::IOscCalculator* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override;
    virtual Spectrum PredictComponentSyst(osc::IOscCalculator* calc,
                                          const SystShifts& shift,
                                          Flavors::Flavors_t flav,
                                          Current::Current_t curr,
                                          Sign::Sign_t sign) const override;

    virtual void Derivative(osc::IOscCalculator* calc,
                            const SystShifts& shift,
                            double pot,
                            std::unordered_map<const ISyst*, std::vector<double>>& dp) const override;

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<PredictionInterp> LoadFrom(TDirectory* dir);

    /// After calling this DebugPlots won't work fully and SaveTo won't work at
    /// all.
    void MinimizeMemory();

    void DebugPlot(const ISyst* syst,
                   osc::IOscCalculator* calc,
                   Flavors::Flavors_t flav = Flavors::kAll,
                   Current::Current_t curr = Current::kBoth,
                   Sign::Sign_t sign = Sign::kBoth) const;

    // If \a savePattern is not empty, print each pad. Must contain a "%s" to
    // contain the name of the systematic.
    void DebugPlots(osc::IOscCalculator* calc,
		    const std::string& savePattern = "",
		    Flavors::Flavors_t flav = Flavors::kAll,
		    Current::Current_t curr = Current::kBoth,
		    Sign::Sign_t sign = Sign::kBoth) const;

    void SetOscSeed(osc::IOscCalculator* oscSeed);

    void DebugPlotColz(const ISyst* syst,
                       osc::IOscCalculator* calc,
                       Flavors::Flavors_t flav = Flavors::kAll,
                       Current::Current_t curr = Current::kBoth,
                       Sign::Sign_t sign = Sign::kBoth) const;

    void DebugPlotsColz(osc::IOscCalculator* calc,
                        const std::string& savePattern = "",
                        Flavors::Flavors_t flav = Flavors::kAll,
                        Current::Current_t curr = Current::kBoth,
                        Sign::Sign_t sign = Sign::kBoth) const;

    bool SplitBySign() const {return fSplitBySign;}
    enum CoeffsType{
      kNueApp, kNueSurv, kNumuSurv, kNC,
      kOther, ///< Taus, numu appearance
      kNCoeffTypes
    };

    PredictionInterp() : fBinning(0, {}, {}, 0, 0) {}

    static void LoadFromBody(TDirectory* dir, PredictionInterp* ret,
			     std::vector<const ISyst*> veto = {});

    struct Coeffs{
      Coeffs(double _a, double _b, double _c, double _d)
        : a(_a), b(_b), c(_c), d(_d) {}
      double a, b, c, d;
    };

    /// Find coefficients describing this set of shifts
    std::vector<std::vector<Coeffs>>
    FitRatios(const std::vector<double>& shifts,
              const std::vector<std::unique_ptr<TH1>>& ratios) const;

    /// Find coefficients describing the ratios from this component
    std::vector<std::vector<Coeffs>>
    FitComponent(const std::vector<double>& shifts,
                 const std::vector<IPrediction*>& preds,
                 Flavors::Flavors_t flav,
                 Current::Current_t curr,
                 Sign::Sign_t sign) const;

    Spectrum ShiftSpectrum(const Spectrum& s,
                           CoeffsType type,
                           bool nubar, // try to use fitsNubar if it exists
                           const SystShifts& shift) const;

    /// Helper for PredictComponentSyst
    Spectrum ShiftedComponent(osc::IOscCalculator* calc,
                              const TMD5* hash,
                              const SystShifts& shift,
                              Flavors::Flavors_t flav,
                              Current::Current_t curr,
                              Sign::Sign_t sign,
                              CoeffsType type) const;

    std::unique_ptr<IPrediction> fPredNom; ///< The nominal prediction

    struct ShiftedPreds
    {
      double Stride() const {return shifts.size() > 1 ? shifts[1]-shifts[0] : 1;}

      std::string systName; ///< What systematic we're interpolating
      std::vector<double> shifts; ///< Shift values sampled
      std::vector<IPrediction*> preds;

      int nCoeffs; // Faster than calling size()

      /// Indices: [type][histogram bin][shift bin]
      std::vector<std::vector<std::vector<Coeffs>>> fits;
      /// Will be filled if signs are separated, otherwise not
      std::vector<std::vector<std::vector<Coeffs>>> fitsNubar;
    };

  protected:
    mutable std::unordered_map<const ISyst*, ShiftedPreds> fPreds;

    /// The oscillation values we assume when evaluating the coefficients
    osc::IOscCalculator* fOscOrigin;

    mutable Spectrum fBinning; ///< Dummy spectrum to provide binning

    struct Key_t
    {
      Flavors::Flavors_t flav;
      Current::Current_t curr;
      Sign::Sign_t sign;
      bool operator<(const Key_t& rhs) const
      {
        return (std::make_tuple(flav, curr, sign) <
                std::make_tuple(rhs.flav, rhs.curr, rhs.sign));
      }
    };
    struct Val_t
    {
      TMD5 hash;
      Spectrum nom;
    };
    mutable std::map<Key_t, Val_t> fNomCache;

    bool fSplitBySign;

    void InitFits() const;

    void InitFitsHelper(ShiftedPreds& sp,
                        std::vector<std::vector<std::vector<Coeffs>>>& fits,
                        Sign::Sign_t sign) const;

     /// Helper for \ref Derivative
    void ComponentDerivative(osc::IOscCalculator* calc,
                             Flavors::Flavors_t flav,
                             Current::Current_t curr,
                             Sign::Sign_t sign,
                             CoeffsType type,
                             const SystShifts& shift,
                             double pot,
                             std::unordered_map<const ISyst*, std::vector<double>>& dp) const;
  };

}
