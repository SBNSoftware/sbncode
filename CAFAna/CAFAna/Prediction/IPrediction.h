#pragma once

#include "CAFAna/Core/Spectrum.h"

#include "CAFAna/Core/OscillatableSpectrum.h"

namespace osc{class IOscCalculator;}

class TDirectory;

namespace ana
{
  /// Enumeration of neutrino transition modes
  namespace Flavors{
    enum Flavors_t{
      kNuEToNuE    = 1<<0, ///< \f$\nu_e\to\nu_e\f$ ('beam \f$\nu_e \f$')
      kNuEToNuMu   = 1<<1, ///< \f$\nu_e\to\nu_\mu\f$ ('\f$\nu_\mu\f$ appearance')
      kNuEToNuTau  = 1<<2, ///< \f$\nu_e\to\nu_\tau\f$
      kNuMuToNuE   = 1<<3, ///< \f$\nu_\mu\to\nu_e\f$ ('\f$\nu_e\f$ appearance')
      kNuMuToNuMu  = 1<<4, ///< \f$\nu_\mu\to\nu_\mu\f$ ('\f$\nu_\mu\f$ survival')
      kNuMuToNuTau = 1<<5, ///< \f$\nu_\mu\to\nu_\tau\f$

      kAllNuE   = kNuEToNuE   | kNuMuToNuE,   ///< All \f$\nu_e\f$
      kAllNuMu  = kNuEToNuMu  | kNuMuToNuMu,  ///< All \f$\nu_\mu\f$
      kAllNuTau = kNuEToNuTau | kNuMuToNuTau, ///< All \f$\nu_\tau\f$

      kAll = kAllNuE | kAllNuMu | kAllNuTau   ///< All neutrinos, any flavor
    };

    inline Flavors_t operator|(Flavors_t a, Flavors_t b)
    {
      // The default definition returns an int. We don't want that
      return Flavors_t(int(a) | int(b));
    }
  }

  /// Enumeration for interaction currents (CC/NC)
  namespace Current{
    enum Current_t{
      kCC = 1<<0, ///< Charged-current interactions
      kNC = 1<<1, ///< Neutral-current interactions

      kBoth = kCC | kNC ///< Interactions of both types
    };
  }

  /// Enumeration for neutrino sign (neutrino/antineutrino)
  namespace Sign{
    enum Sign_t{
      kNu     = 1<<0, ///< Neutrinos-only
      kAntiNu = 1<<1, ///< Antineutrinos-only

      kBoth = kNu | kAntiNu ///< Both neutrinos and antineutrinos
    };
  }

  class SystShifts;

  /// Standard interface to all prediction techniques
  class IPrediction
  {
  public:
    virtual ~IPrediction(){}
    virtual Spectrum PredictUnoscillated() const;
    virtual Spectrum Predict(osc::IOscCalculator* calc) const = 0;
    virtual Spectrum PredictSyst(osc::IOscCalculator* calc,
                                 const SystShifts& syst) const;

    virtual Spectrum PredictComponent(osc::IOscCalculator* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const = 0;
    virtual Spectrum PredictComponentSyst(osc::IOscCalculator* calc,
                                          const SystShifts& syst,
                                          Flavors::Flavors_t flav,
                                          Current::Current_t curr,
                                          Sign::Sign_t sign) const;

    virtual void Derivative(osc::IOscCalculator* calc,
                            const SystShifts& shift,
                            double pot,
                            std::unordered_map<const ISyst*, std::vector<double>>& dchi) const
    {
      // Implementing this function is optional. If you don't implement it,
      // this default implementation will be used, which signals to callers
      // that your Prediction doesn't implement this feature.
      dchi.clear();
    }

    virtual OscillatableSpectrum ComponentCC(int from, int to) const
    {std::cout << "OscillatableSpectrum::ComponentCC() unimplemented" << std::endl; abort();}
    virtual Spectrum ComponentNC() const
    {std::cout << "OscillatableSpectrum::ComponentNC() unimplemented" << std::endl; abort();}

    virtual void SaveTo(TDirectory* dir) const;
  };
}
