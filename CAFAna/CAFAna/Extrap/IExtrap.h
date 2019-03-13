#pragma once

#include "CAFAna/Core/OscillatableSpectrum.h"

namespace ana
{
  /// Interface to extrapolation procedures
  class IExtrap
  {
  public:
    virtual ~IExtrap() {};

    /// Charged current electron neutrino survival (\f$\nu_e\to\nu_e\f$)
    virtual OscillatableSpectrum NueSurvComponent() = 0;
    /// Charged current electron antineutrino survival (\f$\bar\nu_e\to\bar\nu_e\f$)
    virtual OscillatableSpectrum AntiNueSurvComponent() = 0;

    /// Charged current muon neutrino survival (\f$\nu_\mu\to\nu_\mu\f$)
    virtual OscillatableSpectrum NumuSurvComponent() = 0;
    /// Charged current muon antineutrino survival (\f$\bar\nu_\mu\to\bar\nu_\mu\f$)
    virtual OscillatableSpectrum AntiNumuSurvComponent() = 0;

    /// Charged current electron neutrino appearance (\f$\nu_\mu\to\nu_e\f$)
    virtual OscillatableSpectrum NueAppComponent() = 0;
    /// Charged current electron antineutrino appearance (\f$\bar\nu_\mu\to\bar\nu_e\f$)
    virtual OscillatableSpectrum AntiNueAppComponent() = 0;

    /// Charged current muon neutrino appearance (\f$\nu_e\to\nu_\mu\f$)
    virtual OscillatableSpectrum NumuAppComponent() = 0;
    /// Charged current muon antineutrino appearance (\f$\bar\nu_e\to\bar\nu_\mu\f$)
    virtual OscillatableSpectrum AntiNumuAppComponent() = 0;

    /// Charged current tau neutrino appearance from electron neutrino (\f$\nu_e\to\nu_\tau\f$)
    virtual OscillatableSpectrum TauFromEComponent() = 0;
    /// Charged current tau antineutrino appearance from electron antineutrino (\f$\bar\nu_e\to\bar\nu_\tau\f$)
    virtual OscillatableSpectrum AntiTauFromEComponent() = 0;

    /// Charged current tau neutrino appearance from muon neutrino (\f$\nu_\mu\to\nu_\tau\f$)
    virtual OscillatableSpectrum TauFromMuComponent() = 0;
    /// Charged current tau antineutrino appearance from muon antineutrino (\f$\bar\nu_\mu\to\bar\nu_\tau\f$)
    virtual OscillatableSpectrum AntiTauFromMuComponent() = 0;

    /// Neutral currents
    virtual Spectrum NCComponent() = 0;

    virtual void SaveTo(TDirectory* dir) const;
  };
}
