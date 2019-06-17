#pragma once

#include "CAFAna/Core/Spectrum.h"

namespace ana
{
  /// Represent the ratio between two spectra
  class Ratio
  {
  public:
    friend class Spectrum;

    /// \param num Numerator of the ratio
    /// \param denom Denominator of the ratio
    /// \param purOrEffErrs Does this ratio represent a purity or efficiency
    ///                     plot? If so, error bars are calculated differently.
    Ratio(const Spectrum& num, const Spectrum& denom,
	  bool purOrEffErrs = false);

    /// Don't use this constructor unless you REALLY KNOW what you're doing.
    Ratio(TH1* h, std::string varName = "");
    virtual ~Ratio();

    Ratio(const Ratio& rhs);
    Ratio& operator=(const Ratio& rhs);

    Ratio& operator*=(const Ratio& rhs);
    Ratio operator*(const Ratio& rhs) const;

    Ratio& operator/=(const Ratio& rhs);
    Ratio operator/(const Ratio& rhs) const;

    TH1D* ToTH1(Color_t col = kBlack,
                Style_t style = kSolid) const;
  protected:
    TH1D* fHist;
  };

  inline Ratio operator/(const Spectrum& lhs, const Spectrum& rhs){return Ratio(lhs, rhs);}
} // namespace
