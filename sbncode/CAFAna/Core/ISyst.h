#pragma once

#include <string>

#include "StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  /// \brief Encapsulate code to systematically shift a \ref
  /// caf::StandardRecord
  ///
  /// The Shift() function alters the \ref caf::StandardRecord or the weight
  /// associated with the event.
  class ISyst
  {
  public:
    ISyst(const std::string& shortName,
          const std::string& latexName,
	  bool applyPenalty = true,
	  double min = -3,
	  double max = +3);
    ISyst(const ISyst &) = delete;   // no copying.
    ISyst(ISyst && rhs) = delete;    // no moving either.
    virtual ~ISyst();

    void operator=(const ISyst &) = delete;  // still no copying.
    void operator=(ISyst &&)      = delete;  // etc.

    /// The name printed out to the screen
    virtual std::string ShortName() const final {return fShortName;}

    /// The name used on plots (ROOT's TLatex syntax)
    virtual std::string LatexName() const final {return fLatexName;}

    virtual double Penalty(double x) const;

    /// Should a penalty be applied for this shift?
    virtual bool ApplyPenalty() const {return fApplyPenalty;}

    /// Return the min/max value for this syst
    virtual double Min() const{return fMin;}
    virtual double Max() const{return fMax;}

    /// \brief Perform the systematic shift
    ///
    /// \param sigma   Number of sigma to shift record by
    /// \param sr      The record to inspect and alter
    /// \param weight  Scale this weight for reweighting systematics
    virtual void Shift(double sigma,
                       caf::SRSliceProxy* sr,
                       double& weight) const = 0;

    /// PredictionInterp normally interpolates between spectra made at
    /// +/-1,2,3sigma. For some systematics that's overkill. Override this
    /// function to specify different behaviour for this systematic.
    virtual int PredInterpMaxNSigma() const
    {
      return 3;
    }

  private:
    std::string fShortName;
    std::string fLatexName;
    bool fApplyPenalty;
    double fMin;
    double fMax;
  };

} // namespace
