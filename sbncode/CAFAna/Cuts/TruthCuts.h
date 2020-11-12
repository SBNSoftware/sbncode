#pragma once

#include <cassert>
#include "CAFAna/Core/Cut.h"
#include "StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  const Cut kHasMatchedNu([](const caf::SRSliceProxy* slc)
                          {
                            return slc->truth.index >= 0;
                          });

  /// \brief Is this a Neutral %Current event?
  ///
  /// In this case the selection function is simple enough that we can include
  /// it inline as a lambda function.
  const Cut kIsNC([](const caf::SRSliceProxy* slc)
                  {
                    return !slc->truth.iscc;
                  });

  //----------------------------------------------------------------------
  /// Helper for defining true CC event cuts
  class CCFlavSel
  {
  public:
    CCFlavSel(int pdg, int pdgorig) : fPdg(pdg), fPdgOrig(pdgorig)
    {
    }

    bool operator()(const caf::SRSliceProxy* slc) const
    {
      return kHasMatchedNu(slc) &&
        slc->truth.iscc &&
        abs(slc->truth.initpdg) == fPdgOrig &&
        abs(slc->truth.pdg) == fPdg;
    }
  protected:
    int fPdg, fPdgOrig;
  };

  /// Helper for defining true CC event cuts
  class NCFlavOrig
  {
  public:
    NCFlavOrig(int pdgorig) : fPdgOrig(pdgorig)
    {
    }

    bool operator()(const caf::SRSliceProxy* slc) const
    {
      return kHasMatchedNu(slc) && slc->truth.isnc && abs(slc->truth.initpdg) == fPdgOrig;
    }
  protected:
    int fPdgOrig;
  };

  // Finally, the function argument to the Cut constructor can be a "functor"
  // object (one with operator()). This allows similar logic but with different
  // constants to be easily duplicated.

  /// Select CC \f$ \nu_\mu\to\nu_e \f$
  const Cut kIsNueApp (CCFlavSel(12, 14));
  /// Select CC \f$ \nu_\mu\to\nu_\mu \f$
  const Cut kIsNumuCC (CCFlavSel(14, 14));
  /// Select CC \f$ \nu_e\to\nu_e \f$
  const Cut kIsBeamNue(CCFlavSel(12, 12));
  /// Select CC \f$ \nu_e\to\nu_\mu \f$
  const Cut kIsNumuApp(CCFlavSel(14, 12));
  /// Select CC \f$ \nu_\mu\to\nu_\tau \f$
  const Cut kIsTauFromMu(CCFlavSel(16, 14));
  /// Select CC \f$ \nu_e\to\nu_\tau \f$
  const Cut kIsTauFromE(CCFlavSel(16, 12));


  // NC from Numu
  const Cut kIsNCFromNumu(NCFlavOrig(14));
  // NC from Nue
  const Cut kIsNCFromNue(NCFlavOrig(12));

  /// Is this truly an antineutrino?
  const Cut kIsAntiNu([](const caf::SRSliceProxy* slc)
                      {
                        return kHasMatchedNu(slc) && slc->truth.pdg < 0;
                      });
}
