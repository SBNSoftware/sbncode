#include "OscCalculatorPMNSOpt.h"

#include <cassert>
#include <cstdlib>

namespace osc
{
  OscCalculatorPMNSOpt::OscCalculatorPMNSOpt()
    : fMixIdx(0), fDmIdx(0), fLRIdx(0)
  {
  }

  OscCalculatorPMNSOpt::~OscCalculatorPMNSOpt()
  {
    for(int i = 0; i < 2; ++i)
        for(auto it: fPMNSOpt[i])
          delete it.second.pmns;
  }

  IOscCalculatorAdjustable* OscCalculatorPMNSOpt::Copy() const
  {
    OscCalculatorPMNSOpt* ret = new OscCalculatorPMNSOpt(*this);
    // Raw pointers were blindly copied, so we'd be in trouble when the
    // destructors are called. More importantly, having two calculators sharing
    // one PMNS object will lead to both of them getting confused as to whether
    // they're up-to-date or not. Just clear out the cache entirely in the copy
    // and let it be repopulated.
    ret->fPMNSOpt[0].clear();
    ret->fPMNSOpt[1].clear();
    return ret;
  }

  double OscCalculatorPMNSOpt::P(int flavBefore, int flavAfter, double E)
  {
    // Normal usage of a calculator in a fit is to configure one set of
    // oscillation parameters and then calculate probabilities for all flavours
    // and for a range of energies. The underlying PMNSOpt object has to redo
    // its expensive diagonalization every time the energy changes, but not
    // when simply switching to a different flavour. Unfortunately callers in
    // CAFAna wind up making the calls in the wrong order to take advantage of
    // this. So, we maintain a map of calculators for different energies (also
    // for neutrinos and antineutrinos) and select the right one if available,
    // on the hope that it's already configured to do what we're asked, since
    // it might already have done a different flavour. This works because
    // callers often pass the same set of energies (bin centers) over and over
    // again.

    // If the caches get too large clear them. The caller is probably giving
    // the precise energy of each event?
    for(int i = 0; i < 2; ++i){
      if(fPMNSOpt[i].size() > 10000){
        for(auto it: fPMNSOpt[i]) delete it.second.pmns;
        fPMNSOpt[i].clear();
      }
    }

    const int anti = (flavBefore > 0) ? +1 : -1;
    assert(flavAfter/anti > 0);

    Val_t& calc = fPMNSOpt[(1+anti)/2][E];
    if(!calc.pmns) calc.pmns = new PMNSOpt;

    int i = -1, j = -1;
    if(abs(flavBefore) == 12) i = 0;
    if(abs(flavBefore) == 14) i = 1;
    if(abs(flavBefore) == 16) i = 2;
    if(abs(flavAfter) == 12) j = 0;
    if(abs(flavAfter) == 14) j = 1;
    if(abs(flavAfter) == 16) j = 2;
    assert(i >= 0 && j >= 0);

    const bool dirty = (fMixIdx > calc.mixIdx ||
                        fDmIdx > calc.dmIdx ||
                        fLRIdx > calc.lrIdx);

    // If the calculator has out of date parameters update them
    if(fMixIdx > calc.mixIdx){
      calc.pmns->SetMix(fTh12, fTh23, fTh13, fdCP);
      calc.mixIdx = fMixIdx;
    }
    if(fDmIdx > calc.dmIdx){
      calc.pmns->SetDeltaMsqrs(fDmsq21, fDmsq32);
      calc.dmIdx = fDmIdx;
    }

    if(dirty){
      // Cache results for all nine flavour combinations
      for(int ii = 0; ii < 3; ++ii){
        calc.pmns->ResetToFlavour(ii);
        // Assume Z/A=0.5
        const double Ne = fRho/2;
        calc.pmns->PropMatter(fL, E, Ne, anti);
        for(int jj = 0; jj < 3; ++jj){
          calc.P[ii][jj] = calc.pmns->P(jj);
        }
      }

      calc.lrIdx = fLRIdx;
    }

    return calc.P[i][j];
  }
} // namespace
