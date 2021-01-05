#include "SBNAna/Vars/TruthVars.h"
#include "SBNAna/Vars/Vars.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <cmath>

namespace ana
{
    const Var kHasTruthMatch(
      [](const caf::SRSliceProxy *slc)
      {
        return ( slc->truth.index != -1);
      }
      );

    const Var kTruthEnergy(
      [](const caf::SRSliceProxy *slc)
      {
        return ( kHasTruthMatch(slc) ? (float)slc->truth.E : -5.f );
      }
      );

    const Var kTruthVtxX(
      [](const caf::SRSliceProxy *slc)
      {
        return ( kHasTruthMatch(slc) ? (float)slc->truth.position.x : -999.f );
      }
      );

    const Var kTruthVtxY(
      [](const caf::SRSliceProxy *slc)
      {
        return ( kHasTruthMatch(slc) ? (float)slc->truth.position.y : -999.f );
      }
      );

    const Var kTruthVtxZ(
      [](const caf::SRSliceProxy *slc)
      {
        return ( kHasTruthMatch(slc) ? (float)slc->truth.position.z : -999.f );
      }
      );

    const Var kTruthVtxDistX(
      [](const caf::SRSliceProxy *slc)
      {
        return ( kHasTruthMatch(slc) ? (float)std::abs( slc->truth.position.x - kSlcVtxX(slc) ) : -999.f );
      }
      );

    const Var kTruthVtxDistY(
      [](const caf::SRSliceProxy *slc)
      {
        return ( kHasTruthMatch(slc) ? (float)std::abs( slc->truth.position.y - kSlcVtxY(slc) ) : -999.f );
      }
      );

    const Var kTruthVtxDistZ(
      [](const caf::SRSliceProxy *slc)
      {
        return ( kHasTruthMatch(slc) ? (float)std::abs( slc->truth.position.z - kSlcVtxZ(slc) ) : -999.f );
      }
      );

    const Var kTruthVtxDistMag(
      [](const caf::SRSliceProxy *slc)
      {
        return ( kHasTruthMatch(slc) ? (float)std::hypot( kTruthVtxDistX(slc), kTruthVtxDistY(slc), kTruthVtxDistZ(slc) ) : -5.f );
      }
      );
}
