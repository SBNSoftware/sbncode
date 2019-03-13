#pragma once

#include <cassert>
#include "CAFAna/Core/Cut.h"
#include "StandardRecord/StandardRecord.h"

namespace ana
{

  /// \brief Is this a Neutral %Current event?
  ///
  /// We use uniform-initializer syntax to concisely pass the list of necessary
  /// branches. In this case the selection function is simple enough that we
  /// can include it inline as a lambda function.
  const Cut kIsNC({},
                  [](const caf::StandardRecord* sr)
                  {
                    return !sr->sbn.truth.neutrino[0].iscc;
                  });

  //----------------------------------------------------------------------
  /// Helper for defining true CC event cuts
  class CCFlavSel
  {
  public:
    CCFlavSel(int pdg, int pdgorig) : fPdg(pdg), fPdgOrig(pdgorig)
    {
    }

    bool operator()(const caf::StandardRecord* sr) const
    {
      return sr->sbn.truth.neutrino[0].iscc && abs(sr->sbn.truth.neutrino[0].pdg) == fPdgOrig
             && abs(sr->sbn.truth.neutrino[0].pdg) == fPdg;
    }
  protected:
    int fPdg, fPdgOrig;
  };

  // Finally, the function argument to the Cut constructor can be a "functor"
  // object (one with operator()). This allows similar logic but with different
  // constants to be easily duplicated.

  /// Select CC \f$ \nu_\mu\to\nu_e \f$
  const Cut kIsSig    ({}, CCFlavSel(12, 14));
  /// Select CC \f$ \nu_\mu\to\nu_\mu \f$
  const Cut kIsNumuCC ({}, CCFlavSel(14, 14));
  /// Select CC \f$ \nu_e\to\nu_e \f$
  const Cut kIsBeamNue({}, CCFlavSel(12, 12));
  /// Select CC \f$ \nu_e\to\nu_\mu \f$
  const Cut kIsNumuApp({}, CCFlavSel(14, 12));
  /// Select CC \f$ \nu_\mu\to\nu_\tau \f$
  const Cut kIsTauFromMu({}, CCFlavSel(16, 14));
  /// Select CC \f$ \nu_e\to\nu_\tau \f$
  const Cut kIsTauFromE({}, CCFlavSel(16, 12));

  /// Is this truly an antineutrino?
  const Cut kIsAntiNu({},
                      [](const caf::StandardRecord* sr)
                      {
                        return sr->sbn.truth.neutrino[0].pdg < 0;
                      });

  // const Cut kIsTrueFV({},
  //                     [](const caf::StandardRecord* sr)
  //                     {
		// 	if (sr->dune.isFD){
		// 	  return ( 
		// 		  abs(sr->dune.vtx_x) < 310 &&
		// 		  abs(sr->dune.vtx_y) < 550 &&
  // 				  sr->dune.vtx_z >  50      &&
  // 				  sr->dune.vtx_z < 1244 );
		// 	}
		// 	else{
		// 	  return (
		// 		  abs(sr->dune.vtx_x) < 300 &&
		// 		  abs(sr->dune.vtx_y) < 100 &&
		// 		  sr->dune.vtx_z > 50 &&
		// 		  sr->dune.vtx_z < 350
		// 		  );
		// 	}
		//       });

  //ETW 11/5/2018 Fiducial cut using MVA variable
  //Should use the previous one (kIsTrueFV) for nominal analysis
  // const Cut kPassFid_MVA({},
		// 	[](const caf::StandardRecord* sr)
		// 	{
		// 	  return ( sr->dune.mvanumu > -1 );
		// 	});

}
