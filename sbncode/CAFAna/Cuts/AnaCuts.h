/*
#pragma once

#include <cassert>
#include "CAFAna/Core/Cut.h"
#include "StandardRecord/StandardRecord.h"

namespace ana
{

  const Cut kPassFD_CVN_NUE({},
                  [](const caf::StandardRecord* sr)
                  {
                    return (sr->dune.cvnnue > 0.85 && sr->dune.cvnnumu < 0.5);
                  });

  const Cut kPassFD_CVN_NUMU({},
                  [](const caf::StandardRecord* sr)
                  {
                    return (sr->dune.cvnnumu > 0.5 && sr->dune.cvnnue < 0.85);
                  });

  const Cut kPassND_FHC_NUMU({},
                  [](const caf::StandardRecord* sr)
                  {
                    return (
			    sr->dune.reco_numu && 
			    (sr->dune.muon_contained || sr->dune.muon_tracker) &&
			    sr->dune.reco_q == -1 && 
			    sr->dune.Ehad_veto<30);
		      });

    const Cut kPassND_RHC_NUMU({},
                  [](const caf::StandardRecord* sr)
                  {
                    return (
			    sr->dune.reco_numu && 
			    (sr->dune.muon_contained || sr->dune.muon_tracker) &&
			    sr->dune.reco_q == +1 && 
			    sr->dune.Ehad_veto<30);
                  });


}
*/
