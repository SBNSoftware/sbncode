#ifndef MERGESIMSOURCESLITE_H
#define MERGESIMSOURCESLITE_H

#include "larsim/MergeSimSources/MergeSimSources.h"
#include "lardataobj/MCBase/MCParticleLite.h"

namespace sbn {
  class MergeSimSourcesLiteUtility : sim::MergeSimSourcesUtility {

    public:
      MergeSimSourcesLiteUtility(std::vector<int> const&);

      void MergeMCParticleLites(std::vector<sim::MCParticleLite>&,
                                const std::vector<sim::MCParticleLite>&,
                                size_t);
    private:
      std::vector<int> fG4TrackIDOffsets;
  };
}
#endif
