#include "../Data/RecoEvent.h"
#include "core/Event.hh"
#include "Cuts.h"

namespace ana::SBNOsc {
  void SetEvent(numu::RecoEvent &event, const event::Event &core, const ana::SBNOsc::Cuts &cuts, numu::MCType file_type, bool true_particle_id=true);
}
