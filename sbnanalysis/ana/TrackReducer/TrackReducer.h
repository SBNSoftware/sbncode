#ifndef _sbnana_TrackReducer_HH_
#define _sbnana_TrackReducer_HH_

#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/ProviderManager.hh"

#include "Data.h"

namespace sbnana {
class TrackReducer : public core::SelectionBase {
public:
  void Initialize(fhicl::ParameterSet* config=NULL);
  void Finalize() {}
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco);

private:
  Tracks fTracks; 

};

} // end namespace

#endif
