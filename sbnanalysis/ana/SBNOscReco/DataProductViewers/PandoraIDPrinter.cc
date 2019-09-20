/**
 * \file PandoraIDPrinter.cc
 *
 *
 * Author:
 */

#include <iostream>
#include <array>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/Experiment.hh"
#include "core/ProviderManager.hh"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class PandoraIDPrinter
 * \brief Electron neutrino event selection
 */
class PandoraIDPrinter : public core::SelectionBase {
public:
  /** Constructor. */
  PandoraIDPrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fPandoraTags = config ? config->get<std::vector<std::string>>("PandoraTags", {"pandora"}) : std::vector<std::string>(1,"pandora");
    fPandoraTrackTags = config ? config->get<std::vector<std::string>>("PandoraTags", {"pandoraTrack"}) : std::vector<std::string>(1,"pandoraTrack");
    event_ind = 0;
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {
  }

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
    for (unsigned j = 0; j < fPandoraTags.size(); j++) {
      auto const &tracks_handle = ev.getValidHandle<std::vector<recob::Track>>(fPandoraTrackTags[j]);;
      art::FindManyP<recob::PFParticle, void> tracks_to_particles(tracks_handle, ev, fPandoraTrackTags[j]);
      const std::vector<recob::Track> &tracks = *tracks_handle;
      std::cout << "TAG: " << fPandoraTrackTags[j] << std::endl;
      for (unsigned i = 0; i < tracks.size(); i++) {
        unsigned pfp_id = tracks_to_particles.at(i).at(0)->Self();
        std::cout << "Index: " << i << " track ID: " << tracks[i].ID() << " pfp id: " << pfp_id  <<std::endl;
      }
    }
    return false; 
  }

protected:
  std::vector<std::string> fPandoraTags;
  std::vector<std::string> fPandoraTrackTags;
  unsigned event_ind;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::PandoraIDPrinter)

