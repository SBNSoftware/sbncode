/**
 * \file MCParticleTreePrinter.cc
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
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "TH1D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class MCParticleTreePrinter
 * \brief Electron neutrino event selection
 */
class MCParticleTreePrinter : public core::SelectionBase {
public:
  /** Constructor. */
  MCParticleTreePrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fMCParticleTag = config ? config->get<std::string>("MCParticleTag", "largeant") : "largeant";
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
    std::cout << "New Event!\n";
    auto const &mcparticle_handle = ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);;
    const std::vector<simb::MCParticle> &mcparticles = *mcparticle_handle;
    art::FindManyP<simb::MCTruth, sim::GeneratedParticleInfo> particles_to_truth(mcparticle_handle, ev, fMCParticleTag);
    for (unsigned i = 0; i < mcparticles.size(); i++) {
      const simb::MCParticle &part = mcparticles[i];
      std::cout << "TrackID: " << part.TrackId() << " Mother: " << part.Mother() << " index: " << i<< std::endl;
      const sim::GeneratedParticleInfo &info = *particles_to_truth.data(i).at(0);
      const simb::MCTruth &mctruth = *particles_to_truth.at(i).at(0);
      std::cout << "Info index: " << info.generatedParticleIndex() << std::endl;
      if (info.generatedParticleIndex() < mctruth.NParticles()) {
        int gen_track_id = mctruth.GetParticle(info.generatedParticleIndex()).TrackId();
        std::cout << "Gen track ID: " << gen_track_id << std::endl;
      }
    }
    return false; 
  }

protected:
  std::string fMCParticleTag;
  unsigned event_ind;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::MCParticleTreePrinter)

