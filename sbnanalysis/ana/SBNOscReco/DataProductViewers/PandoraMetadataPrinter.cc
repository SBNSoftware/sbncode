/**
 * \file PandoraMetadataPrinter.cc
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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class PandoraMetadataPrinter
 * \brief Electron neutrino event selection
 */
class PandoraMetadataPrinter : public core::SelectionBase {
public:
  /** Constructor. */
  PandoraMetadataPrinter() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fPandoraTags = config ? config->get<std::vector<std::string>>("PandoraTags", {"pandora"}) : std::vector<std::string>(1,"pandora");
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
    std::cout << "\nNew Event!\n";
    for (unsigned j = 0; j < fPandoraTags.size(); j++) {
      auto const &particles_handle = ev.getValidHandle<std::vector<recob::PFParticle>>(fPandoraTags[j]); 
      art::FindManyP<larpandoraobj::PFParticleMetadata, void> particles_to_metadata(particles_handle, ev, fPandoraTags[j]);
      std::cout << "TAG: " << fPandoraTags[j] << std::endl;
      for (unsigned i = 0; i < particles_handle->size(); i++) {
        std::cout << "Particle PDG: " << particles_handle->at(i).PdgCode() << std::endl; 
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadatas = particles_to_metadata.at(i);
        if (metadatas.size() > 0) {
          assert(metadatas.size() == 1);
          const larpandoraobj::PFParticleMetadata &metadata = *metadatas[0];
          std::cout << "New PFParticle!\n";
          for (auto const &prop_pair: metadata.GetPropertiesMap()) {
            std::cout << prop_pair.first << " " << prop_pair.second << std::endl;
          }
        }
      }
    }
    return false; 
  }

protected:
  std::vector<std::string> fPandoraTags;
  unsigned event_ind;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::PandoraMetadataPrinter)

