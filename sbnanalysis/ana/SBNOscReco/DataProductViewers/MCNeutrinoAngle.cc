/**
 * \file MCNeutrinoAngle.cc
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
 * \class MCNeutrinoAngle
 * \brief Electron neutrino event selection
 */
class MCNeutrinoAngle : public core::SelectionBase {
public:
  /** Constructor. */
  MCNeutrinoAngle() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fLeptonAngle = new TH1D("lepton_angle", "lepton_angle", 100, 0, 1);
    fGenTag = config ? config->get<std::string>("GenTag", "generator") : "generator";
    event_ind = 0;
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {
    fOutputFile->cd();
    fLeptonAngle->Write();
  }

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
    const simb::MCTruth &mctruth = ev.getValidHandle<std::vector<simb::MCTruth>>(fGenTag)->at(0);
    if (mctruth.GetNeutrino().CCNC() == 0) {
      fLeptonAngle->Fill(mctruth.GetNeutrino().Lepton().Momentum().Vect().CosTheta());
    }
    return false; 
  }

protected:
  std::string fGenTag;
  TH1D *fLeptonAngle;
  unsigned event_ind;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::MCNeutrinoAngle)

