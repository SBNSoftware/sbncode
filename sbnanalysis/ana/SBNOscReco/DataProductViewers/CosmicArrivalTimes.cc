/**
 * \file CosmicArrivalTimes.cc
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

#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class CosmicArrivalTimes
 * \brief Electron neutrino event selection
 */
class CosmicArrivalTimes : public core::SelectionBase {
public:
  /** Constructor. */
  CosmicArrivalTimes() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fMCParticleTag = config ? config->get<std::string>("MCParticleTag", "largeant") : "largeant";
    fCosmicEnterTimes = new TH1D("cosmic_enter", "cosmic_enter", 700, -60., 10.); 
    fCosmicEnterTimesPE = new TH2D("cosmic_enter_v_pe", "cosmic_enter_v_pe", 700, -60., 10., 100, 0, 10000);
    event_ind = 0;
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {
    fOutputFile->cd();
    fCosmicEnterTimes->Write();
    fCosmicEnterTimesPE->Write();
  }

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& event, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
    std::cout << "New Event!\n";
    auto const &intime = event.getValidHandle<std::vector<simb::MCTruth>>({"GenInTimeSorter", "intime"});
    //auto const &outtime = event.getValidHandle<std::vector<simb::MCTruth>>({"GenInTimeSorter", "outtime"});

    art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> intime_to_particles(intime, event, fMCParticleTag);
    //art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> outtime_to_particles(outtime, event, fMCParticleTag);

    double in_time_photon_energy = 0.;
    unsigned in_time_photons = 0;
    const std::vector<sim::SimPhotons> &photon_list = *event.getValidHandle<std::vector<sim::SimPhotons>>("larg4intime");
    for (const sim::SimPhotons &photons: photon_list) {
      for (const sim::OnePhoton &photon: photons) {
        if (photon.Time > -200. && photon.Time < 1800.) {
          in_time_photon_energy += photon.Energy;
          in_time_photons ++;
        }
      }
    }
    gallery::Handle<std::vector<sim::SimPhotons>> reflected;
    event.getByLabel({"larg4intime", "Reflected"}, reflected);
    if (reflected.isValid()) {
    std::cout << "Reflected!\n";
    const std::vector<sim::SimPhotons> &photon_list = *event.getValidHandle<std::vector<sim::SimPhotons>>({"larg4intime", "Reflected"});
    for (const sim::SimPhotons &photons: photon_list) {
      for (const sim::OnePhoton &photon: photons) {
        if (photon.Time > -200. && photon.Time < 1800.) {
          in_time_photon_energy += photon.Energy;
          in_time_photons ++;
        }
      }
    }
    }

    for (const art::Ptr<simb::MCParticle> part: intime_to_particles.at(0)) {
      if (abs(part->PdgCode()) != 13) continue;
      unsigned n_traj = part->NumberTrajectoryPoints();
      bool next = false;
      for (unsigned i = 0; i < n_traj; i++) {
        for (const geo::BoxBoundedGeo &vol: fActiveVolumes) {
          if (vol.ContainsPosition(part->Position(i).Vect())) {
            fCosmicEnterTimes->Fill(part->Position(i).T() / 1000.) /* ns -> us*/;
            fCosmicEnterTimesPE->Fill(part->Position(i).T() / 1000., in_time_photons); 
            next = true; 
            break;
          }
        }
        if (next) break;
      }
    }

    /*
    for (const art::Ptr<simb::MCParticle> part: outtime_to_particles.at(0)) {
      unsigned n_traj = part->NumberTrajectoryPoints();
      bool next = false;
      for (unsigned i = 0; i < n_traj; i++) {
        for (const geo::BoxBoundedGeo &vol: fActiveVolumes) {
          if (vol.ContainsPosition(part->Position(i).Vect())) {
            fCosmicEnterTimes->Fill(part->Position(i).T() / 1000.);
            next = true; 
            break;
          }
        }
        if (next) break;
      }
    }*/

    return false; 
  }

protected:
  std::string fMCParticleTag;
  unsigned event_ind;
  TH1D *fCosmicEnterTimes;
  TH2D *fCosmicEnterTimesPE;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::CosmicArrivalTimes)

