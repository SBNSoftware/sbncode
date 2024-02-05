////////////////////////////////////////////////////////////////////////
// Class:       AdjustSimForTrigger
// Plugin Type: producer (Unknown Unknown)
// File:        AdjustSimForTrigger_module.cc
//
// Generated at December 2023 by Bruce Howard (howard@fnal.gov) using
// cetskelgen.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include <memory>

class AdjustSimForTrigger;

class AdjustSimForTrigger : public art::EDProducer {
public:
  explicit AdjustSimForTrigger(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AdjustSimForTrigger(AdjustSimForTrigger const&) = delete;
  AdjustSimForTrigger(AdjustSimForTrigger&&) = delete;
  AdjustSimForTrigger& operator=(AdjustSimForTrigger const&) = delete;
  AdjustSimForTrigger& operator=(AdjustSimForTrigger&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  art::InputTag fInputSimEnergyDepositLabel;
  art::InputTag fInputSimPhotonsLabel;
  art::InputTag fInputTriggerLabel;
  art::InputTag fInputWaveformLabel;
  art::InputTag fInputBeamGateInfoLabel;
};

AdjustSimForTrigger::AdjustSimForTrigger(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fInputSimEnergyDepositLabel{p.get<art::InputTag>("InputSimEnergyDepositLabel")}
  , fInputSimPhotonsLabel{p.get<art::InputTag>("InputSimPhotonsLabel")}
  , fInputTriggerLabel{p.get<art::InputTag>("InputTriggerLabel")}
  , fInputWaveformLabel(p.get<std::string>("InputWaveformLabel"))
  , fInputBeamGateInfoLabel{p.get<art::InputTag>("InputBeamGateInfoLabel")}
{
  produces<std::vector<sim::SimEnergyDeposit>>();
  produces<std::vector<sim::SimPhotons>>();
  produces<std::vector<raw::OpDetWaveform>>();
  produces<std::vector<sim::BeamGateInfo>>();
}

void AdjustSimForTrigger::produce(art::Event& e)
{
  auto const& simEDeps =
    e.getProduct<std::vector<sim::SimEnergyDeposit>>(fInputSimEnergyDepositLabel);
  auto const& simPhotons = e.getProduct<std::vector<sim::SimPhotons>>(fInputSimPhotonsLabel);
  auto const& triggers = e.getProduct<std::vector<raw::Trigger>>(fInputTriggerLabel);
  auto const& waveforms = e.getProduct<std::vector<raw::OpDetWaveform>>(fInputWaveformLabel);
  auto const& beamGates = e.getProduct<std::vector<sim::BeamGateInfo>>(fInputBeamGateInfoLabel);

  if (triggers.size() != 1) {
    if (triggers.empty()) {
      throw art::Exception(art::errors::EventProcessorFailure) << "NO TRIGGER IDENTIFIED!\n";
    }
    throw art::Exception(art::errors::EventProcessorFailure)
      << "MORE THAN ONE TRIGGER IN EVENT... why?\n";
  }

  if (beamGates.size() != 1) {
    if (beamGates.empty()) {
      throw art::Exception(art::errors::EventProcessorFailure) << "THERE IS NO BEAM GATE INFO!\n";
    }
    throw art::Exception(art::errors::EventProcessorFailure) << "MORE THAN ONE BEAM GATE?\n";
  }

  // Assuming there is a trigger, get time shift
  auto const& trigger = triggers[0];

  const bool hasValidTriggerTime =
    trigger.TriggerTime() >
      (std::numeric_limits<double>::min() + std::numeric_limits<double>::epsilon()) &&
    trigger.TriggerTime() <
      (std::numeric_limits<double>::max() - std::numeric_limits<double>::epsilon());

  const double timeShiftForTrigger_us =
    hasValidTriggerTime ? trigger.TriggerTime() - trigger.BeamGateTime() : 0.;
  const double timeShiftForTrigger_ns = 1000. * timeShiftForTrigger_us;

  mf::LogInfo("AdjustSimForTrigger")
    << "FOR THIS EVENT THE TIME SHIFT BEING ASSUMED IS " << timeShiftForTrigger_ns << " ns ...\n";

  // Loop over the SimEnergyDeposit objects and shift time BACK by the TRIGGER
  auto pSimEDeps = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  pSimEDeps->reserve(simEDeps.size());

  for (auto const& inSimEDep : simEDeps) {
    const int numphotons = inSimEDep.NumPhotons();
    const int numelectrons = inSimEDep.NumElectrons();
    const double syratio = inSimEDep.ScintYieldRatio();
    const double energy = inSimEDep.Energy();
    const geo::Point_t start = {
      inSimEDep.Start().X(), inSimEDep.Start().Y(), inSimEDep.Start().Z()};
    const geo::Point_t end = {inSimEDep.End().X(), inSimEDep.End().Y(), inSimEDep.End().Z()};
    const double startT = inSimEDep.StartT() - timeShiftForTrigger_ns;
    const double endT = inSimEDep.EndT() - timeShiftForTrigger_ns;
    const int thisID = inSimEDep.TrackID();
    const int thisPDG = inSimEDep.PdgCode();
    const int origID = inSimEDep.OrigTrackID();

    pSimEDeps->push_back(sim::SimEnergyDeposit(numphotons,
                                               numelectrons,
                                               syratio,
                                               energy,
                                               start,
                                               end,
                                               startT,
                                               endT,
                                               thisID,
                                               thisPDG,
                                               origID));
  }

  // Repeat for sim::BeamGateInfo
  const auto& beamGate = beamGates[0];

  const double shiftedBeamGateStart = beamGate.Start() - timeShiftForTrigger_ns;
  const double gateWidth = beamGate.Width();
  const sim::BeamType_t beam = beamGate.BeamType();

  const sim::BeamGateInfo shiftedBeamGate(shiftedBeamGateStart, gateWidth, beam);

  auto pBeamGateInfos = std::make_unique<std::vector<sim::BeamGateInfo>>(1, shiftedBeamGate);

  // Repeat for sim::SimPhotons
  auto pSimPhotonss = std::make_unique<std::vector<sim::SimPhotons>>(simPhotons);

  for (auto& photons : *pSimPhotonss) {
    for (auto& photon : photons) {
      photon.Time -= timeShiftForTrigger_ns;
    }
  }

  // Repeat for raw::OpDetWaveform
  auto pWaveforms = std::make_unique<std::vector<raw::OpDetWaveform>>(waveforms);

  for (auto& waveform : *pWaveforms) {
    waveform.SetTimeStamp(waveform.TimeStamp() - timeShiftForTrigger_us);
  }

  e.put(std::move(pSimEDeps));
  e.put(std::move(pBeamGateInfos));
  e.put(std::move(pSimPhotonss));
  e.put(std::move(pWaveforms));
}

DEFINE_ART_MODULE(AdjustSimForTrigger)
