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
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include <lardata/DetectorInfoServices/DetectorClocksService.h>
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
  art::InputTag fInputAuxDetSimChannelLabel;
  art::InputTag fInputBeamGateInfoLabel;
  double fAdditionalOffset;
  bool fShiftSimEnergyDeposits;
  bool fShiftSimPhotons;
  bool fShiftWaveforms;
  bool fShiftAuxDetIDEs;
  bool fShiftBeamGateInfo;
};

AdjustSimForTrigger::AdjustSimForTrigger(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fInputSimEnergyDepositLabel{p.get<art::InputTag>("InputSimEnergyDepositLabel", "")}
  , fInputSimPhotonsLabel{p.get<art::InputTag>("InputSimPhotonsLabel", "")}
  , fInputTriggerLabel{p.get<art::InputTag>("InputTriggerLabel", "")}
  , fInputWaveformLabel(p.get<art::InputTag>("InputWaveformLabel", ""))
  , fInputAuxDetSimChannelLabel(p.get<art::InputTag>("InputAuxDetSimChannelLabel", ""))
  , fInputBeamGateInfoLabel{p.get<art::InputTag>("InputBeamGateInfoLabel", "")}
  , fAdditionalOffset{p.get<double>("AdditionalOffset")}
  , fShiftSimEnergyDeposits{p.get<bool>("ShiftSimEnergyDeposits", false)}
  , fShiftSimPhotons{p.get<bool>("ShiftSimPhotons", false)}
  , fShiftWaveforms{p.get<bool>("ShiftWaveforms", false)}
  , fShiftAuxDetIDEs{p.get<bool>("ShiftAuxDetIDEs", false)}
  , fShiftBeamGateInfo{p.get<bool>("ShiftBeamGateInfo", false)}
{
  if (!(fShiftSimEnergyDeposits | fShiftSimPhotons | fShiftWaveforms | fShiftAuxDetIDEs | fShiftBeamGateInfo)) {
    throw art::Exception(art::errors::EventProcessorFailure)
      << "NO SHIFTS ENABLED!\n"
      << "SHIFTING SIMENERGYDEPOSITS? " << fShiftSimEnergyDeposits << '\n'
      << "SHIFTING SIMPHOTONS? " << fShiftSimPhotons << '\n'
      << "SHIFTING OPDETWAVEFORMS? " << fShiftWaveforms << '\n'
      << "SHIFTING AUXDETIDES? " << fShiftAuxDetIDEs << '\n'
      << "SHIFTING BEAMGATEINFO? " << fShiftBeamGateInfo << '\n';
  }

  if (fShiftSimEnergyDeposits) { produces<std::vector<sim::SimEnergyDeposit>>(); }
  if (fShiftSimPhotons) { produces<std::vector<sim::SimPhotons>>(); }
  if (fShiftWaveforms) { produces<std::vector<raw::OpDetWaveform>>(); }
  if (fShiftAuxDetIDEs) { produces<std::vector<sim::AuxDetSimChannel>>(); }
  if (fShiftBeamGateInfo) { produces<std::vector<sim::BeamGateInfo>>(); }
}

void AdjustSimForTrigger::produce(art::Event& e)
{
  auto const& clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const& triggers = e.getProduct<std::vector<raw::Trigger>>(fInputTriggerLabel);

  if (triggers.size() != 1) {
    if (triggers.empty()) {
      throw art::Exception(art::errors::EventProcessorFailure) << "NO TRIGGER IDENTIFIED!\n";
    }
    throw art::Exception(art::errors::EventProcessorFailure)
      << "MORE THAN ONE TRIGGER IN EVENT... why?\n";
  }

  // Assuming there is a trigger, get time shift
  auto const& trigger = triggers[0];

  const bool hasValidTriggerTime =
    trigger.TriggerTime() >
      (std::numeric_limits<double>::min() + std::numeric_limits<double>::epsilon()) &&
    trigger.TriggerTime() <
      (std::numeric_limits<double>::max() - std::numeric_limits<double>::epsilon());

  const double timeShiftForTrigger_us =
    hasValidTriggerTime ? clock_data.TriggerTime() - trigger.TriggerTime() + fAdditionalOffset : 0.;
  const double timeShiftForTrigger_ns = 1000. * timeShiftForTrigger_us;

  mf::LogInfo("AdjustSimForTrigger")
    << "FOR THIS EVENT THE TIME SHIFT BEING ASSUMED IS " << timeShiftForTrigger_ns << " ns ...\n";

  // Loop over the SimEnergyDeposit objects and shift time BACK by the TRIGGER
  if (fShiftSimEnergyDeposits) {
    auto const& simEDeps =
      e.getProduct<std::vector<sim::SimEnergyDeposit>>(fInputSimEnergyDepositLabel);

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
      const double startT = inSimEDep.StartT() + timeShiftForTrigger_ns;
      const double endT = inSimEDep.EndT() + timeShiftForTrigger_ns;
      const int thisID = inSimEDep.TrackID();
      const int thisPDG = inSimEDep.PdgCode();
      const int origID = inSimEDep.OrigTrackID();

      pSimEDeps->emplace_back(sim::SimEnergyDeposit(numphotons,
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
    e.put(std::move(pSimEDeps));
  }

  // Repeat for sim::SimPhotons
  if (fShiftSimPhotons) {
    auto const& simPhotons = e.getProduct<std::vector<sim::SimPhotons>>(fInputSimPhotonsLabel);
    auto pSimPhotonss = std::make_unique<std::vector<sim::SimPhotons>>(simPhotons);

    for (auto& photons : *pSimPhotonss) {
      for (auto& photon : photons) {
        photon.Time += timeShiftForTrigger_ns;
      }
    }
    e.put(std::move(pSimPhotonss));
  }

  // Repeat for raw::OpDetWaveform
  if (fShiftWaveforms) {
    auto const& waveforms = e.getProduct<std::vector<raw::OpDetWaveform>>(fInputWaveformLabel);
    auto pWaveforms = std::make_unique<std::vector<raw::OpDetWaveform>>(waveforms);

    for (auto& waveform : *pWaveforms) {
      waveform.SetTimeStamp(waveform.TimeStamp() + timeShiftForTrigger_us);
    }
    e.put(std::move(pWaveforms));
  }

  // Repeat for sim::AuxDetIDE
  if (fShiftAuxDetIDEs) {
    auto const& simChannels = e.getProduct<std::vector<sim::AuxDetSimChannel>>(fInputAuxDetSimChannelLabel);
    auto pSimChannels = std::make_unique<std::vector<sim::AuxDetSimChannel>>();

    pSimChannels->reserve(simChannels.size());

    for (auto const& simChannel : simChannels) {
      std::vector<sim::AuxDetIDE> shiftedAuxDetIDEs = simChannel.AuxDetIDEs();
      for (auto& auxDetIDE : shiftedAuxDetIDEs) {
        auxDetIDE.entryT += timeShiftForTrigger_ns;
        auxDetIDE.exitT += timeShiftForTrigger_ns;
      }
      pSimChannels->emplace_back(sim::AuxDetSimChannel(simChannel.AuxDetID(),
                                 shiftedAuxDetIDEs,
                                 simChannel.AuxDetSensitiveID()));
    }
    e.put(std::move(pSimChannels));
  }

  // Repeat for sim::BeamGateInfo
  if (fShiftBeamGateInfo) {
    auto const& beamGates = e.getProduct<std::vector<sim::BeamGateInfo>>(fInputBeamGateInfoLabel);

    if (beamGates.size() != 1) {
      if (beamGates.empty()) {
        throw art::Exception(art::errors::EventProcessorFailure) << "THERE IS NO BEAM GATE INFO!\n";
      }
      throw art::Exception(art::errors::EventProcessorFailure) << "MORE THAN ONE BEAM GATE?\n";
    }

    const auto& beamGate = beamGates[0];

    const double shiftedBeamGateStart = beamGate.Start() + timeShiftForTrigger_ns;
    const double gateWidth = beamGate.Width();
    const sim::BeamType_t beam = beamGate.BeamType();

    const sim::BeamGateInfo shiftedBeamGate(shiftedBeamGateStart, gateWidth, beam);

    auto pBeamGateInfos = std::make_unique<std::vector<sim::BeamGateInfo>>(1, shiftedBeamGate);

    e.put(std::move(pBeamGateInfos));
  }
}

DEFINE_ART_MODULE(AdjustSimForTrigger)
