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
  art::InputTag fInputTriggerLabel;
  art::InputTag fInitAuxDetSimChannelLabel;
  art::InputTag fInitBeamGateInfoLabel;
  art::InputTag fInitSimEnergyDepositLabel;
  art::InputTag fInitSimPhotonsLabel;
  art::InputTag fInitWaveformLabel;
  bool fShiftAuxDetIDEs;
  bool fShiftBeamGateInfo;
  bool fShiftSimEnergyDeposits;
  bool fShiftSimPhotons;
  bool fShiftWaveforms;
  double fAdditionalOffset;
  static constexpr auto& kModuleName = "AdjustSimForTrigger";
};

AdjustSimForTrigger::AdjustSimForTrigger(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fInputTriggerLabel{p.get<art::InputTag>("InputTriggerLabel", "undefined")}
  , fInitAuxDetSimChannelLabel(p.get<art::InputTag>("InitAuxDetSimChannelLabel", "undefined"))
  , fInitBeamGateInfoLabel{p.get<art::InputTag>("InitBeamGateInfoLabel", "undefined")}
  , fInitSimEnergyDepositLabel{p.get<art::InputTag>("InitSimEnergyDepositLabel", "undefined")}
  , fInitSimPhotonsLabel{p.get<art::InputTag>("InitSimPhotonsLabel", "undefined")}
  , fInitWaveformLabel(p.get<art::InputTag>("InitWaveformLabel", "undefined"))
  , fShiftAuxDetIDEs{p.get<bool>("ShiftAuxDetIDEs", false)}
  , fShiftBeamGateInfo{p.get<bool>("ShiftBeamGateInfo", false)}
  , fShiftSimEnergyDeposits{p.get<bool>("ShiftSimEnergyDeposits", false)}
  , fShiftSimPhotons{p.get<bool>("ShiftSimPhotons", false)}
  , fShiftWaveforms{p.get<bool>("ShiftWaveforms", false)}
  , fAdditionalOffset{p.get<double>("AdditionalOffset", 0.)}
{
  if (!(fShiftSimEnergyDeposits || fShiftSimPhotons || fShiftWaveforms || fShiftAuxDetIDEs ||
        fShiftBeamGateInfo)) {
    throw art::Exception(art::errors::EventProcessorFailure)
      << kModuleName << ": NO SHIFTS ENABLED!\n";
  }
  mf::LogInfo(kModuleName) << std::boolalpha << "SHIFTING AUXDETIDES? " << fShiftAuxDetIDEs << '\n'
                           << "SHIFTING BEAMGATEINFO? " << fShiftBeamGateInfo << '\n'
                           << "SHIFTING SIMENERGYDEPOSITS? " << fShiftSimEnergyDeposits << '\n'
                           << "SHIFTING SIMPHOTONS? " << fShiftSimPhotons << '\n'
                           << "SHIFTING OPDETWAVEFORMS? " << fShiftWaveforms;

  if (fShiftAuxDetIDEs) { produces<std::vector<sim::AuxDetSimChannel>>(); }
  if (fShiftBeamGateInfo) { produces<std::vector<sim::BeamGateInfo>>(); }
  if (fShiftSimEnergyDeposits) { produces<std::vector<sim::SimEnergyDeposit>>(); }
  if (fShiftSimPhotons) { produces<std::vector<sim::SimPhotons>>(); }
  if (fShiftWaveforms) { produces<std::vector<raw::OpDetWaveform>>(); }
}

void AdjustSimForTrigger::produce(art::Event& e)
{
  auto const& clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const& triggers = e.getProduct<std::vector<raw::Trigger>>(fInputTriggerLabel);

  if (triggers.size() != 1) {
    if (triggers.empty()) {
      throw art::Exception(art::errors::EventProcessorFailure)
        << kModuleName << ": NO TRIGGER IDENTIFIED!";
    }
    throw art::Exception(art::errors::EventProcessorFailure)
      << kModuleName << ": MORE THAN ONE TRIGGER IN EVENT... why?";
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

  mf::LogInfo(kModuleName) << "FOR THIS EVENT THE TIME SHIFT BEING ASSUMED IS "
                           << timeShiftForTrigger_ns << " ns ...";

  // Loop over the sim::AuxDetIDE and shift time BACK by the TRIGGER
  if (fShiftAuxDetIDEs) {
    auto const& simChannels =
      e.getProduct<std::vector<sim::AuxDetSimChannel>>(fInitAuxDetSimChannelLabel);
    auto pSimChannels = std::make_unique<std::vector<sim::AuxDetSimChannel>>();

    pSimChannels->reserve(simChannels.size());

    for (auto const& simChannel : simChannels) {
      std::vector<sim::AuxDetIDE> shiftedAuxDetIDEs = simChannel.AuxDetIDEs();
      for (auto& auxDetIDE : shiftedAuxDetIDEs) {
        auxDetIDE.entryT += timeShiftForTrigger_ns;
        auxDetIDE.exitT += timeShiftForTrigger_ns;
      }
      pSimChannels->emplace_back(sim::AuxDetSimChannel(
        simChannel.AuxDetID(), shiftedAuxDetIDEs, simChannel.AuxDetSensitiveID()));
    }
    e.put(std::move(pSimChannels));
  }

  // Repeat for sim::BeamGateInfo
  if (fShiftBeamGateInfo) {
    auto const& beamGates = e.getProduct<std::vector<sim::BeamGateInfo>>(fInitBeamGateInfoLabel);

    if (beamGates.size() != 1) {
      if (beamGates.empty()) {
        throw art::Exception(art::errors::EventProcessorFailure)
          << kModuleName << ": THERE IS NO BEAM GATE INFO!\n";
      }
      throw art::Exception(art::errors::EventProcessorFailure)
        << kModuleName << ": MORE THAN ONE BEAM GATE?\n";
    }

    const auto& beamGate = beamGates[0];

    const double shiftedBeamGateStart = beamGate.Start() + timeShiftForTrigger_ns;
    const double gateWidth = beamGate.Width();
    const sim::BeamType_t beam = beamGate.BeamType();

    const sim::BeamGateInfo shiftedBeamGate(shiftedBeamGateStart, gateWidth, beam);

    auto pBeamGateInfos = std::make_unique<std::vector<sim::BeamGateInfo>>(1, shiftedBeamGate);

    e.put(std::move(pBeamGateInfos));
  }

  // Repeat for sim::SimEnergyDeposit
  if (fShiftSimEnergyDeposits) {
    auto const& simEDeps =
      e.getProduct<std::vector<sim::SimEnergyDeposit>>(fInitSimEnergyDepositLabel);

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
    auto const& simPhotons = e.getProduct<std::vector<sim::SimPhotons>>(fInitSimPhotonsLabel);
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
    auto const& waveforms = e.getProduct<std::vector<raw::OpDetWaveform>>(fInitWaveformLabel);
    auto pWaveforms = std::make_unique<std::vector<raw::OpDetWaveform>>(waveforms);

    for (auto& waveform : *pWaveforms) {
      waveform.SetTimeStamp(waveform.TimeStamp() + timeShiftForTrigger_us);
    }
    e.put(std::move(pWaveforms));
  }
}

DEFINE_ART_MODULE(AdjustSimForTrigger)
