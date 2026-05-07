/**
 * @file    sbncode/DetSim/AdjustSimForTrigger_module.cc
 * @author  Bruce Howard (howard@fnal.gov)
 * @date    December 2023
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimEnergyDepositLite.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "sbncode/Utilities/AssnsUtils.h" // sbn::RebindAssociatedProducts()

#include <ios> // std::boolalpha
#include <limits>
#include <memory>
#include <utility> // std::move()
#include <vector>


/**
 * @brief Applies a time shift to selected data products for simulation.
 * 
 * This module produces a selection of simulation data products shifting their
 * time reference to a new "event" (in the sense of something happening, not in
 * the DAQ sense).
 * 
 * The time point of the new reference event is specified in
 * @ref DetectorClocksElectronicsTime "electronics time scale".
 * Technically, this module applies a time shift so that the time of that
 * reference event becomes the previous reference value in electronics time scale.
 * 
 * More practically, we start from an electronics time scale where the reference
 * time is set at `DetectorClocksData::TriggerTime()` (typically `1500.0`
 * &micro;s).
 * That time reference may represent a hardware trigger, as the function name
 * suggest, or anything else. At this point we introduce a new trigger time
 * (from the simulation; e.g. `1502.0` &micro;s) and we want that time to become
 * the new reference (and "hardware trigger") of the electronics time scale.
 * In this example, this module will add a shift of -2 &micro;s to all the
 * selected data products, so that the time of the new trigger is `1500.0`
 * and all follows.
 * 
 * @note The new trigger time is going to be whatever value is in
 *       `DetectorClocksData::TriggerTime()`, rather than an hard-coded value
 *       like `1500.0`. Therefore some care needs to be taken in the workflow
 *       to make sure that at the time of the execution of this module
 *       `DetectorClocksService` is yielding the desired new-reference value.
 *       Also, after this module is called `DetectorClocksService` itself may
 *       report obsolete values. Specifically, `DetectorClocksServiceStandard`
 *       will report the old `BeamGateTime()`, which was from ether a
 *       configuration parameter or from a trigger data product.
 * 
 * This module reads the new reference time from a `raw::Trigger` object
 * (assumed to hold a time on the current electronics scale), and, unless
 * configured not to (`DropTriggerProduct`) it produces a new trigger data
 * product with the shifted trigger information. This data product is suitable
 * for configuring services like `DetectorClocksServiceStandard` for the
 * following stages of the workflow.
 * 
 * If the reference time is not valid, the applied shift is `0`.
 * A reference time is valid if it is neither the maximum nor the minimum value
 * of a double (`std::numeric_limits<double>::max()` and
 * `std::numeric_limits<double>::min()`).
 * 
 * @note This module does not provide indication on whether a shift was
 *       performed or not. The only way to know is to check the reference time
 *       input used and verify whether it's considered valid by this module
 *       according to the logic described above.
 * 
 * 
 * ### Optional shift
 * 
 * It is possible to specify a fixed additional shift to the time reference.
 * This may be useful for example if the simulation of the trigger yields
 * exactly the time at which the triggering activity is detected, but the
 * trigger hardware would take still some time in order to tag that time.
 * Note, however, that the additional shift is also applied to the beam gate
 * time.
 * This additional time is specified via the `AdditionalOffset` configuration
 * parameter.
 * If the new reference time is invalid, no shift is applied and this offset is
 * also ignored.
 * 
 * 
 * ### Multiple shifts on the same event
 * 
 * In general, using as input a sample that has already undergone a previous
 * shift is surprisingly robust: as long as the input is all consistently
 * shifted, offsets and delays that were already applied are not re-applied
 * (this was explicitly verified on `sim::SimPhotons`).
 * However, shifting is not "undone": if an input event was shifted and now the
 * new proposed time reference is invalid, the event will be applied no
 * additional shift, leaving it to the time reference after the first shift.
 * 
 * The chances of the old time reference being completely invalidated depend on
 * how the changes intervened between the first and the second shift.
 * The module _could_ be adapted to attempt a shifting reversal, if provided
 * with the necessary additional information on the first shift.
 * 
 * 
 * Input
 * ------
 * 
 * * `std::vector<raw::Trigger>` (`InputTriggerLabel`): the first of the
 *   triggers in the collection will be used as a new reference.
 *   If the trigger time value is not valid, no shift at all will be performed.
 *   However, the beam gate time is still expected to be valid.
 *   An empty collection is not allowed, even when there is no valid trigger.
 * 
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>`
 *   (`BindWaveformBaselines`): the association of the original optical detector
 *   waveforms to their baselines.
 * 
 * 
 * Output
 * -------
 * 
 * For each enabled data product to be shifted, the corresponding shifted data
 * product is produced, with elements in the same order as in the original
 * collections. Normally, no associations are ported on.
 * 
 * * `std::vector<sim::AuxDetSimChannel>` (if `ShiftAuxDetIDEs` is set):
 *   all `sim::AuxDetIDE` entry and exit times are shifted. Order of channels
 *   and IDE is preserved.
 * * `std::vector<sim::SimEnergyDeposit>` (if `ShiftSimEnergyDeposits` is set):
 *   all `sim::SimEnergyDeposit` entry and exit times are shifted.
 *   Order of deposits is preserved.
 * * `std::vector<sim::SimEnergyDepositLite>` (if `ShiftSimEnergyDepositLites`
 *   is set): all `sim::SimEnergyDepositLite` entry and exit times are shifted.
 *   Order of deposits is preserved.
 * * `std::vector<sim::SimPhotons>` (if `ShiftSimPhotons` is set):
 *   all `sim::SimPhotons` times in all channels are shifted.
 *   Order of channels and photons is preserved.
 * * `std::vector<raw::Trigger>` (unless `DropTriggerProduct` is set):
 *   a collection of shifted triggers is produced; the first one is the shifted
 *   version of the reference trigger, which can then be used as new trigger
 *   data product e.g. for `DetectorClocksServiceStandard`. If the reference
 *   trigger time is not valid, the trigger time will be overwritten: the
 *   trigger time will be set to the value from `DetectorClocksService`, and the
 *   beam gate time will be the same as the input trigger object.
 * * `std::vector<sim::BeamGateInfo>` (if `ShiftBeamGateInfo` is set):
 *   the beam gate used by the event generator, shifted. Generator times and
 *   particles from the detector simulation (GEANT4) are not shifted, but pretty
 *   much everything else may be, including scintillation photons and energy
 *   depositions. Depending on which aspect of the simulation is being
 *   investigated, either the unshifted (input) or shifted (output of this
 *   module) gate needs to be used.
 * * `std::vector<raw::OpDetWaveform>` (if `ShiftWaveforms` is set): the input
 *   PMT waveforms, from which presumably the new trigger was extracted, are
 *   shifted by simply modifying their timestamps.
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>`: enabled only
 *   if the optical waveforms are being shifted _and_ an association data
 *   product name is specified (`BindWaveformBaselines`), it rebinds the
 *   original waveform baselines to the shifted waveforms. Note that neither the
 *   baseline nor the order nor the content of the waveforms change: only their
 *   start time (and implicitly the end time) does.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * * `InputTriggerLabel` (input tag, mandatory): the tag of the trigger data
 *   product with the new reference time. It must be available and not empty.
 * * `ShiftAuxDetIDE` (bool, default: `false`): enable the shifting of auxiliary
 *   detector simulation data product at `InitAuxDetSimChannelLabel`.
 * * `InitAuxDetSimChannelLabel` (input tag): tag of the auxiliary detector
 *   simulation `sim::AuxDetSimChannel` product to be shifted.
 * * `ShiftBeamGateInfo` (bool, default: `false`): enable the shifting of
 *   beam gate data product at `InitBeamGateInfoLabel`.
 * * `InitBeamGateInfoLabel` (input tag): tag of the simulation beam gate data
 *   product to be shifted. This can be produced by a LArSoft generation module
 *   or by a trigger module.
 * * `ShiftSimEnergyDeposits` (bool, default: `false`): enable the shifting of
 *   full energy deposit data product at `InitSimEnergyDepositLabel`.
 * * `InitSimEnergyDepositLabel` (input tag): tag of the simulated energy
 *   deposition data product to be shifted.
 * * `ShiftSimEnergyDepositLites` (bool, default: `false`): enable the shifting
 *   of lightweight energy deposit data product at
 *   `InitSimEnergyDepositLiteLabel`.
 * * `InitSimEnergyDepositLiteLabel` (input tag): tag of the lightweight
 *   simulated energy deposition data product to be shifted. This reduced
 *   version is typically kept around for tracking back to truth information.
 * * `ShiftSimPhotons` (bool, default: `false`): enable the shifting of
 *   scintillation photon data product at `InitSimPhotonsLabel`.
 * * `InitSimPhotonsLabel` (input tag): tag of the simulated scintillation photon
 *   data product to be shifted.
 * * `ShiftWaveforms` (bool, default: `false`): enable the shifting of optical
 *   detector waveform data product at `InitWaveformLabel`.
 * * `InitWaveformLabel` (input tag): tag of the simulated optical detector
 *   waveform data product to be shifted. These waveforms may already be
 *   available if the worflow extracted the trigger time (new reference) out of
 *   them.
 * * `BindWaveformBaselines` (input tag, default: `""`): tag of the association
 *   between the optical detector waveforms being shifted and their baselines.
 *   If empty (default), baseline associations will not be produced.
 * * `AdditionalOffset` (real value, default: `0`): additional offset in
 *   microseconds to be added to the new reference trigger.
 * * `DropTriggerProduct` (flag, default: `false`): if set, no shifted trigger
 *   data product will be produced.
 * 
 */
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
  art::InputTag fInitSimEnergyDepositLiteLabel;
  art::InputTag fInitSimPhotonsLabel;
  art::InputTag fInitWaveformLabel;
  bool fShiftAuxDetIDEs;
  bool fShiftBeamGateInfo;
  bool fShiftSimEnergyDeposits;
  bool fShiftSimEnergyDepositLites;
  bool fShiftSimPhotons;
  bool fShiftWaveforms;
  art::InputTag fBindWaveformBaselines; ///< Tag of OpDetWaveform-baseline associations to be rebound.
  double fAdditionalOffset;
  bool fDropTriggerProduct; ///< Do not put the shifted trigger data product into the event.
  static constexpr auto& kModuleName = "AdjustSimForTrigger";
};

AdjustSimForTrigger::AdjustSimForTrigger(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fInputTriggerLabel{p.get<art::InputTag>("InputTriggerLabel")}
  , fInitAuxDetSimChannelLabel(p.get<art::InputTag>("InitAuxDetSimChannelLabel", "undefined"))
  , fInitBeamGateInfoLabel{p.get<art::InputTag>("InitBeamGateInfoLabel", "undefined")}
  , fInitSimEnergyDepositLabel{p.get<art::InputTag>("InitSimEnergyDepositLabel", "undefined")}
  , fInitSimEnergyDepositLiteLabel{p.get<art::InputTag>("InitSimEnergyDepositLiteLabel", "undefined")}
  , fInitSimPhotonsLabel{p.get<art::InputTag>("InitSimPhotonsLabel", "undefined")}
  , fInitWaveformLabel(p.get<art::InputTag>("InitWaveformLabel", "undefined"))
  , fShiftAuxDetIDEs{p.get<bool>("ShiftAuxDetIDEs", false)}
  , fShiftBeamGateInfo{p.get<bool>("ShiftBeamGateInfo", false)}
  , fShiftSimEnergyDeposits{p.get<bool>("ShiftSimEnergyDeposits", false)}
  , fShiftSimEnergyDepositLites{p.get<bool>("ShiftSimEnergyDepositLites", false)}
  , fShiftSimPhotons{p.get<bool>("ShiftSimPhotons", false)}
  , fShiftWaveforms{p.get<bool>("ShiftWaveforms", false)}
  , fBindWaveformBaselines{p.get<art::InputTag>("BindWaveformBaselines", "")}
  , fAdditionalOffset{p.get<double>("AdditionalOffset", 0.)}
  , fDropTriggerProduct{p.get<bool>("DropTriggerProduct", false)}
{
  if (!(fShiftSimEnergyDeposits || fShiftSimPhotons || fShiftWaveforms || fShiftAuxDetIDEs ||
        fShiftBeamGateInfo || fShiftSimEnergyDepositLites)) {
    throw art::Exception(art::errors::EventProcessorFailure)
      << kModuleName << ": NO SHIFTS ENABLED!\n";
  }
  bool const doWaveformBaselines = fShiftWaveforms && !fBindWaveformBaselines.empty();
  mf::LogInfo(kModuleName) << std::boolalpha << "SHIFTING AUXDETIDES? " << fShiftAuxDetIDEs << '\n'
                           << "SHIFTING BEAMGATEINFO? " << fShiftBeamGateInfo << '\n'
                           << "SHIFTING SIMENERGYDEPOSITS? " << fShiftSimEnergyDeposits << '\n'
                           << "SHIFTING SIMENERGYDEPOSITLITES? " << fShiftSimEnergyDepositLites << '\n'
                           << "SHIFTING SIMPHOTONS? " << fShiftSimPhotons << '\n'
                           << "SHIFTING OPDETWAVEFORMS? " << fShiftWaveforms << '\n'
                           << "   ASSNS OPDETWAVEFORM-BASELINES? " << doWaveformBaselines
                           << (doWaveformBaselines? (" ('" + fBindWaveformBaselines.encode() + ")"): "");

  if (!fDropTriggerProduct) produces<std::vector<raw::Trigger>>();
  if (fShiftAuxDetIDEs) { produces<std::vector<sim::AuxDetSimChannel>>(); }
  if (fShiftBeamGateInfo) { produces<std::vector<sim::BeamGateInfo>>(); }
  if (fShiftSimEnergyDeposits) { produces<std::vector<sim::SimEnergyDeposit>>(); }
  if (fShiftSimEnergyDepositLites) { produces<std::vector<sim::SimEnergyDepositLite>>(); }
  if (fShiftSimPhotons) { produces<std::vector<sim::SimPhotons>>(); }
  if (fShiftWaveforms) { 
    produces<std::vector<raw::OpDetWaveform>>(); 
    if (!fBindWaveformBaselines.empty())
      produces<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>(); 
  }
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

  const double newReferenceTime =
    hasValidTriggerTime ? trigger.TriggerTime() - fAdditionalOffset : clock_data.TriggerTime();
  const double timeShiftForTrigger_us = clock_data.TriggerTime() - newReferenceTime;
  const double timeShiftForTrigger_ns = 1000. * timeShiftForTrigger_us;

  mf::LogInfo(kModuleName) << "FOR THIS EVENT THE TIME SHIFT BEING ASSUMED IS "
                           << timeShiftForTrigger_ns << " ns ...";

  // Shifted trigger (beam gate info, optional, is later)
  // sbn::ExtraTriggerInfo and raw::ExternalTrigger are not shifted (so far);
  // it's debatable if they should be, since they only hold absolute timestamps
  if (!fDropTriggerProduct) {
    auto pShiftedTriggers = std::make_unique<std::vector<raw::Trigger>>();
    for (raw::Trigger const& unshiftedTrigger: triggers) { // ok, we required there is just one
      pShiftedTriggers->emplace_back(
        unshiftedTrigger.TriggerNumber(),
        hasValidTriggerTime  // trigger_time
          ? (unshiftedTrigger.TriggerTime() + timeShiftForTrigger_us)
          : clock_data.TriggerTime(),
        unshiftedTrigger.BeamGateTime() + timeShiftForTrigger_us,// beamgate_time
        unshiftedTrigger.TriggerBits()
        );
    } // for all triggers
    e.put(std::move(pShiftedTriggers));
  } // if produce shifted trigger

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
  // and SimEnergyDepositLite's
  if (fShiftSimEnergyDepositLites) {
    auto const& simEDepLites =
      e.getProduct<std::vector<sim::SimEnergyDepositLite>>(fInitSimEnergyDepositLiteLabel);

    auto pSimEDepLites = std::make_unique<std::vector<sim::SimEnergyDepositLite>>();

    for (auto const& inSimEDepLite : simEDepLites) {
      double energy = inSimEDepLite.Energy(); 
      geo::Point_t middlePos = inSimEDepLite.Position();
      double middleTime = inSimEDepLite.Time() + timeShiftForTrigger_ns;
      int ID = inSimEDepLite.TrackID();
      
      pSimEDepLites->emplace_back(energy, middlePos, middleTime, ID);
    }

    e.put(std::move(pSimEDepLites));
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
    
    if (!fBindWaveformBaselines.empty()) {
      // given that the shifting is one-to-one, rebinding is just replacing
      // each existing waveform pointer with one to the new waveform in the same position
      auto const& waveformBaselineAssns
       = e.getProduct<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>(fBindWaveformBaselines);
      art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr{ e };
      
      auto pWaveformBaselineAssns
        = std::make_unique<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>
        (sbn::RebindAssociatedProducts(waveformBaselineAssns, makeWaveformPtr));
      
      e.put(std::move(pWaveformBaselineAssns));
    } // if rebinding associations
    
  }
}

DEFINE_ART_MODULE(AdjustSimForTrigger)
