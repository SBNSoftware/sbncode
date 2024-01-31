////////////////////////////////////////////////////////////////////////
// Class:       AdjustSimForTrigger
// Plugin Type: producer (Unknown Unknown)
// File:        AdjustSimForTrigger_module.cc
//
// Generated at December 2023 by Bruce Howard (howard@fnal.gov) using cetskelgen.
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

#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/OpDetWaveform.h"
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
  art::InputTag fWaveformLabel;
};

AdjustSimForTrigger::AdjustSimForTrigger(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fInputSimEnergyDepositLabel{p.get<art::InputTag>("InputSimEnergyDepositLabel")}
  , fInputSimPhotonsLabel{p.get<art::InputTag>("InputSimPhotonsLabel")}
  , fInputTriggerLabel{p.get<art::InputTag>("InputTriggerLabel")}
  ,  fWaveformLabel(p.get<std::string>("WaveformLabel"))
{
  produces<std::vector<sim::SimEnergyDeposit>>();
  produces<std::vector<sim::SimPhotons>>();
  products<std::vector<raw::OpDetWaveform>>();
}

void AdjustSimForTrigger::produce(art::Event& e)
{
  auto const& simEDeps =
    e.getProduct<std::vector<sim::SimEnergyDeposit>>(fInputSimEnergyDepositLabel);

  auto const& simPhotons = e.getProduct<std::vector<sim::SimPhotons>>(fInputSimPhotonsLabel);

  auto const& triggers = e.getProduct<std::vector<raw::Trigger>>(fInputTriggerLabel);

  // Loop over the SimEnergyDeposit objects and shift time BACK by the TRIGGER time, assuming there is a trigger
  double timeShiftForTrigger = 0.;
  if (triggers.size() != 1) {
    throw art::Exception(art::errors::EventProcessorFailure)
      << "MORE THAN ONE TRIGGER IN EVENT... why?";
  }
  else if (triggers[0].TriggerTime() >
             (std::numeric_limits<double>::min() + std::numeric_limits<double>::epsilon()) &&
           triggers[0].TriggerTime() <
             (std::numeric_limits<double>::max() - std::numeric_limits<double>::epsilon())) {
    timeShiftForTrigger = 1000. * (triggers[0].TriggerTime() - triggers[0].BeamGateTime());
  }

  mf::LogInfo("AdjustSimForTrigger") << "FOR THIS EVENT THE TIME SHIFT BEING ASSUMED IS "
                                     << timeShiftForTrigger << " ns ..." << std::endl;

  auto simEDepVec = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  simEDepVec->reserve(simEDeps.size());

  for (auto const& inSimEDep : simEDeps) {
    const int numphotons = inSimEDep.NumPhotons();
    const int numelectrons = inSimEDep.NumElectrons();
    const double syratio = inSimEDep.ScintYieldRatio();
    const double energy = inSimEDep.Energy();
    const geo::Point_t start = {
      inSimEDep.Start().X(), inSimEDep.Start().Y(), inSimEDep.Start().Z()};
    const geo::Point_t end = {inSimEDep.End().X(), inSimEDep.End().Y(), inSimEDep.End().Z()};
    const double startT = inSimEDep.StartT() - timeShiftForTrigger;
    const double endT = inSimEDep.EndT() - timeShiftForTrigger;
    const int thisID = inSimEDep.TrackID();
    const int thisPDG = inSimEDep.PdgCode();
    const int origID = inSimEDep.OrigTrackID();

    simEDepVec->push_back(sim::SimEnergyDeposit(numphotons,
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

  // Repeat for SimPhotons
  auto simPhotonsVec = std::make_unique<std::vector<sim::SimPhotons>>(simPhotons);

  for (auto& photons : *simPhotonsVec) {
    for (auto& photon : photons) {
      photon.Time -= timeShiftForTrigger;
    }
  }

  e.put(std::move(simEDepVec));
  e.put(std::move(simPhotonsVec));
}

DEFINE_ART_MODULE(AdjustSimForTrigger)
