////////////////////////////////////////////////////////////////////////
// Class:       FilterSimEnergyDeposits
// Plugin Type: producer (Unknown Unknown)
// File:        FilterSimEnergyDeposits_module.cc
//
// Generated at Thu Oct 17 20:16:12 2024 by Jacob Zettlemoyer using cetskelgen
// from cetlib version 3.18.02.
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

#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "TVector3.h"

#include <memory>

class FilterSimEnergyDeposits;


class FilterSimEnergyDeposits : public art::EDProducer {
public:
  explicit FilterSimEnergyDeposits(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilterSimEnergyDeposits(FilterSimEnergyDeposits const&) = delete;
  FilterSimEnergyDeposits(FilterSimEnergyDeposits&&) = delete;
  FilterSimEnergyDeposits& operator=(FilterSimEnergyDeposits const&) = delete;
  FilterSimEnergyDeposits& operator=(FilterSimEnergyDeposits&&) = delete;
  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

  art::InputTag fInitSimEnergyDepositLabel;
  std::string fOutSimEnergyDepositLabel;
  double fP1_x;
  double fP1_y;
  double fP1_z;
  double fP2_x;
  double fP2_y;
  double fP2_z;
  static constexpr auto& kModuleName = "FilterSimEnergyDeposit";

  bool isWithinVolume(TVector3 p1, TVector3 p2, TVector3 sedep);

};


FilterSimEnergyDeposits::FilterSimEnergyDeposits(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fInitSimEnergyDepositLabel{p.get<art::InputTag>("InitSimEnergyDepositLabel", "undefined")}
  , fOutSimEnergyDepositLabel{p.get<std::string>("OutSimEnergyDepositLabel", "undefined")}
  , fP1_x{p.get<double>("P1_X", 0.0)}
  , fP1_y{p.get<double>("P1_Y", 0.0)}
  , fP1_z{p.get<double>("P1_Z", 0.0)}
  , fP2_x{p.get<double>("P2_X", 0.0)}
  , fP2_y{p.get<double>("P2_Y", 0.0)}
  , fP2_z{p.get<double>("P2_Z", 0.0)}
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  produces<std::vector<sim::SimEnergyDeposit>>();
}

bool FilterSimEnergyDeposits::isWithinVolume(TVector3 p_min, TVector3 p_max, TVector3 edep)
{
  bool ingap = false;
  if(edep.X() >= p_min.X() && edep.X() <= p_max.X() &&
     edep.Y() >= p_min.Y() && edep.Y() <= p_max.Y() &&
     edep.Z() >= p_min.Z() && edep.Z() <= p_max.Z())
    {
      ingap = true;
    }
  return ingap;
}

void FilterSimEnergyDeposits::produce(art::Event& e)
{
  auto const& simEDeps =
    e.getProduct<std::vector<sim::SimEnergyDeposit>>(fInitSimEnergyDepositLabel);
  auto pSimEDeps = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  pSimEDeps->reserve(simEDeps.size());

  TVector3 p1(fP1_x,fP1_y,fP1_z);
  TVector3 p2(fP2_x,fP2_y,fP2_z);

  for(auto const& sedep : simEDeps)
  {
    TVector3 dep(sedep.StartX(), sedep.StartY(), sedep.StartZ());
    if(!(isWithinVolume(p1, p2, dep)))
    {
      pSimEDeps->emplace_back(sedep);
    }
  }
  e.put(std::move(pSimEDeps));
}

DEFINE_ART_MODULE(FilterSimEnergyDeposits)
