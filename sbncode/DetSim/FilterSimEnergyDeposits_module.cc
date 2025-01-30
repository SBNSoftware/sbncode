/**
 * @file   sbncode/DetSim/FilterSimEnergyDeposits_module.cc
 * @date   October 17, 2024
 * @author Jacob Zettlemoyer
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

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
  geo::BoxBoundedGeo fBox;
  static constexpr auto kModuleName = "FilterSimEnergyDeposit";
  double ShiftX(double z) const;
  double fA;
  double fB;
  double fC;

};


FilterSimEnergyDeposits::FilterSimEnergyDeposits(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fInitSimEnergyDepositLabel{p.get<art::InputTag>("InitSimEnergyDepositLabel", "undefined")}
  , fBox{p.get<double>("P1_X"), p.get<double>("P2_X"),
	 p.get<double>("P1_Y"), p.get<double>("P2_Y"),
	 p.get<double>("P1_Z"), p.get<double>("P2_Z")}
  , fA{p.get<double>("A")}
  , fB{p.get<double>("B")}
  , fC{p.get<double>("C")}
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  produces<std::vector<sim::SimEnergyDeposit>>();
}


double FilterSimEnergyDeposits::ShiftX(double z) const
{
  //known as a "Hill equation" from fitting distribution of angle corrected tracks looking at the deflection as it gets closer to z = 0
  double a = fA;
  double b = fB;
  double c = fC;
  double num = a*std::pow(z, b);
  double denom = c + std::pow(z, b);
  return num/denom;
}

void FilterSimEnergyDeposits::produce(art::Event& e)
{
  auto const& simEDeps =
    e.getProduct<std::vector<sim::SimEnergyDeposit>>(fInitSimEnergyDepositLabel);
  auto pSimEDeps = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
  pSimEDeps->reserve(simEDeps.size());

  for(auto const& sedep : simEDeps)
  {
    if(fBox.ContainsPosition(sedep.Start()))
      continue;
    art::ServiceHandle<geo::Geometry const> geom;
    geo::TPCGeo const* TPC = geom->PositionToTPCptr(sedep.MidPoint());
    if(!TPC) {
      mf::LogVerbatim(kModuleName) << "TPC ID not found! Not performing a shift!";
    }
    const int numphotons = sedep.NumPhotons();
    const int numelectrons = sedep.NumElectrons();
    const double syratio = sedep.ScintYieldRatio();
    const double energy = sedep.Energy();
    geo::Point_t start = sedep.Start();
    geo::Point_t end = sedep.End();
    if (TPC) {
      auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();
      auto const& referencePlane = wireReadout.FirstPlane(TPC->ID());
      referencePlane.DriftPoint(start, ShiftX(std::abs(sedep.Start().Z())));
      referencePlane.DriftPoint(end, ShiftX(std::abs(sedep.End().Z())));
    }
    const double startT = sedep.StartT();
    const double endT = sedep.EndT();
    const int thisID = sedep.TrackID();
    const int thisPDG = sedep.PdgCode();
    const int origID = sedep.OrigTrackID();
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
  e.put(std::move(pSimEDeps));
}

DEFINE_ART_MODULE(FilterSimEnergyDeposits)
