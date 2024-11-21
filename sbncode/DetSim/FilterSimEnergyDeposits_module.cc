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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

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
  double ShiftX(double z);

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

double FilterSimEnergyDeposits::ShiftX(double z)
{
  //known as a "Hill equation" from fitting distribution of angle corrected tracks looking at the deflection as it gets closer to z = 0
  double a = -0.431172;
  double b = -2.52724;
  double c = 0.22043;
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

  TVector3 p1(fP1_x,fP1_y,fP1_z);
  TVector3 p2(fP2_x,fP2_y,fP2_z);

  for(auto const& sedep : simEDeps)
  {
    TVector3 dep(sedep.StartX(), sedep.StartY(), sedep.StartZ());
    if(!(isWithinVolume(p1, p2, dep)))
    {
      art::ServiceHandle<geo::Geometry const> fGeometry;
      geo::TPCID tpcid = fGeometry->PositionToTPCID(sedep.MidPoint());
      short int sign = 1;
      //if not in TPC, don't do anything to be sure
      if(bool(tpcid))
      {
	const geo::TPCGeo& tpcGeo = fGeometry->TPC(tpcid);
	sign = -1*tpcGeo.DetectDriftDirection(); //this gives us direction towards the anode, I want to shift towards the cathode so take the inverse of what I get 
      }
      else
      {
	std::cout << "TPC ID not found! Not performing a shift!" << std::endl;
	sign = 0;
      }
      const int numphotons = sedep.NumPhotons();
      const int numelectrons = sedep.NumElectrons();
      const double syratio = sedep.ScintYieldRatio();
      const double energy = sedep.Energy();
      //need to shift in -x for x positions > cathode position and +x for x positions < cathode position in each TPC    
      double shift_start_x = sedep.Start().X() + sign*ShiftX(std::abs(sedep.Start().Z()));
      const geo::Point_t start = {shift_start_x, sedep.Start().Y(), sedep.Start().Z()};
      double shift_end_x = sedep.End().X() + sign*ShiftX(std::abs(sedep.End().Z()));
      const geo::Point_t end = {shift_end_x, sedep.End().Y(), sedep.End().Z()};
      const double startT = sedep.StartT();
      const double endT = sedep.EndT();
      const int thisID = sedep.TrackID();
      const int thisPDG = sedep.PdgCode();
      const int origID = sedep.OrigTrackID();
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
  }
  e.put(std::move(pSimEDeps));
}

DEFINE_ART_MODULE(FilterSimEnergyDeposits)
