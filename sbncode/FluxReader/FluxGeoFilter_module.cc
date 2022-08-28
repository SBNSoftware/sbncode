////////////////////////////////////////////////////////////////////////
// Class:       FluxGeoFilter
// Plugin Type: filter (art v3_01_02)
// File:        FluxGeoFilter_module.cc
//
// Generated at Fri Apr  3 16:21:46 2020 by Zarko Pavlovic using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

//geometry
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "TGeoManager.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <memory>


class FluxGeoFilter;


class FluxGeoFilter : public art::EDFilter {
public:
  explicit FluxGeoFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FluxGeoFilter(FluxGeoFilter const&) = delete;
  FluxGeoFilter(FluxGeoFilter&&) = delete;
  FluxGeoFilter& operator=(FluxGeoFilter const&) = delete;
  FluxGeoFilter& operator=(FluxGeoFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:
  // Declare member data here.
  std::set<std::string> fVolFlux;
};

FluxGeoFilter::FluxGeoFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  std::cout<<"Configuring flux filter."<<std::endl;
  std::vector<std::string> tmp=p.get<std::vector<std::string> >("volumes");
  fVolFlux.insert(tmp.begin(),tmp.end());

  std::cout<<"Filtering flux through volumes: "<<std::endl;
  for (auto s: fVolFlux) std::cout<<"\t"<<s<<std::endl;
}

bool FluxGeoFilter::filter(art::Event& e)
{

  // Implementation of required member function here.
  bool result=false; //false filters out, true passes

  //check if neutrino goes through volTPCActive
  geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
  TGeoManager* rgeo=geom->ROOTGeoManager();

  // art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
  // e.getByLabel("flux",mcFluxHandle);
  // std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;

  art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
  e.getByLabel("flux",mctruthHandle);
  std::vector<simb::MCTruth> const& mclist = *mctruthHandle;

  for(unsigned int inu = 0; inu < mclist.size(); inu++){
    simb::MCParticle nu = mclist[inu].GetNeutrino().Nu();
    rgeo->SetCurrentPoint(nu.Vx(),nu.Vy(),nu.Vz());
    rgeo->SetCurrentDirection(nu.Px(),nu.Py(),nu.Pz());
    TGeoNode* node=rgeo->FindNode();
    while (node) {
      std::string volname=node->GetVolume()->GetName();
      rgeo->FindNextBoundary();
      node=gGeoManager->Step();
      if (fVolFlux.find(volname)!=fVolFlux.end()) {
        result=true;
        break;
      }
    }
  }
  return result;
}

DEFINE_ART_MODULE(FluxGeoFilter)
