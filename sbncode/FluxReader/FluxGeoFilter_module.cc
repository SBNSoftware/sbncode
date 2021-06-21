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
#include "TGeoManager.h"
#include "TVector3.h"
#include "TTree.h"

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

  float _x_shift = 45.7; // cm
  float _baseline = 11000; // cm
  float _baseline_icarus = 60000; // cm

  /// Returns the intersection point with the front face of the TPC
  TVector3 GetIntersection(TVector3 nu_pos, TVector3 nu_dir, float z_location=0);

  TTree* _tree;
  bool _nu_hit; /// True if the neutrino hit the requested volumes
  int _nu_pdg; /// PDG of neutrino
  float _nu_e; /// Energy of neutrino
  float _nu_x; /// X poisition of neutrino at the front face of the TPC
  float _nu_y; /// Y poisition of neutrino at the front face of the TPC
  float _nu_z; /// Z poisition of neutrino at the front face of the TPC
  float _nu_px; /// X momentum of neutrino
  float _nu_py; /// Y momentum of neutrino
  float _nu_pz; /// Z momentum of neutrino
  float _nu_p_angle; /// Angle between neutrino and parent direction
  int _nu_p_type; /// Neutrino parent pdg
  float _nu_p_dpx; /// Neutrino parent momentum x
  float _nu_p_dpy; /// Neutrino parent momentum x
  float _nu_p_dpz; /// Neutrino parent momentum x
  float _nu_r; /// Neutrino r
  float _nu_oaa; /// Neutrino off axis angle

  float _nu_icarus_x; /// X poisition of neutrino at the front face of the TPC (ICARUS)
  float _nu_icarus_y; /// Y poisition of neutrino at the front face of the TPC (ICARUS)
  float _nu_icarus_z; /// Z poisition of neutrino at the front face of the TPC (ICARUS)
  float _nu_icarus_r; /// Neutrino r (ICARUS)
  float _nu_icarus_oaa; /// Neutrino off axis angle (ICARUS)
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

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("tree", "");
  _tree->Branch("nu_hit", &_nu_hit, "nu_hit/O");
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("nu_x", &_nu_x, "nu_x/F");
  _tree->Branch("nu_y", &_nu_y, "nu_y/F");
  _tree->Branch("nu_z", &_nu_z, "nu_z/F");
  _tree->Branch("nu_px", &_nu_px, "nu_px/F");
  _tree->Branch("nu_py", &_nu_py, "nu_py/F");
  _tree->Branch("nu_pz", &_nu_pz, "nu_pz/F");
  _tree->Branch("nu_p_angle", &_nu_p_angle, "nu_p_angle/F");
  _tree->Branch("nu_p_type", &_nu_p_type, "nu_p_type/I");
  _tree->Branch("nu_p_dpx", &_nu_p_dpx, "nu_p_dpx/F");
  _tree->Branch("nu_p_dpy", &_nu_p_dpy, "nu_p_dpy/F");
  _tree->Branch("nu_p_dpz", &_nu_p_dpz, "nu_p_dpz/F");
  _tree->Branch("nu_r", &_nu_r, "nu_r/F");
  _tree->Branch("nu_oaa", &_nu_oaa, "nu_oaa/F");

  _tree->Branch("nu_icarus_x", &_nu_icarus_x, "nu_icarus_x/F");
  _tree->Branch("nu_icarus_y", &_nu_icarus_y, "nu_icarus_y/F");
  _tree->Branch("nu_icarus_z", &_nu_icarus_z, "nu_icarus_z/F");
  _tree->Branch("nu_icarus_r", &_nu_icarus_r, "nu_icarus_r/F");
  _tree->Branch("nu_icarus_oaa", &_nu_icarus_oaa, "nu_icarus_oaa/F");
}

bool FluxGeoFilter::filter(art::Event& e)
{

  // Implementation of required member function here.
  bool result=false; //false filters out, true passes

  //check if neutrino goes through volTPCActive
  geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
  TGeoManager* rgeo=geom->ROOTGeoManager();

  //  std::vector< TGeoNode const * > geoNodeVec=geom->FindAllVolumes(volFlux);

  art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
  e.getByLabel("flux",mcFluxHandle);
  std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;

  art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
  e.getByLabel("flux",mctruthHandle);
  std::vector<simb::MCTruth> const& mclist = *mctruthHandle;

  for(unsigned int inu = 0; inu < mclist.size(); inu++){
    simb::MCParticle nu = mclist[inu].GetNeutrino().Nu();
    simb::MCFlux flux = fluxlist[inu];
    rgeo->SetCurrentPoint(nu.Vx(),nu.Vy(),nu.Vz());
    rgeo->SetCurrentDirection(nu.Px(),nu.Py(),nu.Pz());
    // std::cout << "This neutrino has vtx " << nu.Vx() << ", " << nu.Vy() << ", " << nu.Vz() << std::endl;
    // std::cout << "This neutrino has dir " << nu.Px() << ", " << nu.Py() << ", " << nu.Pz() << std::endl;
    TVector3 intersection = GetIntersection(TVector3(nu.Vx(),nu.Vy(),nu.Vz()),
                                            TVector3(nu.Px(),nu.Py(),nu.Pz()));

    TVector3 intersection_icarus = GetIntersection(TVector3(nu.Vx(),nu.Vy(),nu.Vz()),
                                                   TVector3(nu.Px(),nu.Py(),nu.Pz()),
                                                   49000.); // 600 - 110

    _nu_pdg = nu.PdgCode();
    _nu_e = nu.E();
    _nu_x = intersection.X();
    _nu_y = intersection.Y();
    _nu_z = intersection.Z();
    _nu_px = nu.Px();
    _nu_py = nu.Py();
    _nu_pz = nu.Pz();
    _nu_p_type = flux.fptype;
    _nu_p_dpx = flux.fpdpx;
    _nu_p_dpy = flux.fpdpy;
    _nu_p_dpz = flux.fpdpz;
    _nu_p_angle = TVector3(_nu_px, _nu_py, _nu_pz).Angle(TVector3(_nu_p_dpx, _nu_p_dpy, _nu_p_dpz));
    _nu_r = std::sqrt((_nu_x - _x_shift) * (_nu_x - _x_shift) + _nu_y * _nu_y);
    _nu_oaa = std::atan(_nu_r * _nu_r / _baseline);

    _nu_icarus_x = intersection_icarus.X();
    _nu_icarus_y = intersection_icarus.Y();
    _nu_icarus_z = intersection_icarus.Z();
    _nu_icarus_r = std::sqrt((_nu_icarus_x - _x_shift) * (_nu_icarus_x - _x_shift) + _nu_icarus_y * _nu_icarus_y);
    _nu_icarus_oaa = std::atan(_nu_icarus_r * _nu_icarus_r / _baseline_icarus);

    // if (_nu_r > 1500) {
    //   result = false;
    // } else {
    //   result = true;
    // }

    bool in_tpc_active_area = false;
    if (_nu_x >= -200 && _nu_x <= 200 && _nu_y >= -200 && _nu_y <= 200) {
      in_tpc_active_area = true;
    }

    TGeoNode* node=rgeo->FindNode();
    // const double *xyz = rgeo->GetCurrentPoint();
    // std::cout << "Start position (x, y, z) " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << std::endl;
    _nu_hit = false;
    while (node) {
      std::string volname=node->GetVolume()->GetName();
      rgeo->FindNextBoundary();
      node=gGeoManager->Step();
      if (fVolFlux.find(volname)!=fVolFlux.end()) {
        result=true;
        // const double *xyz = rgeo->GetCurrentPoint();
        // std::cout << "-> It intersects the active volume. E = " << nu.E() << std::endl;
        // std::cout << "-> It intersects the active volume. (x, y, z) " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << std::endl;
        // std::cout << "-> My intersection is (x, y, z) " << intersection.X() << " " << intersection.Y() << " " << intersection.Z() << " " << std::endl;
        // std::cout << "-> Step is  " << gGeoManager->GetStep() << std::endl;
        // std::cout << "------------------------------" << std::endl;
        _nu_hit = true;
        break;
      }
    }
    if (_nu_hit || in_tpc_active_area) {
      _tree->Fill();
    }
  }
  return result;
}

TVector3 FluxGeoFilter::GetIntersection(TVector3 nu_pos, TVector3 nu_dir, float z_location) {

  TVector3 plane_point(0, 0, z_location);
  TVector3 plane_normal(0, 0, 1);

  TVector3 diff = nu_pos - plane_point;
  double prod1 = diff.Dot(plane_normal);
  double prod2 = nu_dir.Dot(plane_normal);
  double prod3 = prod1 / prod2;
  return nu_pos - nu_dir * prod3;

}

DEFINE_ART_MODULE(FluxGeoFilter)
