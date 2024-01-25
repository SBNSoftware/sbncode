/**
 * @file FluxReaderAna_module.cc
 * @author Marco Del Tutto
 * @date 3 Dec 2021
 * @brief An analyzer module to make flux histograms
 *
 */

////////////////////////////////////////////////////////////////////////
// Class:       FluxReaderAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        FluxReaderAna_module.cc
//
// Generated at Mon Oct 11 19:10:14 2021 by Marco Del Tutto using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/TriggerResults.h"

#include "larcore/Geometry/Geometry.h"
#include "TGeoManager.h"
#include "TVector3.h"
#include "TTree.h"
#include "TDatabasePDG.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"

class FluxReaderAna;


class FluxReaderAna : public art::EDAnalyzer {
public:
  explicit FluxReaderAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FluxReaderAna(FluxReaderAna const&) = delete;
  FluxReaderAna(FluxReaderAna&&) = delete;
  FluxReaderAna& operator=(FluxReaderAna const&) = delete;
  FluxReaderAna& operator=(FluxReaderAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  virtual void beginSubRun(art::SubRun const& sr) override;

private:

  std::string _flux_label; ///< Label for flux dataproduct (to be set via fcl)
  std::string _flux_eventweight_multisim_producer; ///< Label for flux weights (to be set via fcl)
  float _x_shift; // cm (to be set via fcl)
  float _baseline; // cm (to be set via fcl)
  float _nu_intersection_z; ///< Z position where to evaluate the flux (to be set via fcl)
  float _nu_other_intersection_z; ///< Additional Z position where to evaluate the flux (to be set via fcl)
  std::vector<int> _requested_nu_pdgs; ///< Only save neutrinos with PDGs in this list (to be set via fcl)

  bool _apply_position_cuts; // wheter or not to apply position cuts (to be set via fcl)
  float _x_cut; // cm, only saves neutrinos with x in [-_xcut, +_xcut] (to be set via fcl)
  float _y_cut; // cm, only saves neutrinos with y in [-_ycut, +_ycut] (to be set via fcl)

  bool _add_systs; // wheter or not to add systs

  const std::vector<std::string> weight_names = { "expskin_Flux",
                                                  "horncurrent_Flux",
                                                  "kminus_Flux",
                                                  "kplus_Flux",
                                                  "kzero_Flux",
                                                  "nucleoninexsec_Flux",
                                                  "nucleonqexsec_Flux",
                                                  "nucleontotxsec_Flux",
                                                  "piminus_Flux",
                                                  "pioninexsec_Flux",
                                                  "pionqexsec_Flux",
                                                  "piontotxsec_Flux",
                                                  "piplus_Flux"
  };

  const TDatabasePDG *_pdg_database = TDatabasePDG::Instance();

  /// Returns the intersection point with the front face of the TPC
  TVector3 GetIntersection(TVector3 nu_pos, TVector3 nu_dir, float z_location=0);

  TTree* _tree;
  bool _nu_hit; /// True if the neutrino hit the requested volumes
  int _nu_pdg; /// PDG of neutrino
  float _nu_e; /// Energy of the neutrino
  float _nu_t; /// Time of the neutrino
  float _nu_w; /// Neutrino weight
  float _nu_x; /// X poisition of neutrino at the front face of the TPC
  float _nu_y; /// Y poisition of neutrino at the front face of the TPC
  float _nu_z; /// Z poisition of neutrino at the front face of the TPC
  float _nu_vx; /// X poisition of neutrino at neutrino production point
  float _nu_vy; /// Y poisition of neutrino at neutrino production point
  float _nu_vz; /// Z poisition of neutrino at neutrino production point
  float _nu_px; /// X momentum of neutrino
  float _nu_py; /// Y momentum of neutrino
  float _nu_pz; /// Z momentum of neutrino
  int _nu_decay; /// Neutrino parent decay code
  float _nu_dk2gen; /// Distance from decay to ray origin
  float _nu_p_angle; /// Angle between neutrino and parent direction
  int _nu_p_type; /// Neutrino parent pdg
  float _nu_p_dpx; /// Neutrino parent momentum x
  float _nu_p_dpy; /// Neutrino parent momentum x
  float _nu_p_dpz; /// Neutrino parent momentum x
  float _nu_p_e; /// Neutrino parent energy
  float _nu_r; /// Neutrino r
  float _nu_oaa; /// Neutrino off axis angle

  float _nu_other_x; /// X poisition of neutrino at the front face of the TPC (other)
  float _nu_other_y; /// Y poisition of neutrino at the front face of the TPC (other)
  float _nu_other_z; /// Z poisition of neutrino at the front face of the TPC (other)
  float _nu_other_t; /// Time of the neutrino (other)

  std::vector<float> _evtwgt_flux_oneweight; ///< Weights for FLUX reweighting (multisim) (combines all variations)
  std::vector<std::vector<float>> _evtwgt_flux_weights = std::vector<std::vector<float>>(weight_names.size(), std::vector<float>()); ///< Weights for FLUX reweighting (multisim)

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< Number of POTs per subrun
};


FluxReaderAna::FluxReaderAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{

  _flux_label = p.get<std::string>("FluxLabel", "flux");
  _flux_eventweight_multisim_producer = p.get<std::string>("FluxEventWeightProducer", "fluxweight");
  _requested_nu_pdgs = p.get<std::vector<int>>("RequestedNuPDGs", std::vector<int>());

  _baseline = p.get<float>("Baseline"); // cm, distance from detector to target position
  _x_shift = p.get<float>("XShift", 0.); // cm, detector shift along X w.r.t. beamline coordinate system

  _nu_intersection_z = p.get<float>("NuIntersectionZ", 0.);
  _nu_other_intersection_z = p.get<float>("NuOtherIntersectionZ");
  // 49000 is ICARUS location in SBND coordinate system (600 - 110)
  // 36000 is MicroBooNE location in SBND coordinate system (470 - 110)

  _apply_position_cuts = p.get<bool>("ApplyPositionCuts", false);
  _x_cut = p.get<float>("XCut"); // cm
  _y_cut = p.get<float>("YCut"); // cm

  _add_systs = p.get<bool>("AddSysts");

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("tree", "");
  _tree->Branch("nu_hit", &_nu_hit, "nu_hit/O");
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("nu_t", &_nu_t, "nu_t/F");
  _tree->Branch("nu_w", &_nu_w, "nu_w/F");
  _tree->Branch("nu_x", &_nu_x, "nu_x/F");
  _tree->Branch("nu_y", &_nu_y, "nu_y/F");
  _tree->Branch("nu_z", &_nu_z, "nu_z/F");
  _tree->Branch("nu_vx", &_nu_vx, "nu_vx/F");
  _tree->Branch("nu_vy", &_nu_vy, "nu_vy/F");
  _tree->Branch("nu_vz", &_nu_vz, "nu_vz/F");
  _tree->Branch("nu_px", &_nu_px, "nu_px/F");
  _tree->Branch("nu_py", &_nu_py, "nu_py/F");
  _tree->Branch("nu_pz", &_nu_pz, "nu_pz/F");
  _tree->Branch("nu_decay", &_nu_decay, "nu_decay/I");
  _tree->Branch("nu_dk2gen", &_nu_dk2gen, "nu_dk2gen/F");
  _tree->Branch("nu_p_angle", &_nu_p_angle, "nu_p_angle/F");
  _tree->Branch("nu_p_type", &_nu_p_type, "nu_p_type/I");
  _tree->Branch("nu_p_dpx", &_nu_p_dpx, "nu_p_dpx/F");
  _tree->Branch("nu_p_dpy", &_nu_p_dpy, "nu_p_dpy/F");
  _tree->Branch("nu_p_dpz", &_nu_p_dpz, "nu_p_dpz/F");
  _tree->Branch("nu_p_e", &_nu_p_e, "nu_p_e/F");
  _tree->Branch("nu_r", &_nu_r, "nu_r/F");
  _tree->Branch("nu_oaa", &_nu_oaa, "nu_oaa/F");

  _tree->Branch("nu_other_x", &_nu_other_x, "nu_other_x/F");
  _tree->Branch("nu_other_y", &_nu_other_y, "nu_other_y/F");
  _tree->Branch("nu_other_z", &_nu_other_z, "nu_other_z/F");
  _tree->Branch("nu_other_t", &_nu_other_t, "nu_other_t/F");

  _tree->Branch("evtwgt_flux_oneweight", "std::vector<float>", &_evtwgt_flux_oneweight);

  int weight_i = 0;

  for(auto const& name : weight_names)
    {
      _tree->Branch(Form("evtwgt_flux_weight_%s", name.c_str()), &_evtwgt_flux_weights[weight_i]);
      ++weight_i;
    }

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");
}

void FluxReaderAna::analyze(art::Event const& e)
{

  art::Handle<art::TriggerResults> trigger_h;
  e.getByLabel("TriggerResults", trigger_h);
  art::TriggerResults const& trigger = *trigger_h;
  if (!trigger.accept()) {
    std::cout << "[FluxReaderAna] Event was not accepted." << std::endl;
    return;
  }

  art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
  e.getByLabel(_flux_label, mcFluxHandle);
  std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;

  art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
  e.getByLabel(_flux_label, mctruthHandle);
  std::vector<simb::MCTruth> const& mclist = *mctruthHandle;

  art::FindManyP<sbn::evwgh::EventWeightMap> mct_to_fluxewm(mctruthHandle, e, _flux_eventweight_multisim_producer);

  for(unsigned int inu = 0; inu < mclist.size(); inu++) {
    simb::MCParticle nu = mclist[inu].GetNeutrino().Nu();
    simb::MCFlux flux = fluxlist[inu];

    if (_requested_nu_pdgs.size() &&
        std::find(_requested_nu_pdgs.begin(), _requested_nu_pdgs.end(), nu.PdgCode()) == _requested_nu_pdgs.end()) {
      std::cout << "[FluxReaderAna] Skipping neutrino with PDG " << nu.PdgCode() << std::endl;
      continue;
    }

    // std::cout << "This neutrino has vtx " << nu.Vx() << ", " << nu.Vy() << ", " << nu.Vz() << std::endl;
    // std::cout << "This neutrino has dir " << nu.Px() << ", " << nu.Py() << ", " << nu.Pz() << std::endl;

    // Calculate the position of the neutrino at Z = _nu_intersection_z
    TVector3 intersection = GetIntersection(TVector3(nu.Vx(),nu.Vy(),nu.Vz()),
                                            TVector3(nu.Px(),nu.Py(),nu.Pz()),
                                            _nu_intersection_z);
    // Calculate the distance from the original to the new position
    float dist = (TVector3(nu.Vx(),nu.Vy(),nu.Vz()) - intersection).Mag();

    // Do the same for an additional Z
    TVector3 intersection_other = GetIntersection(TVector3(nu.Vx(),nu.Vy(),nu.Vz()),
                                                  TVector3(nu.Px(),nu.Py(),nu.Pz()),
                                                  _nu_other_intersection_z);
    float other_dist = (TVector3(nu.Vx(),nu.Vy(),nu.Vz()) - intersection_other).Mag();

    float light_speed = TMath::C() / 1e7; // cm / ns

    _nu_pdg = nu.PdgCode();
    _nu_e = nu.E();
    _nu_t = nu.T() + dist / light_speed;
    _nu_w = flux.fnimpwt;
    _nu_x = intersection.X();
    _nu_y = intersection.Y();
    _nu_z = intersection.Z();
    _nu_vx = nu.Vx();
    _nu_vy = nu.Vy();
    _nu_vz = nu.Vz();
    _nu_px = nu.Px();
    _nu_py = nu.Py();
    _nu_pz = nu.Pz();
    _nu_decay = flux.fndecay;
    _nu_dk2gen = flux.fdk2gen;
    _nu_p_type = flux.fptype;
    _nu_p_dpx = flux.fpdpx;
    _nu_p_dpy = flux.fpdpy;
    _nu_p_dpz = flux.fpdpz;
    float parent_mass = _pdg_database->GetParticle(_nu_p_type)->Mass();
    _nu_p_e = std::sqrt(flux.fpdpx * flux.fpdpx +
                        flux.fpdpy * flux.fpdpy +
                        flux.fpdpz * flux.fpdpz +
                        parent_mass * parent_mass);
    _nu_p_angle = TVector3(_nu_px, _nu_py, _nu_pz).Angle(TVector3(_nu_p_dpx, _nu_p_dpy, _nu_p_dpz));
    _nu_r = std::sqrt((_nu_x - _x_shift) * (_nu_x - _x_shift) + _nu_y * _nu_y);
    _nu_oaa = std::atan(_nu_r / _baseline);

    _nu_other_x = intersection_other.X();
    _nu_other_y = intersection_other.Y();
    _nu_other_z = intersection_other.Z();
    _nu_other_t = nu.T() + other_dist / light_speed;

    if (_apply_position_cuts) {
      if (_nu_x < -_x_cut || _nu_x > _x_cut || _nu_y < -_y_cut || _nu_y > _y_cut) {
        continue;
      }
    }

    // Flux weights
    if (_add_systs) {
      std::vector<art::Ptr<sbn::evwgh::EventWeightMap>> flux_ewm_v = mct_to_fluxewm.at(inu);
      if (flux_ewm_v.size() != 1) {
        std::cout << "[FluxReaderAna] EventWeightMap of " << _flux_eventweight_multisim_producer << " bigger than 1?" << std::endl;
      }
      std::map<std::string, std::vector<float>> evtwgt_map = *(flux_ewm_v[0]);

      _evtwgt_flux_oneweight.clear();

      std::vector<float> previous_weights;
      std::vector<float> final_weights;

      int countFunc = 0;
      for(auto it : evtwgt_map) {
        std::string func_name = it.first;
        std::vector<float> weight_v = it.second;

        if (previous_weights.size() == 0) {
          previous_weights.resize(weight_v.size(), 1.);
          final_weights.resize(weight_v.size(), 1.);
        }

        _evtwgt_flux_weights[countFunc] = weight_v;

        countFunc++;

        // Construct a single weight
        std::transform(previous_weights.begin(), previous_weights.end(),
                       weight_v.begin(),
                       final_weights.begin(),
                       std::multiplies<float>());
        previous_weights = final_weights;
      }

      _evtwgt_flux_oneweight = final_weights;
    }

    _tree->Fill();
  }
}

TVector3 FluxReaderAna::GetIntersection(TVector3 nu_pos, TVector3 nu_dir, float z_location) {

  TVector3 plane_point(0, 0, z_location);
  TVector3 plane_normal(0, 0, 1);

  TVector3 diff = nu_pos - plane_point;
  double prod1 = diff.Dot(plane_normal);
  double prod2 = nu_dir.Dot(plane_normal);
  double prod3 = prod1 / prod2;
  return nu_pos - nu_dir * prod3;

}


void FluxReaderAna::beginSubRun(art::SubRun const& sr) {

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> pot_handle;
  sr.getByLabel(_flux_label, pot_handle);

  if (pot_handle.isValid()) {
    _sr_pot = pot_handle->totpot;
  } else {
    _sr_pot = 0.;
  }
  std::cout << "POT for this subrun: " << _sr_pot << std::endl;

  _sr_tree->Fill();

}

DEFINE_ART_MODULE(FluxReaderAna)
