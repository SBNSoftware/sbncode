#ifndef SBN_ANALYSIS_CALORIMETRYANALYSIS_CXX
#define SBN_ANALYSIS_CALORIMETRYANALYSIS_CXX

#include <iostream>
#include <map>

#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"


#include "lardataobj/AnalysisBase/T0.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcorealg/GeoAlgo/GeoAlgo.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// #include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"


namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       CalorimetryAnalysis
// File:        CalorimetryAnalysis.cc
//
//              A basic analysis example
//
// Configuration parameters:
//
// TBD
//
// Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
//
////////////////////////////////////////////////////////////////////////

// declare helper functions
void TrkDirectionAtXYZ(const recob::Track trk, const double x, const double y, const double z, float out[3]);
void True2RecoMappingXYZ(float t, float x, float y, float z, float out[3]);
float XYZtoPlanecoordinate(const float x, const float y, const float z, const int plane);
float distance3d(const float& x1, const float& y1, const float& z1,
                  const float& x2, const float& y2, const float& z2);
float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                      const std::vector<geoalgo::AABox> &boxes);
void FillHits(const detinfo::DetectorClocksData &clockData,
              const std::vector<art::Ptr<recob::Hit>> &hits,
              std::vector<float> &charge_u,
              std::vector<float> &charge_v,
              std::vector<float> &charge_y,
              std::vector<unsigned> &wire_u,
              std::vector<unsigned> &wire_v,
              std::vector<unsigned> &wire_y,
              std::vector<unsigned> &channel_u,
              std::vector<unsigned> &channel_v,
              std::vector<unsigned> &channel_y,
              std::vector<int> &multiplicity_u,
              std::vector<int> &multiplicity_v,
              std::vector<int> &multiplicity_y,
              std::vector<float> &width_u,
              std::vector<float> &width_v,
              std::vector<float> &width_y,
              std::vector<float> &time_u,
              std::vector<float> &time_v,
              std::vector<float> &time_y,
              std::vector<unsigned> &nhit_u,
              std::vector<unsigned> &nhit_v,
              std::vector<unsigned> &nhit_y,
              bool use_integral=true);

class CalorimetryAnalysis : public art::EDAnalyzer
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  CalorimetryAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~CalorimetryAnalysis(){};

  /**
     * @brief Analysis function
     */
  void analyze(art::Event const &e) override;

  /**
     * @brief Fill Default info for event associated to neutrino
     */
  void fillDefault();

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree);

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree);

  void respondToOpenInputFile(const art::FileBlock& fb) override {
    (void) fb;
    _fileno ++;
  }

private:

  // function that given a track and its calo info fills NTuple variables
  void FillCalorimetry(art::Event const &e,
            const std::vector<art::Ptr<recob::Wire>> &allWires,
            const std::vector<art::Ptr<recob::Hit>> &allHits,
            const art::Ptr<recob::PFParticle> &pfp,
            const art::Ptr<recob::Track> &trk,
            const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints,
            const std::vector<art::Ptr<anab::Calorimetry>> &calos,
            const std::vector<art::Ptr<anab::ParticleID>> &pid,
            const std::vector<art::Ptr<recob::Hit>> &areaHits,
            const std::vector<art::Ptr<recob::Hit>> &trkHits,
            const std::vector<art::Ptr<recob::Hit>> &caloHits,
            const recob::MCSFitResult &muon_mcs,
            const sbn::RangeP &muon_range,
            const sbn::RangeP &proton_range,
            const bool fData,
            const std::vector<std::pair<int, float>> &particleMatches,
            const art::Ptr<simb::MCParticle> &trueParticle,
            const std::array<std::map<unsigned, std::vector<const sim::IDE*>>, 3> particle_ide_map,
            const std::array<std::vector<const sim::IDE*>,3> &allIDEs,
            const std::vector<geo::BoxBoundedGeo> &activeVolumes);


  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

  art::InputTag fPFPproducer;
  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;
  art::InputTag fTRKproducer;
  art::InputTag fT0producer;
  art::InputTag fHitFilterproducer;
  art::InputTag fAreaHitproducer;
  std::vector<art::InputTag> fHitProducers;
  art::InputTag fWireProducer;
  std::string fMCSproducer;
  std::string fRangeproducer;
  bool fAllTrueEnergyDeposits;

  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;

  TTree* _calo_tree;

  int _fileno;
  int _run, _sub, _evt;
  int _pfpno;
  int _trkid;
  // backtracking information
  int _backtracked_pdg;            // PDG code of backtracked particle
  float _backtracked_e;            // energy of backtracked particle
  float _backtracked_deposited_e_u;            // energy of backtracked particle
  float _backtracked_deposited_e_v;            // energy of backtracked particle
  float _backtracked_deposited_e_y;            // energy of backtracked particle
  float _backtracked_purity;       // purity of backtracking
  float _backtracked_completeness; // completeness of backtracking
  float _backtracked_overlay_purity; // purity of overlay
  float _backtracked_q_u;
  float _backtracked_q_v;
  float _backtracked_q_y;
  float _backtracked_qtotal_u;
  float _backtracked_qtotal_v;
  float _backtracked_qtotal_y;

  std::vector<float> _backtracked_pc_q_u;
  std::vector<float> _backtracked_pc_q_v;
  std::vector<float> _backtracked_pc_q_y;

  std::vector<float> _backtracked_pc_e_u;
  std::vector<float> _backtracked_pc_e_v;
  std::vector<float> _backtracked_pc_e_y;

  std::vector<unsigned> _backtracked_c_u;
  std::vector<unsigned> _backtracked_c_v;
  std::vector<unsigned> _backtracked_c_y;

  std::vector<int>   _backtracked_nloc_c_u;
  std::vector<float> _backtracked_locx_c_u;
  std::vector<float> _backtracked_locy_c_u;
  std::vector<float> _backtracked_locz_c_u;
  std::vector<float> _backtracked_dirx_c_u;
  std::vector<float> _backtracked_diry_c_u;
  std::vector<float> _backtracked_dirz_c_u;
  std::vector<float> _backtracked_de_c_u;
  std::vector<float> _backtracked_pitch_c_u;
  std::vector<float> _backtracked_scepitch_c_u;

  std::vector<int>   _backtracked_nloc_c_v;
  std::vector<float> _backtracked_locx_c_v;
  std::vector<float> _backtracked_locy_c_v;
  std::vector<float> _backtracked_locz_c_v;
  std::vector<float> _backtracked_dirx_c_v;
  std::vector<float> _backtracked_diry_c_v;
  std::vector<float> _backtracked_dirz_c_v;
  std::vector<float> _backtracked_de_c_v;
  std::vector<float> _backtracked_pitch_c_v;
  std::vector<float> _backtracked_scepitch_c_v;

  std::vector<int>   _backtracked_nloc_c_y;
  std::vector<float> _backtracked_locx_c_y;
  std::vector<float> _backtracked_locy_c_y;
  std::vector<float> _backtracked_locz_c_y;
  std::vector<float> _backtracked_dirx_c_y;
  std::vector<float> _backtracked_diry_c_y;
  std::vector<float> _backtracked_dirz_c_y;
  std::vector<float> _backtracked_de_c_y;
  std::vector<float> _backtracked_pitch_c_y;
  std::vector<float> _backtracked_scepitch_c_y;

  std::vector<float> _backtracked_x;
  std::vector<float> _backtracked_y;
  std::vector<float> _backtracked_z;

  std::vector<float> _backtracked_dirx;
  std::vector<float> _backtracked_diry;
  std::vector<float> _backtracked_dirz;

  float _backtracked_theta;
  float _backtracked_phi;

  float _backtracked_px;
  float _backtracked_py;
  float _backtracked_pz;

  float _backtracked_end_x;
  float _backtracked_end_y;
  float _backtracked_end_z;

  float _backtracked_start_x;
  float _backtracked_start_y;
  float _backtracked_start_z;
  float _backtracked_start_t;
  float _backtracked_start_U;
  float _backtracked_start_V;
  float _backtracked_start_Y;
  float _backtracked_sce_start_x;
  float _backtracked_sce_start_y;
  float _backtracked_sce_start_z;
  float _backtracked_sce_start_U;
  float _backtracked_sce_start_V;
  float _backtracked_sce_start_Y;

  float _backtracked_length;
  float _backtracked_length_start_to_end;

  bool _backtracked_process_is_stopping;
  bool _backtracked_end_in_tpc;

  uint _generation;    // generation, 1 is primary
  uint _shr_daughters; // number of shower daughters
  uint _trk_daughters; // number of track daughters
  uint _daughters;
  uint _n_pfp;

  // track information
  float _trk_score;

  float _trk_theta;
  float _trk_phi;

  float _trk_dir_x;
  float _trk_dir_y;
  float _trk_dir_z;

  float _trk_len;

  float _trk_calo_range_u;
  float _trk_calo_range_v;
  float _trk_calo_range_y;

  float _trk_start_x;
  float _trk_start_y;
  float _trk_start_z;

  float _trk_sce_start_x;
  float _trk_sce_start_y;
  float _trk_sce_start_z;

  float _trk_end_x;
  float _trk_end_y;
  float _trk_end_z;

  float _trk_sce_end_x;
  float _trk_sce_end_y;
  float _trk_sce_end_z;

  float _trk_bragg_p_y;
  float _trk_bragg_mu_y;
  float _trk_bragg_mip_y;
  float _trk_pid_chipr_y;
  float _trk_pid_chika_y;
  float _trk_pid_chipi_y;
  float _trk_pid_chimu_y;
  float _trk_pida_y;

  float _trk_bragg_p_u;
  float _trk_bragg_mu_u;
  float _trk_bragg_mip_u;
  float _trk_pid_chipr_u;
  float _trk_pid_chika_u;
  float _trk_pid_chipi_u;
  float _trk_pid_chimu_u;
  float _trk_pida_u;

  float _trk_bragg_p_v;
  float _trk_bragg_mu_v;
  float _trk_bragg_mip_v;
  float _trk_pid_chipr_v;
  float _trk_pid_chika_v;
  float _trk_pid_chipi_v;
  float _trk_pid_chimu_v;
  float _trk_pida_v;

  float _trk_bragg_p_three_planes;

  float _trk_mcs_muon_mom;
  float _trk_energy_proton;
  float _trk_energy_muon;

  int _longest; // longest track in slice?

  std::vector<float> _integrated_charge_u;
  std::vector<float> _integrated_charge_v;
  std::vector<float> _integrated_charge_y;

  std::vector<float> _allhit_charge_u;
  std::vector<float> _allhit_charge_v;
  std::vector<float> _allhit_charge_y;
  std::vector<unsigned> _allhit_wire_u;
  std::vector<unsigned> _allhit_wire_v;
  std::vector<unsigned> _allhit_wire_y;
  std::vector<unsigned> _allhit_channel_u;
  std::vector<unsigned> _allhit_channel_v;
  std::vector<unsigned> _allhit_channel_y;
  std::vector<int> _allhit_multiplicity_u;
  std::vector<int> _allhit_multiplicity_v;
  std::vector<int> _allhit_multiplicity_y;
  std::vector<float> _allhit_width_u;
  std::vector<float> _allhit_width_v;
  std::vector<float> _allhit_width_y;
  std::vector<float> _allhit_time_u;
  std::vector<float> _allhit_time_v;
  std::vector<float> _allhit_time_y;
  std::vector<unsigned> _allhit_nhit_u;
  std::vector<unsigned> _allhit_nhit_v;
  std::vector<unsigned> _allhit_nhit_y;

  std::vector<float> _trkhit_charge_u;
  std::vector<float> _trkhit_charge_v;
  std::vector<float> _trkhit_charge_y;
  std::vector<unsigned> _trkhit_wire_u;
  std::vector<unsigned> _trkhit_wire_v;
  std::vector<unsigned> _trkhit_wire_y;
  std::vector<unsigned> _trkhit_channel_u;
  std::vector<unsigned> _trkhit_channel_v;
  std::vector<unsigned> _trkhit_channel_y;
  std::vector<int> _trkhit_multiplicity_u;
  std::vector<int> _trkhit_multiplicity_v;
  std::vector<int> _trkhit_multiplicity_y;
  std::vector<float> _trkhit_width_u;
  std::vector<float> _trkhit_width_v;
  std::vector<float> _trkhit_width_y;
  std::vector<float> _trkhit_time_u;
  std::vector<float> _trkhit_time_v;
  std::vector<float> _trkhit_time_y;
  std::vector<unsigned> _trkhit_nhit_u;
  std::vector<unsigned> _trkhit_nhit_v;
  std::vector<unsigned> _trkhit_nhit_y;

  std::vector<float> _calohit_charge_u;
  std::vector<float> _calohit_charge_v;
  std::vector<float> _calohit_charge_y;
  std::vector<unsigned> _calohit_wire_u;
  std::vector<unsigned> _calohit_wire_v;
  std::vector<unsigned> _calohit_wire_y;
  std::vector<unsigned> _calohit_channel_u;
  std::vector<unsigned> _calohit_channel_v;
  std::vector<unsigned> _calohit_channel_y;
  std::vector<int> _calohit_multiplicity_u;
  std::vector<int> _calohit_multiplicity_v;
  std::vector<int> _calohit_multiplicity_y;
  std::vector<float> _calohit_width_u;
  std::vector<float> _calohit_width_v;
  std::vector<float> _calohit_width_y;
  std::vector<float> _calohit_time_u;
  std::vector<float> _calohit_time_v;
  std::vector<float> _calohit_time_y;
  std::vector<unsigned> _calohit_nhit_u;
  std::vector<unsigned> _calohit_nhit_v;
  std::vector<unsigned> _calohit_nhit_y;

  std::vector<float> _sumhit_charge_u;
  std::vector<float> _sumhit_charge_v;
  std::vector<float> _sumhit_charge_y;
  std::vector<unsigned> _sumhit_wire_u;
  std::vector<unsigned> _sumhit_wire_v;
  std::vector<unsigned> _sumhit_wire_y;
  std::vector<unsigned> _sumhit_channel_u;
  std::vector<unsigned> _sumhit_channel_v;
  std::vector<unsigned> _sumhit_channel_y;
  std::vector<int> _sumhit_multiplicity_u;
  std::vector<int> _sumhit_multiplicity_v;
  std::vector<int> _sumhit_multiplicity_y;
  std::vector<float> _sumhit_width_u;
  std::vector<float> _sumhit_width_v;
  std::vector<float> _sumhit_width_y;
  std::vector<float> _sumhit_time_u;
  std::vector<float> _sumhit_time_v;
  std::vector<float> _sumhit_time_y;
  std::vector<unsigned> _sumhit_nhit_u;
  std::vector<unsigned> _sumhit_nhit_v;
  std::vector<unsigned> _sumhit_nhit_y;

  std::vector<float> _areahit_charge_u;
  std::vector<float> _areahit_charge_v;
  std::vector<float> _areahit_charge_y;
  std::vector<unsigned> _areahit_wire_u;
  std::vector<unsigned> _areahit_wire_v;
  std::vector<unsigned> _areahit_wire_y;
  std::vector<unsigned> _areahit_channel_u;
  std::vector<unsigned> _areahit_channel_v;
  std::vector<unsigned> _areahit_channel_y;
  std::vector<int> _areahit_multiplicity_u;
  std::vector<int> _areahit_multiplicity_v;
  std::vector<int> _areahit_multiplicity_y;
  std::vector<float> _areahit_width_u;
  std::vector<float> _areahit_width_v;
  std::vector<float> _areahit_width_y;
  std::vector<float> _areahit_time_u;
  std::vector<float> _areahit_time_v;
  std::vector<float> _areahit_time_y;
  std::vector<unsigned> _areahit_nhit_u;
  std::vector<unsigned> _areahit_nhit_v;
  std::vector<unsigned> _areahit_nhit_y;

  // dedx vector
  std::vector<float> _dqdx_u;
  std::vector<float> _dqdx_v;
  std::vector<float> _dqdx_y;

  std::vector<float> _dedx_u;
  std::vector<float> _dedx_v;
  std::vector<float> _dedx_y;

  std::vector<unsigned> _dedx_channel_u;
  std::vector<unsigned> _dedx_channel_v;
  std::vector<unsigned> _dedx_channel_y;

  std::vector<float> _rr_u;
  std::vector<float> _rr_v;
  std::vector<float> _rr_y;

  std::vector<float> _pitch_u;
  std::vector<float> _pitch_v;
  std::vector<float> _pitch_y;

  std::vector<float> _sx;
  std::vector<float> _sy;
  std::vector<float> _sz;

  std::vector<float> _w_u;
  std::vector<float> _w_v;
  std::vector<float> _w_y;

  std::vector<float> _t_u;
  std::vector<float> _t_v;
  std::vector<float> _t_y;

  std::vector<float> _x_u;
  std::vector<float> _x_v;
  std::vector<float> _x_y;

  std::vector<float> _y_u;
  std::vector<float> _y_v;
  std::vector<float> _y_y;

  std::vector<float> _z_u;
  std::vector<float> _z_v;
  std::vector<float> _z_y;

  std::vector<float> _dir_x_u;
  std::vector<float> _dir_x_v;
  std::vector<float> _dir_x_y;

  std::vector<float> _dir_y_u;
  std::vector<float> _dir_y_v;
  std::vector<float> _dir_y_y;

  std::vector<float> _dir_z_u;
  std::vector<float> _dir_z_v;
  std::vector<float> _dir_z_y;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CalorimetryAnalysis::CalorimetryAnalysis(const fhicl::ParameterSet &p)
  : art::EDAnalyzer{p}

{
  fPFPproducer  = p.get< art::InputTag > ("PFPproducer","pandoraTrackGausCryo0");
  fCALOproducer = p.get< art::InputTag > ("CALOproducer");
  fPIDproducer  = p.get< art::InputTag > ("PIDproducer" );
  fTRKproducer  = p.get< art::InputTag > ("TRKproducer" );
  fMCSproducer = p.get<std::string>("MCSproducer", "pandoraTrackMCS");
  fRangeproducer = p.get<std::string>("Rangeproducer", "pandoraTrackRange");
  fT0producer  = p.get< art::InputTag > ("T0producer", "" );
  fAreaHitproducer = p.get<art::InputTag> ("AreaHitproducer", "areahit");
  fAllTrueEnergyDeposits = p.get<bool>("AllTrueEnergyDeposits", true);
  fHitFilterproducer = p.get<art::InputTag>("HitFilterproducer", "filtgoodhit"); 
  fHitProducers = p.get<std::vector<art::InputTag>>("HitProducers");
  fWireProducer = p.get<art::InputTag>("WireProducer", "decon1DroiTPC0");

  art::ServiceHandle<art::TFileService> tfs;

  _calo_tree = tfs->make<TTree>("CalorimetryAnalyzer", "Calo Tree");

  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // first the TPC volumes 
  for (auto const &cryo: geometry->IterateCryostats()) {
    geo::GeometryCore::TPC_iterator iTPC = geometry->begin_TPC(cryo.ID()),
                                    tend = geometry->end_TPC(cryo.ID());
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    while (iTPC != tend) {
      geo::TPCGeo const& TPC = *iTPC;
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
      iTPC++;
    }
     fTPCVolumes.push_back(std::move(this_tpc_volumes));
  }

  // then combine them into active volumes
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: fTPCVolumes) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
  }

  _fileno = 0;

  setBranches(NULL);
}

void CalorimetryAnalysis::analyze(art::Event const &e)
{
  bool fData = e.isRealData();
  (void) fData;

  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  std::cout << "[CalorimetryAnalysis::analyzeEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: "<< _evt << std::endl;

  // Build larpandora info:
  //lar_pandora::LArPandoraHelper larpandora;
  //lar_pandora::PFParticleVector pfparticles;
  //lar_pandora::PFParticleMap particleMap;
  //larpandora.CollectPFParticles(e, "pandora", pfparticles);
  //larpandora.BuildPFParticleMap(pfparticles, particleMap);

  // recob Hits
  std::vector<art::Ptr<recob::Hit>> hitList;
  for (const art::InputTag &t: fHitProducers) {
    art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(t);
    art::fill_ptr_vector(hitList, hitHandle);
  }
  // recob Wires:w

  art::Handle<std::vector<recob::Wire>> wireHandle;
  e.getByLabel(fWireProducer, wireHandle);

  std::vector<art::Ptr<recob::Wire>> wireList;
  if (wireHandle.isValid()) {
    art::fill_ptr_vector(wireList, wireHandle);
  }

  // true particles
  art::ValidHandle<std::vector<simb::MCParticle>> particleHandle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  std::vector<art::Ptr<simb::MCParticle>> particleList;
  art::fill_ptr_vector(particleList, particleHandle);

  art::ValidHandle<std::vector<sim::SimChannel>> simchannelHandle = e.getValidHandle<std::vector<sim::SimChannel>>("largeant");
  std::vector<art::Ptr<sim::SimChannel>> simchannelList;
  art::fill_ptr_vector(simchannelList, simchannelHandle);

  art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
  std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
  art::fill_ptr_vector(PFParticleList, pfparticles);
  _n_pfp = pfparticles->size();

  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(fTRKproducer); 

  art::FindManyP<recob::SpacePoint> fmSpacePoint(PFParticleList, e, fPFPproducer);

  art::FindManyP<recob::Track> fmTracks(pfparticles, e, fTRKproducer);
  // art::FindManyP<recob::SpacePoint> fmSpacepoints(pfparticles, e, fPFPproducer); 
  // also get hits...
  art::FindManyP<recob::Hit> fmtrkHits(tracks, e, fTRKproducer);
  art::FindManyP<recob::Hit> fmareaHits(tracks, e, fAreaHitproducer);
  art::FindManyP<recob::Hit> fmcaloHits(tracks, e, fHitFilterproducer);

  art::FindManyP<anab::Calorimetry> fmCalo(tracks, e, fCALOproducer);
  art::FindManyP<anab::ParticleID> fmPID(tracks, e, fPIDproducer);

  art::InputTag mcs_muon  {fMCSproducer, "muon"};
  art::InputTag range_muon {fRangeproducer, "muon"};
  art::InputTag range_proton {fRangeproducer, "proton"};

  art::FindManyP<recob::MCSFitResult> fmMuonMCS(tracks, e, mcs_muon);
  art::FindManyP<sbn::RangeP> fmMuonRange(tracks, e, range_muon);
  art::FindManyP<sbn::RangeP> fmProtonRange(tracks, e, range_proton);

  // service data
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  for (art::Ptr<recob::PFParticle> p_pfp: PFParticleList) {
    const recob::PFParticle &pfp = *p_pfp;

    fillDefault();

    // grab associated tracks
    const std::vector<art::Ptr<recob::Track>> thisTrack = fmTracks.at(p_pfp.key());
    if (thisTrack.size() != 1)
      continue;

    size_t parent_ind = pfp.Parent();
    std::cout << "Parent ind: " << parent_ind << std::endl;
    if (parent_ind == recob::PFParticle::kPFParticlePrimary || parent_ind >= PFParticleList.size()) {
      std::cout << "Parent doesn't exist.\n";
      continue;
    }
    // check if parent is the primary
    const recob::PFParticle *p_parent = NULL;
    for (unsigned i_p = 0; i_p < PFParticleList.size(); i_p++) {
      if (PFParticleList[i_p]->Self() == pfp.Parent()) {
        p_parent = &*PFParticleList[i_p];
        break;
      }
    }

    if (p_parent == NULL || !p_parent->IsPrimary()) {
      std::cout << "Parent not primary\n";
      continue;
    }

    art::Ptr<recob::Track> trkPtr = thisTrack.at(0);
    const recob::Track &track = *trkPtr;

    // get all the data

    // track data
    std::vector<art::Ptr<recob::SpacePoint>> emptySPVector;
    const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints = fmSpacePoint.isValid() ? fmSpacePoint.at(p_pfp.key()) : emptySPVector;

    std::vector<art::Ptr<anab::Calorimetry>> emptyCaloVector;
    const std::vector<art::Ptr<anab::Calorimetry>> &calo = fmCalo.isValid() ? fmCalo.at(trkPtr.key()) : emptyCaloVector;

    std::vector<art::Ptr<anab::ParticleID>> emptyPIDVector;
    const std::vector<art::Ptr<anab::ParticleID>> &pid = fmPID.isValid() ? fmPID.at(trkPtr.key()) : emptyPIDVector;

    recob::MCSFitResult muon_mcs;
    if (fmMuonMCS.isValid() && fmMuonMCS.at(trkPtr.key()).size()) muon_mcs = *fmMuonMCS.at(trkPtr.key()).at(0);

    sbn::RangeP muon_range;
    if (fmMuonRange.isValid() && fmMuonRange.at(trkPtr.key()).size()) muon_range = *fmMuonRange.at(trkPtr.key()).at(0);

    sbn::RangeP proton_range; 
    if (fmProtonRange.isValid() && fmProtonRange.at(trkPtr.key()).size()) proton_range = *fmProtonRange.at(trkPtr.key()).at(0);

    std::vector<art::Ptr<recob::Hit>> emptyHitVector;
    const std::vector<art::Ptr<recob::Hit>> &trkHits  = fmtrkHits.isValid() ? fmtrkHits.at(trkPtr.key()) : emptyHitVector;
    const std::vector<art::Ptr<recob::Hit>> &areaHits = fmareaHits.isValid() ?  fmareaHits.at(trkPtr.key()) : emptyHitVector;
    const std::vector<art::Ptr<recob::Hit>> &caloHits = fmcaloHits.isValid() ? fmcaloHits.at(trkPtr.key()) : emptyHitVector;
    
    // Get the true matching MC particle
    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, trkHits, true);
    // get the match with the most energy
    int match_id = matches.size() ? std::max_element(matches.begin(), matches.end(),
                                                     [](const auto &a, const auto &b) { return a.second < b.second; })->first : -1;
    float matchE = matches.size() ? std::max_element(matches.begin(), matches.end(),
                                                     [](const auto &a, const auto &b) { return a.second < b.second; })->second : -1;
    art::Ptr<simb::MCParticle> true_particle;
    for (art::Ptr<simb::MCParticle> p: particleList) {
      if (p->TrackId() == match_id) {
        true_particle = p;
        break;
      }
    }
    float trueTrackE = 0.;
    for (auto pair: matches) trueTrackE += pair.second;

    // TODO: fix??
    //
    // ignore tracks with no match
    if (!true_particle) continue;
    // only take particle that matches to true instigator
    if (true_particle->Process() != "primary") continue;


    std::cout << "Track: " << track.ID() << " len: " << track.Length() << " matched to particle: " << true_particle->TrackId() << " frac: " << (matchE/trueTrackE) << std::endl;

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    std::array<std::map<unsigned, std::vector<const sim::IDE*>>, 3> particle_ide_map;
    for (const art::Ptr<sim::SimChannel> sc: simchannelList) {
      for (const auto &pair: sc->TDCIDEMap()) {
        for (const sim::IDE &ide: pair.second) {
          if (ide.trackID == true_particle->TrackId() || fAllTrueEnergyDeposits) {
            unsigned plane_id = geometry->ChannelToWire(sc->Channel()).at(0).Plane;
            particle_ide_map[plane_id][sc->Channel()].push_back(&ide);
          }
        }
      }
    } 

    std::array<std::vector<const sim::IDE*>, 3> allIDEs;
    for (unsigned i = 0; i < particleList.size(); i++) {
      int G4ID = particleList[i]->TrackId();
      for (unsigned plane = 0; plane < 3; plane++) {
        // build the plane ID
        geo::PlaneID plane_id(0, 0, plane);
        std::vector<const sim::IDE*> this_ides = bt_serv->TrackIdToSimIDEs_Ps(G4ID, geometry->View(plane_id));
        allIDEs[plane].insert(allIDEs[plane].end(), this_ides.begin(), this_ides.end());
      }
    }

    FillCalorimetry(e,
          wireList,
          hitList,
          p_pfp,
	  thisTrack.at(0),
          spacepoints,
          calo,
          pid,
          areaHits,
          trkHits,
          caloHits,
          muon_mcs,
          muon_range,
          proton_range,
          false,
          matches,
          true_particle,
          particle_ide_map, 
          allIDEs,
          fActiveVolumes);

  }// for all PFParticles
}// analyzeEvent

void CalorimetryAnalysis::fillDefault()
{
  // backtracking information
  _backtracked_pdg = std::numeric_limits<int>::lowest();            // PDG code of backtracked particle
  _backtracked_e = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_deposited_e_u = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_deposited_e_v = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_deposited_e_y = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_q_u = std::numeric_limits<float>::lowest();
  _backtracked_q_v = std::numeric_limits<float>::lowest();
  _backtracked_q_y = std::numeric_limits<float>::lowest();
  _backtracked_qtotal_u = std::numeric_limits<float>::lowest();
  _backtracked_qtotal_v = std::numeric_limits<float>::lowest();
  _backtracked_qtotal_y = std::numeric_limits<float>::lowest();
  _backtracked_purity = std::numeric_limits<float>::lowest();       // purity of backtracking
  _backtracked_completeness = std::numeric_limits<float>::lowest(); // completeness of backtracking
  _backtracked_overlay_purity = std::numeric_limits<float>::lowest(); // purity of overlay
  _backtracked_process_is_stopping = false;
  _backtracked_end_in_tpc = false;

  _backtracked_pc_q_u.clear();
  _backtracked_pc_q_v.clear();
  _backtracked_pc_q_y.clear();

  _backtracked_pc_e_u.clear();
  _backtracked_pc_e_v.clear();
  _backtracked_pc_e_y.clear();

  _backtracked_c_u.clear();
  _backtracked_c_v.clear();
  _backtracked_c_y.clear();

  _backtracked_nloc_c_u.clear();
  _backtracked_locx_c_u.clear();
  _backtracked_locy_c_u.clear();
  _backtracked_locz_c_u.clear();
  _backtracked_dirx_c_u.clear();
  _backtracked_diry_c_u.clear();
  _backtracked_dirz_c_u.clear();
  _backtracked_de_c_u.clear();
  _backtracked_pitch_c_u.clear();
  _backtracked_scepitch_c_u.clear();

  _backtracked_nloc_c_v.clear();
  _backtracked_locx_c_v.clear();
  _backtracked_locy_c_v.clear();
  _backtracked_locz_c_v.clear();
  _backtracked_dirx_c_v.clear();
  _backtracked_diry_c_v.clear();
  _backtracked_dirz_c_v.clear();
  _backtracked_de_c_v.clear();
  _backtracked_pitch_c_v.clear();
  _backtracked_scepitch_c_v.clear();

  _backtracked_nloc_c_y.clear();
  _backtracked_locx_c_y.clear();
  _backtracked_locy_c_y.clear();
  _backtracked_locz_c_y.clear();
  _backtracked_dirx_c_y.clear();
  _backtracked_diry_c_y.clear();
  _backtracked_dirz_c_y.clear();
  _backtracked_de_c_y.clear();
  _backtracked_pitch_c_y.clear();
  _backtracked_scepitch_c_y.clear();

  _backtracked_x.clear();
  _backtracked_y.clear();
  _backtracked_z.clear();

  _backtracked_dirx.clear();
  _backtracked_diry.clear();
  _backtracked_dirz.clear();

  _backtracked_length = std::numeric_limits<float>::lowest();
  _backtracked_length_start_to_end = std::numeric_limits<float>::lowest();

  _backtracked_theta = std::numeric_limits<float>::lowest();
  _backtracked_phi = std::numeric_limits<float>::lowest();

  _backtracked_px = std::numeric_limits<float>::lowest();
  _backtracked_py = std::numeric_limits<float>::lowest();
  _backtracked_pz = std::numeric_limits<float>::lowest();

  _backtracked_end_x = std::numeric_limits<float>::lowest();
  _backtracked_end_y = std::numeric_limits<float>::lowest();
  _backtracked_end_z = std::numeric_limits<float>::lowest();

  _backtracked_start_x = std::numeric_limits<float>::lowest();
  _backtracked_start_y = std::numeric_limits<float>::lowest();
  _backtracked_start_z = std::numeric_limits<float>::lowest();
  _backtracked_start_t = std::numeric_limits<float>::lowest();
  _backtracked_start_U = std::numeric_limits<float>::lowest();
  _backtracked_start_V = std::numeric_limits<float>::lowest();
  _backtracked_start_Y = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_x = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_y = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_z = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_U = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_V = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_Y = std::numeric_limits<float>::lowest();

  _generation = std::numeric_limits<uint>::lowest();
  _shr_daughters = std::numeric_limits<uint>::lowest();
  _trk_daughters = std::numeric_limits<uint>::lowest();
  _daughters = std::numeric_limits<uint>::lowest();

  // track information
  _trk_score = std::numeric_limits<float>::lowest();

  _trk_theta = std::numeric_limits<float>::lowest();
  _trk_phi = std::numeric_limits<float>::lowest();
  _trk_len = std::numeric_limits<float>::lowest();

  _trk_calo_range_u = std::numeric_limits<float>::lowest();
  _trk_calo_range_v = std::numeric_limits<float>::lowest();
  _trk_calo_range_y = std::numeric_limits<float>::lowest();

  _trk_dir_x = std::numeric_limits<float>::lowest();
  _trk_dir_y = std::numeric_limits<float>::lowest();
  _trk_dir_z = std::numeric_limits<float>::lowest();

  _trk_start_x = std::numeric_limits<float>::lowest();
  _trk_start_y = std::numeric_limits<float>::lowest();
  _trk_start_z = std::numeric_limits<float>::lowest();

  _trk_sce_start_x = std::numeric_limits<float>::lowest();
  _trk_sce_start_y = std::numeric_limits<float>::lowest();
  _trk_sce_start_z = std::numeric_limits<float>::lowest();

  _trk_end_x = std::numeric_limits<float>::lowest();
  _trk_end_y = std::numeric_limits<float>::lowest();
  _trk_end_z = std::numeric_limits<float>::lowest();

  _trk_sce_end_x = std::numeric_limits<float>::lowest();
  _trk_sce_end_y = std::numeric_limits<float>::lowest();
  _trk_sce_end_z = std::numeric_limits<float>::lowest();

  _trk_bragg_p_y = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_y = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_y = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_y = std::numeric_limits<float>::lowest();
  _trk_pid_chika_y = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_y = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_y = std::numeric_limits<float>::lowest();
  _trk_pida_y = std::numeric_limits<float>::lowest();

  _trk_bragg_p_u = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_u = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_u = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_u = std::numeric_limits<float>::lowest();
  _trk_pid_chika_u = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_u = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_u = std::numeric_limits<float>::lowest();
  _trk_pida_u = std::numeric_limits<float>::lowest();

  _trk_bragg_p_v = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_v = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_v = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_v = std::numeric_limits<float>::lowest();
  _trk_pid_chika_v = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_v = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_v = std::numeric_limits<float>::lowest();
  _trk_pida_v = std::numeric_limits<float>::lowest();

  _trk_bragg_p_three_planes = std::numeric_limits<float>::lowest();

  _longest = std::numeric_limits<int>::lowest();

  _trk_mcs_muon_mom = std::numeric_limits<float>::lowest();
  _trk_energy_proton = std::numeric_limits<float>::lowest();
  _trk_energy_muon = std::numeric_limits<float>::lowest();


  _integrated_charge_u.clear();
  _integrated_charge_v.clear();
  _integrated_charge_y.clear();

  _trkhit_charge_u.clear();
  _trkhit_charge_v.clear();
  _trkhit_charge_y.clear();
  _trkhit_wire_u.clear();
  _trkhit_wire_v.clear();
  _trkhit_wire_y.clear();
  _trkhit_channel_u.clear();
  _trkhit_channel_v.clear();
  _trkhit_channel_y.clear();
  _trkhit_multiplicity_u.clear();
  _trkhit_multiplicity_v.clear();
  _trkhit_multiplicity_y.clear();
  _trkhit_width_u.clear();
  _trkhit_width_v.clear();
  _trkhit_width_y.clear();
  _trkhit_time_u.clear();
  _trkhit_time_v.clear();
  _trkhit_time_y.clear();
  _trkhit_nhit_u.clear();
  _trkhit_nhit_v.clear();
  _trkhit_nhit_y.clear();

  _areahit_charge_u.clear();
  _areahit_charge_v.clear();
  _areahit_charge_y.clear();
  _areahit_wire_u.clear();
  _areahit_wire_v.clear();
  _areahit_wire_y.clear();
  _areahit_channel_u.clear();
  _areahit_channel_v.clear();
  _areahit_channel_y.clear();
  _areahit_multiplicity_u.clear();
  _areahit_multiplicity_v.clear();
  _areahit_multiplicity_y.clear();
  _areahit_width_u.clear();
  _areahit_width_v.clear();
  _areahit_width_y.clear();
  _areahit_time_u.clear();
  _areahit_time_v.clear();
  _areahit_time_y.clear();
  _areahit_nhit_u.clear();
  _areahit_nhit_v.clear();
  _areahit_nhit_y.clear();

  _allhit_charge_u.clear();
  _allhit_charge_v.clear();
  _allhit_charge_y.clear();
  _allhit_wire_u.clear();
  _allhit_wire_v.clear();
  _allhit_wire_y.clear();
  _allhit_channel_u.clear();
  _allhit_channel_v.clear();
  _allhit_channel_y.clear();
  _allhit_multiplicity_u.clear();
  _allhit_multiplicity_v.clear();
  _allhit_multiplicity_y.clear();
  _allhit_width_u.clear();
  _allhit_width_v.clear();
  _allhit_width_y.clear();
  _allhit_time_u.clear();
  _allhit_time_v.clear();
  _allhit_time_y.clear();
  _allhit_nhit_u.clear();
  _allhit_nhit_v.clear();
  _allhit_nhit_y.clear();

  _calohit_charge_u.clear();
  _calohit_charge_v.clear();
  _calohit_charge_y.clear();
  _calohit_wire_u.clear();
  _calohit_wire_v.clear();
  _calohit_wire_y.clear();
  _calohit_channel_u.clear();
  _calohit_channel_v.clear();
  _calohit_channel_y.clear();
  _calohit_multiplicity_u.clear();
  _calohit_multiplicity_v.clear();
  _calohit_multiplicity_y.clear();
  _calohit_width_u.clear();
  _calohit_width_v.clear();
  _calohit_width_y.clear();
  _calohit_time_u.clear();
  _calohit_time_v.clear();
  _calohit_time_y.clear();
  _calohit_nhit_u.clear();
  _calohit_nhit_v.clear();
  _calohit_nhit_y.clear();

  _sumhit_charge_u.clear();
  _sumhit_charge_v.clear();
  _sumhit_charge_y.clear();
  _sumhit_wire_u.clear();
  _sumhit_wire_v.clear();
  _sumhit_wire_y.clear();
  _sumhit_channel_u.clear();
  _sumhit_channel_v.clear();
  _sumhit_channel_y.clear();
  _sumhit_multiplicity_u.clear();
  _sumhit_multiplicity_v.clear();
  _sumhit_multiplicity_y.clear();
  _sumhit_width_u.clear();
  _sumhit_width_v.clear();
  _sumhit_width_y.clear();
  _sumhit_time_u.clear();
  _sumhit_time_v.clear();
  _sumhit_time_y.clear();
  _sumhit_nhit_u.clear();
  _sumhit_nhit_v.clear();
  _sumhit_nhit_y.clear();

  // dedx vector
  _dqdx_u.clear();
  _dqdx_v.clear();
  _dqdx_y.clear();

  _dedx_u.clear();
  _dedx_v.clear();
  _dedx_y.clear();

  _dedx_channel_u.clear();
  _dedx_channel_v.clear();
  _dedx_channel_y.clear();

  _rr_u.clear();
  _rr_v.clear();
  _rr_y.clear();

  _pitch_u.clear();
  _pitch_v.clear();
  _pitch_y.clear();

  _sx.clear();
  _sy.clear();
  _sz.clear();

  _w_u.clear();
  _w_v.clear();
  _w_y.clear();

  _t_u.clear();
  _t_v.clear();
  _t_y.clear();

  _x_u.clear();
  _x_v.clear();
  _x_y.clear();

  _y_u.clear();
  _y_v.clear();
  _y_y.clear();

  _z_u.clear();
  _z_v.clear();
  _z_y.clear();

  _dir_x_u.clear();
  _dir_x_v.clear();
  _dir_x_y.clear();

  _dir_y_u.clear();
  _dir_y_v.clear();
  _dir_y_y.clear();

  _dir_z_u.clear();
  _dir_z_v.clear();
  _dir_z_y.clear();
}

void CalorimetryAnalysis::setBranches(TTree *_tree)
{
  _calo_tree->Branch("fileno", &_fileno, "fileno/i");
  _calo_tree->Branch("run", &_run, "run/i");
  _calo_tree->Branch("sub", &_sub, "sub/i");
  _calo_tree->Branch("evt", &_evt, "evt/i");
  _calo_tree->Branch("pfpno", &_pfpno, "pfpno/i");
  _calo_tree->Branch("trkid", &_trkid, "trkid/i");

  // backtracking information
  _calo_tree->Branch("backtracked_pdg", &_backtracked_pdg, "backtracked_pdg/I");            // PDG code of backtracked particle
  _calo_tree->Branch("backtracked_e", &_backtracked_e, "backtracked_e/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_deposited_e_u", &_backtracked_deposited_e_u, "backtracked_deposited_e_u/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_deposited_e_v", &_backtracked_deposited_e_v, "backtracked_deposited_e_v/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_deposited_e_y", &_backtracked_deposited_e_y, "backtracked_deposited_e_y/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_q_u", &_backtracked_q_u, "backtracked_q_u/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_q_v", &_backtracked_q_v, "backtracked_q_v/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_q_y", &_backtracked_q_y, "backtracked_q_y/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_qtotal_u", &_backtracked_qtotal_u, "backtracked_qtotal_u/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_qtotal_v", &_backtracked_qtotal_v, "backtracked_qtotal_v/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_qtotal_y", &_backtracked_qtotal_y, "backtracked_qtotal_y/F");            // energy of backtracked particle

  _calo_tree->Branch("backtracked_pc_q_u", "std::vector<float>", &_backtracked_pc_q_u);
  _calo_tree->Branch("backtracked_pc_q_v", "std::vector<float>", &_backtracked_pc_q_v);
  _calo_tree->Branch("backtracked_pc_q_y", "std::vector<float>", &_backtracked_pc_q_y);

  _calo_tree->Branch("backtracked_pc_e_u", "std::vector<float>", &_backtracked_pc_e_u);
  _calo_tree->Branch("backtracked_pc_e_v", "std::vector<float>", &_backtracked_pc_e_v);
  _calo_tree->Branch("backtracked_pc_e_y", "std::vector<float>", &_backtracked_pc_e_y);

  _calo_tree->Branch("backtracked_c_u", "std::vector<unsigned>", &_backtracked_c_u);
  _calo_tree->Branch("backtracked_c_v", "std::vector<unsigned>", &_backtracked_c_v);
  _calo_tree->Branch("backtracked_c_y", "std::vector<unsigned>", &_backtracked_c_y);

  _calo_tree->Branch("backtracked_nloc_c_u", "std::vector<int>", &_backtracked_nloc_c_u);
  _calo_tree->Branch("backtracked_locx_c_u", "std::vector<float>", &_backtracked_locx_c_u);
  _calo_tree->Branch("backtracked_locy_c_u", "std::vector<float>", &_backtracked_locy_c_u);
  _calo_tree->Branch("backtracked_locz_c_u", "std::vector<float>", &_backtracked_locz_c_u);
  _calo_tree->Branch("backtracked_dirx_c_u", "std::vector<float>", &_backtracked_dirx_c_u);
  _calo_tree->Branch("backtracked_diry_c_u", "std::vector<float>", &_backtracked_diry_c_u);
  _calo_tree->Branch("backtracked_dirz_c_u", "std::vector<float>", &_backtracked_dirz_c_u);
  _calo_tree->Branch("backtracked_de_c_u", "std::vector<float>", &_backtracked_de_c_u);
  _calo_tree->Branch("backtracked_pitch_c_u", "std::vector<float>", &_backtracked_pitch_c_u);
  _calo_tree->Branch("backtracked_scepitch_c_u", "std::vector<float>", &_backtracked_scepitch_c_u);

  _calo_tree->Branch("backtracked_nloc_c_v", "std::vector<int>", &_backtracked_nloc_c_v);
  _calo_tree->Branch("backtracked_locx_c_v", "std::vector<float>", &_backtracked_locx_c_v);
  _calo_tree->Branch("backtracked_locy_c_v", "std::vector<float>", &_backtracked_locy_c_v);
  _calo_tree->Branch("backtracked_locz_c_v", "std::vector<float>", &_backtracked_locz_c_v);
  _calo_tree->Branch("backtracked_dirx_c_v", "std::vector<float>", &_backtracked_dirx_c_v);
  _calo_tree->Branch("backtracked_diry_c_v", "std::vector<float>", &_backtracked_diry_c_v);
  _calo_tree->Branch("backtracked_dirz_c_v", "std::vector<float>", &_backtracked_dirz_c_v);
  _calo_tree->Branch("backtracked_de_c_v", "std::vector<float>", &_backtracked_de_c_v);
  _calo_tree->Branch("backtracked_pitch_c_v", "std::vector<float>", &_backtracked_pitch_c_v);
  _calo_tree->Branch("backtracked_scepitch_c_v", "std::vector<float>", &_backtracked_scepitch_c_v);

  _calo_tree->Branch("backtracked_nloc_c_y", "std::vector<int>", &_backtracked_nloc_c_y);
  _calo_tree->Branch("backtracked_locx_c_y", "std::vector<float>", &_backtracked_locx_c_y);
  _calo_tree->Branch("backtracked_locy_c_y", "std::vector<float>", &_backtracked_locy_c_y);
  _calo_tree->Branch("backtracked_locz_c_y", "std::vector<float>", &_backtracked_locz_c_y);
  _calo_tree->Branch("backtracked_dirx_c_y", "std::vector<float>", &_backtracked_dirx_c_y);
  _calo_tree->Branch("backtracked_diry_c_y", "std::vector<float>", &_backtracked_diry_c_y);
  _calo_tree->Branch("backtracked_dirz_c_y", "std::vector<float>", &_backtracked_dirz_c_y);
  _calo_tree->Branch("backtracked_de_c_y", "std::vector<float>", &_backtracked_de_c_y);
  _calo_tree->Branch("backtracked_pitch_c_y", "std::vector<float>", &_backtracked_pitch_c_y);
  _calo_tree->Branch("backtracked_scepitch_c_y", "std::vector<float>", &_backtracked_scepitch_c_y);

  _calo_tree->Branch("backtracked_x", "std::vector<float>", &_backtracked_x);
  _calo_tree->Branch("backtracked_y", "std::vector<float>", &_backtracked_y);
  _calo_tree->Branch("backtracked_z", "std::vector<float>", &_backtracked_z);

  _calo_tree->Branch("backtracked_dirx", "std::vector<float>", &_backtracked_dirx);
  _calo_tree->Branch("backtracked_diry", "std::vector<float>", &_backtracked_diry);
  _calo_tree->Branch("backtracked_dirz", "std::vector<float>", &_backtracked_dirz);

  _calo_tree->Branch("backtracked_theta", &_backtracked_theta, "backtracked_theta/F");
  _calo_tree->Branch("backtracked_phi", &_backtracked_phi, "backtracked_phi/F");

  _calo_tree->Branch("backtracked_px", &_backtracked_px, "backtracked_px/F");
  _calo_tree->Branch("backtracked_py", &_backtracked_py, "backtracked_py/F");
  _calo_tree->Branch("backtracked_pz", &_backtracked_pz, "backtracked_pz/F");

  _calo_tree->Branch("backtracked_length", &_backtracked_length,"backtracked_length/F");
  _calo_tree->Branch("backtracked_length_start_to_end", &_backtracked_length_start_to_end, "backtracked_length_start_to_end/F");

  _calo_tree->Branch("backtracked_purity", &_backtracked_purity, "backtracked_purity/F");       // purity of backtracking
  _calo_tree->Branch("backtracked_completeness", &_backtracked_completeness, "backtracked_completeness/F"); // completeness of backtracking
  _calo_tree->Branch("backtracked_overlay_purity", &_backtracked_overlay_purity, "backtracked_overlay_purity/F"); // purity of overlay

  _calo_tree->Branch("integrated_charge_u", "std::vector<float>", &_integrated_charge_u);
  _calo_tree->Branch("integrated_charge_v", "std::vector<float>", &_integrated_charge_v);
  _calo_tree->Branch("integrated_charge_y", "std::vector<float>", &_integrated_charge_y);

  _calo_tree->Branch("areahit_charge_u", "std::vector<float>", &_areahit_charge_u);
  _calo_tree->Branch("areahit_charge_v", "std::vector<float>", &_areahit_charge_v);
  _calo_tree->Branch("areahit_charge_y", "std::vector<float>", &_areahit_charge_y);
  _calo_tree->Branch("areahit_wire_u", "std::vector<unsigned>", &_areahit_wire_u);
  _calo_tree->Branch("areahit_wire_v", "std::vector<unsigned>", &_areahit_wire_v);
  _calo_tree->Branch("areahit_wire_y", "std::vector<unsigned>", &_areahit_wire_y);
  _calo_tree->Branch("areahit_channel_u", "std::vector<unsigned>", &_areahit_channel_u);
  _calo_tree->Branch("areahit_channel_v", "std::vector<unsigned>", &_areahit_channel_v);
  _calo_tree->Branch("areahit_channel_y", "std::vector<unsigned>", &_areahit_channel_y);
  _calo_tree->Branch("areahit_multiplicity_u", "std::vector<int>", &_areahit_multiplicity_u); 
  _calo_tree->Branch("areahit_multiplicity_v", "std::vector<int>", &_areahit_multiplicity_v); 
  _calo_tree->Branch("areahit_multiplicity_y", "std::vector<int>", &_areahit_multiplicity_y); 
  _calo_tree->Branch("areahit_width_u", "std::vector<float>", &_areahit_width_u);
  _calo_tree->Branch("areahit_width_v", "std::vector<float>", &_areahit_width_v);
  _calo_tree->Branch("areahit_width_y", "std::vector<float>", &_areahit_width_y);
  _calo_tree->Branch("areahit_time_u", "std::vector<float>", &_areahit_time_u);
  _calo_tree->Branch("areahit_time_v", "std::vector<float>", &_areahit_time_v);
  _calo_tree->Branch("areahit_time_y", "std::vector<float>", &_areahit_time_y);
  _calo_tree->Branch("areahit_nhit_u", "std::vector<unsigned>", &_areahit_nhit_u);
  _calo_tree->Branch("areahit_nhit_v", "std::vector<unsigned>", &_areahit_nhit_v);
  _calo_tree->Branch("areahit_nhit_y", "std::vector<unsigned>", &_areahit_nhit_y);

  _calo_tree->Branch("allhit_charge_u", "std::vector<float>", &_allhit_charge_u);
  _calo_tree->Branch("allhit_charge_v", "std::vector<float>", &_allhit_charge_v);
  _calo_tree->Branch("allhit_charge_y", "std::vector<float>", &_allhit_charge_y);
  _calo_tree->Branch("allhit_wire_u", "std::vector<unsigned>", &_allhit_wire_u);
  _calo_tree->Branch("allhit_wire_v", "std::vector<unsigned>", &_allhit_wire_v);
  _calo_tree->Branch("allhit_wire_y", "std::vector<unsigned>", &_allhit_wire_y);
  _calo_tree->Branch("allhit_channel_u", "std::vector<unsigned>", &_allhit_channel_u);
  _calo_tree->Branch("allhit_channel_v", "std::vector<unsigned>", &_allhit_channel_v);
  _calo_tree->Branch("allhit_channel_y", "std::vector<unsigned>", &_allhit_channel_y);
  _calo_tree->Branch("allhit_multiplicity_u", "std::vector<int>", &_allhit_multiplicity_u); 
  _calo_tree->Branch("allhit_multiplicity_v", "std::vector<int>", &_allhit_multiplicity_v); 
  _calo_tree->Branch("allhit_multiplicity_y", "std::vector<int>", &_allhit_multiplicity_y); 
  _calo_tree->Branch("allhit_width_u", "std::vector<float>", &_allhit_width_u);
  _calo_tree->Branch("allhit_width_v", "std::vector<float>", &_allhit_width_v);
  _calo_tree->Branch("allhit_width_y", "std::vector<float>", &_allhit_width_y);
  _calo_tree->Branch("allhit_time_u", "std::vector<float>", &_allhit_time_u);
  _calo_tree->Branch("allhit_time_v", "std::vector<float>", &_allhit_time_v);
  _calo_tree->Branch("allhit_time_y", "std::vector<float>", &_allhit_time_y);
  _calo_tree->Branch("allhit_nhit_u", "std::vector<unsigned>", &_allhit_nhit_u);
  _calo_tree->Branch("allhit_nhit_v", "std::vector<unsigned>", &_allhit_nhit_v);
  _calo_tree->Branch("allhit_nhit_y", "std::vector<unsigned>", &_allhit_nhit_y);

  _calo_tree->Branch("trkhit_charge_u", "std::vector<float>", &_trkhit_charge_u);
  _calo_tree->Branch("trkhit_charge_v", "std::vector<float>", &_trkhit_charge_v);
  _calo_tree->Branch("trkhit_charge_y", "std::vector<float>", &_trkhit_charge_y);
  _calo_tree->Branch("trkhit_wire_u", "std::vector<unsigned>", &_trkhit_wire_u);
  _calo_tree->Branch("trkhit_wire_v", "std::vector<unsigned>", &_trkhit_wire_v);
  _calo_tree->Branch("trkhit_wire_y", "std::vector<unsigned>", &_trkhit_wire_y);
  _calo_tree->Branch("trkhit_channel_u", "std::vector<unsigned>", &_trkhit_channel_u);
  _calo_tree->Branch("trkhit_channel_v", "std::vector<unsigned>", &_trkhit_channel_v);
  _calo_tree->Branch("trkhit_channel_y", "std::vector<unsigned>", &_trkhit_channel_y);
  _calo_tree->Branch("trkhit_multiplicity_u", "std::vector<int>", &_trkhit_multiplicity_u); 
  _calo_tree->Branch("trkhit_multiplicity_v", "std::vector<int>", &_trkhit_multiplicity_v); 
  _calo_tree->Branch("trkhit_multiplicity_y", "std::vector<int>", &_trkhit_multiplicity_y); 
  _calo_tree->Branch("trkhit_width_u", "std::vector<float>", &_trkhit_width_u);
  _calo_tree->Branch("trkhit_width_v", "std::vector<float>", &_trkhit_width_v);
  _calo_tree->Branch("trkhit_width_y", "std::vector<float>", &_trkhit_width_y);
  _calo_tree->Branch("trkhit_time_u", "std::vector<float>", &_trkhit_time_u);
  _calo_tree->Branch("trkhit_time_v", "std::vector<float>", &_trkhit_time_v);
  _calo_tree->Branch("trkhit_time_y", "std::vector<float>", &_trkhit_time_y);
  _calo_tree->Branch("trkhit_nhit_u", "std::vector<unsigned>", &_trkhit_nhit_u);
  _calo_tree->Branch("trkhit_nhit_v", "std::vector<unsigned>", &_trkhit_nhit_v);
  _calo_tree->Branch("trkhit_nhit_y", "std::vector<unsigned>", &_trkhit_nhit_y);

  _calo_tree->Branch("calohit_charge_u", "std::vector<float>", &_calohit_charge_u);
  _calo_tree->Branch("calohit_charge_v", "std::vector<float>", &_calohit_charge_v);
  _calo_tree->Branch("calohit_charge_y", "std::vector<float>", &_calohit_charge_y);
  _calo_tree->Branch("calohit_wire_u", "std::vector<unsigned>", &_calohit_wire_u);
  _calo_tree->Branch("calohit_wire_v", "std::vector<unsigned>", &_calohit_wire_v);
  _calo_tree->Branch("calohit_wire_y", "std::vector<unsigned>", &_calohit_wire_y);
  _calo_tree->Branch("calohit_channel_u", "std::vector<unsigned>", &_calohit_channel_u);
  _calo_tree->Branch("calohit_channel_v", "std::vector<unsigned>", &_calohit_channel_v);
  _calo_tree->Branch("calohit_channel_y", "std::vector<unsigned>", &_calohit_channel_y);
  _calo_tree->Branch("calohit_multiplicity_u", "std::vector<int>", &_calohit_multiplicity_u); 
  _calo_tree->Branch("calohit_multiplicity_v", "std::vector<int>", &_calohit_multiplicity_v); 
  _calo_tree->Branch("calohit_multiplicity_y", "std::vector<int>", &_calohit_multiplicity_y); 
  _calo_tree->Branch("calohit_width_u", "std::vector<float>", &_calohit_width_u);
  _calo_tree->Branch("calohit_width_v", "std::vector<float>", &_calohit_width_v);
  _calo_tree->Branch("calohit_width_y", "std::vector<float>", &_calohit_width_y);
  _calo_tree->Branch("calohit_time_u", "std::vector<float>", &_calohit_time_u);
  _calo_tree->Branch("calohit_time_v", "std::vector<float>", &_calohit_time_v);
  _calo_tree->Branch("calohit_time_y", "std::vector<float>", &_calohit_time_y);
  _calo_tree->Branch("calohit_nhit_u", "std::vector<unsigned>", &_calohit_nhit_u);
  _calo_tree->Branch("calohit_nhit_v", "std::vector<unsigned>", &_calohit_nhit_v);
  _calo_tree->Branch("calohit_nhit_y", "std::vector<unsigned>", &_calohit_nhit_y);

  _calo_tree->Branch("sumhit_charge_u", "std::vector<float>", &_sumhit_charge_u);
  _calo_tree->Branch("sumhit_charge_v", "std::vector<float>", &_sumhit_charge_v);
  _calo_tree->Branch("sumhit_charge_y", "std::vector<float>", &_sumhit_charge_y);
  _calo_tree->Branch("sumhit_wire_u", "std::vector<unsigned>", &_sumhit_wire_u);
  _calo_tree->Branch("sumhit_wire_v", "std::vector<unsigned>", &_sumhit_wire_v);
  _calo_tree->Branch("sumhit_wire_y", "std::vector<unsigned>", &_sumhit_wire_y);
  _calo_tree->Branch("sumhit_channel_u", "std::vector<unsigned>", &_sumhit_channel_u);
  _calo_tree->Branch("sumhit_channel_v", "std::vector<unsigned>", &_sumhit_channel_v);
  _calo_tree->Branch("sumhit_channel_y", "std::vector<unsigned>", &_sumhit_channel_y);
  _calo_tree->Branch("sumhit_multiplicity_u", "std::vector<int>", &_sumhit_multiplicity_u); 
  _calo_tree->Branch("sumhit_multiplicity_v", "std::vector<int>", &_sumhit_multiplicity_v); 
  _calo_tree->Branch("sumhit_multiplicity_y", "std::vector<int>", &_sumhit_multiplicity_y); 
  _calo_tree->Branch("sumhit_width_u", "std::vector<float>", &_sumhit_width_u);
  _calo_tree->Branch("sumhit_width_v", "std::vector<float>", &_sumhit_width_v);
  _calo_tree->Branch("sumhit_width_y", "std::vector<float>", &_sumhit_width_y);
  _calo_tree->Branch("sumhit_time_u", "std::vector<float>", &_sumhit_time_u);
  _calo_tree->Branch("sumhit_time_v", "std::vector<float>", &_sumhit_time_v);
  _calo_tree->Branch("sumhit_time_y", "std::vector<float>", &_sumhit_time_y);
  _calo_tree->Branch("sumhit_nhit_u", "std::vector<unsigned>", &_sumhit_nhit_u);
  _calo_tree->Branch("sumhit_nhit_v", "std::vector<unsigned>", &_sumhit_nhit_v);
  _calo_tree->Branch("sumhit_nhit_y", "std::vector<unsigned>", &_sumhit_nhit_y);

  _calo_tree->Branch("backtracked_end_x", &_backtracked_end_x, "backtracked_end_x/F");
  _calo_tree->Branch("backtracked_end_y", &_backtracked_end_y, "backtracked_end_y/F");
  _calo_tree->Branch("backtracked_end_z", &_backtracked_end_z, "backtracked_end_z/F");

  _calo_tree->Branch("backtracked_start_x", &_backtracked_start_x, "backtracked_start_x/F");
  _calo_tree->Branch("backtracked_start_y", &_backtracked_start_y, "backtracked_start_y/F");
  _calo_tree->Branch("backtracked_start_z", &_backtracked_start_z, "backtracked_start_z/F");
  _calo_tree->Branch("backtracked_start_t", &_backtracked_start_t, "backtracked_start_t/F");
  // _calo_tree->Branch("backtracked_start_U", &_backtracked_start_U, "backtracked_start_U/F");
  // _calo_tree->Branch("backtracked_start_V", &_backtracked_start_V, "backtracked_start_V/F");
  // _calo_tree->Branch("backtracked_start_Y", &_backtracked_start_Y, "backtracked_start_Y/F");
  _calo_tree->Branch("backtracked_sce_start_x", &_backtracked_sce_start_x, "backtracked_sce_start_x/F");
  _calo_tree->Branch("backtracked_sce_start_y", &_backtracked_sce_start_y, "backtracked_sce_start_y/F");
  _calo_tree->Branch("backtracked_sce_start_z", &_backtracked_sce_start_z, "backtracked_sce_start_z/F");
  // _calo_tree->Branch("backtracked_sce_start_U", &_backtracked_sce_start_U, "backtracked_sce_start_U/F");
  // _calo_tree->Branch("backtracked_sce_start_V", &_backtracked_sce_start_V, "backtracked_sce_start_V/F");
  // _calo_tree->Branch("backtracked_sce_start_Y", &_backtracked_sce_start_Y, "backtracked_sce_start_Y/F");

  _calo_tree->Branch("backtracked_process_is_stopping", &_backtracked_process_is_stopping, "backtracked_process_is_stopping/O");
  _calo_tree->Branch("backtracked_end_in_tpc", &_backtracked_end_in_tpc, "backtracked_end_in_tpc/O");

  _calo_tree->Branch("generation", &_generation, "generation/i");
  _calo_tree->Branch("trk_daughters", &_trk_daughters, "trk_daughters/i");
  _calo_tree->Branch("daughters", &_daughters, "daughters/i");
  _calo_tree->Branch("shr_daughters", &_shr_daughters, "shr_daughters/i");
  _calo_tree->Branch("n_pfp", &_n_pfp, "n_pfp/i");

  // track information
  _calo_tree->Branch("trk_score", &_trk_score, "trk_score/F");

  _calo_tree->Branch("trk_theta", &_trk_theta, "trk_theta/F");
  _calo_tree->Branch("trk_phi", &_trk_phi, "trk_phi/F");
  _calo_tree->Branch("trk_len", &_trk_len, "trk_len/F");

  _calo_tree->Branch("trk_calo_range_u", &_trk_calo_range_u, "trk_calo_range_u/F");
  _calo_tree->Branch("trk_calo_range_v", &_trk_calo_range_v, "trk_calo_range_v/F");
  _calo_tree->Branch("trk_calo_range_y", &_trk_calo_range_y, "trk_calo_range_y/F");

  _calo_tree->Branch("trk_dir_x", &_trk_dir_x, "trk_dir_x/F");
  _calo_tree->Branch("trk_dir_y", &_trk_dir_y, "trk_dir_y/F");
  _calo_tree->Branch("trk_dir_z", &_trk_dir_z, "trk_dir_z/F");

  _calo_tree->Branch("longest", &_longest, "longest/I");

  _calo_tree->Branch("trk_start_x", &_trk_start_x, "trk_start_x/F");
  _calo_tree->Branch("trk_start_y", &_trk_start_y, "trk_start_y/F");
  _calo_tree->Branch("trk_start_z", &_trk_start_z, "trk_start_z/F");

  _calo_tree->Branch("trk_sce_start_x", &_trk_sce_start_x, "trk_sce_start_x/F");
  _calo_tree->Branch("trk_sce_start_y", &_trk_sce_start_y, "trk_sce_start_y/F");
  _calo_tree->Branch("trk_sce_start_z", &_trk_sce_start_z, "trk_sce_start_z/F");

  _calo_tree->Branch("trk_end_x", &_trk_end_x, "trk_end_x/F");
  _calo_tree->Branch("trk_end_y", &_trk_end_y, "trk_end_y/F");
  _calo_tree->Branch("trk_end_z", &_trk_end_z, "trk_end_z/F");

  _calo_tree->Branch("trk_sce_end_x", &_trk_sce_end_x, "trk_sce_end_x/F");
  _calo_tree->Branch("trk_sce_end_y", &_trk_sce_end_y, "trk_sce_end_y/F");
  _calo_tree->Branch("trk_sce_end_z", &_trk_sce_end_z, "trk_sce_end_z/F");

  _calo_tree->Branch("trk_bragg_p_u", &_trk_bragg_p_u, "trk_bragg_p_u/F");
  _calo_tree->Branch("trk_bragg_mu_u", &_trk_bragg_mu_u, "trk_bragg_mu_u/F");
  _calo_tree->Branch("trk_bragg_mip_u", &_trk_bragg_mip_u, "trk_bragg_mip_u/F");
  _calo_tree->Branch("trk_pid_chipr_u", &_trk_pid_chipr_u, "trk_pid_chipr_u/F");
  _calo_tree->Branch("trk_pid_chika_u", &_trk_pid_chika_u, "trk_pid_chika_u/F");
  _calo_tree->Branch("trk_pid_chipi_u", &_trk_pid_chipi_u, "trk_pid_chipi_u/F");
  _calo_tree->Branch("trk_pid_chimu_u", &_trk_pid_chimu_u, "trk_pid_chimu_u/F");
  _calo_tree->Branch("trk_pida_u", &_trk_pida_u, "trk_pida_u/F");

  _calo_tree->Branch("trk_bragg_p_v", &_trk_bragg_p_v, "trk_bragg_p_v/F");
  _calo_tree->Branch("trk_bragg_mu_v", &_trk_bragg_mu_v, "trk_bragg_mu_v/F");
  _calo_tree->Branch("trk_bragg_mip_v", &_trk_bragg_mip_v, "trk_bragg_mip_v/F");
  _calo_tree->Branch("trk_pid_chipr_v", &_trk_pid_chipr_v, "trk_pid_chipr_v/F");
  _calo_tree->Branch("trk_pid_chika_v", &_trk_pid_chika_v, "trk_pid_chika_v/F");
  _calo_tree->Branch("trk_pid_chipi_v", &_trk_pid_chipi_v, "trk_pid_chipi_v/F");
  _calo_tree->Branch("trk_pid_chimu_v", &_trk_pid_chimu_v, "trk_pid_chimu_v/F");
  _calo_tree->Branch("trk_pida_v", &_trk_pida_v, "trk_pida_v/F");

  _calo_tree->Branch("trk_bragg_p_y", &_trk_bragg_p_y, "trk_bragg_p_y/F");
  _calo_tree->Branch("trk_bragg_mu_y", &_trk_bragg_mu_y, "trk_bragg_mu_y/F");
  _calo_tree->Branch("trk_bragg_mip_y", &_trk_bragg_mip_y, "trk_bragg_mip_y/F");
  _calo_tree->Branch("trk_pid_chipr_y", &_trk_pid_chipr_y, "trk_pid_chipr_y/F");
  _calo_tree->Branch("trk_pid_chika_y", &_trk_pid_chika_y, "trk_pid_chika_y/F");
  _calo_tree->Branch("trk_pid_chipi_y", &_trk_pid_chipi_y, "trk_pid_chipi_y/F");
  _calo_tree->Branch("trk_pid_chimu_y", &_trk_pid_chimu_y, "trk_pid_chimu_y/F");
  _calo_tree->Branch("trk_pida_y", &_trk_pida_y, "trk_pida_y/F");

  _calo_tree->Branch("trk_bragg_p_three_planes", &_trk_bragg_p_three_planes, "trk_bragg_p_three_planes/F");

  _calo_tree->Branch("trk_mcs_muon_mom", &_trk_mcs_muon_mom, "trk_mcs_muon_mom/F");
  _calo_tree->Branch("trk_energy_proton", &_trk_energy_proton, "trk_energy_proton/F");
  _calo_tree->Branch("trk_energy_muon", &_trk_energy_muon, "trk_energy_muon/F");

  // dedx vector
  _calo_tree->Branch("dqdx_u", "std::vector<float>", &_dqdx_u);
  _calo_tree->Branch("dqdx_v", "std::vector<float>", &_dqdx_v);
  _calo_tree->Branch("dqdx_y", "std::vector<float>", &_dqdx_y);

  _calo_tree->Branch("dedx_u", "std::vector<float>", &_dedx_u);
  _calo_tree->Branch("dedx_v", "std::vector<float>", &_dedx_v);
  _calo_tree->Branch("dedx_y", "std::vector<float>", &_dedx_y);

  _calo_tree->Branch("dedx_channel_u", "std::vector<unsigned>", &_dedx_channel_u);
  _calo_tree->Branch("dedx_channel_v", "std::vector<unsigned>", &_dedx_channel_v);
  _calo_tree->Branch("dedx_channel_y", "std::vector<unsigned>", &_dedx_channel_y);

  _calo_tree->Branch("rr_u", "std::vector<float>", &_rr_u);
  _calo_tree->Branch("rr_v", "std::vector<float>", &_rr_v);
  _calo_tree->Branch("rr_y", "std::vector<float>", &_rr_y);

  _calo_tree->Branch("pitch_u", "std::vector<float>", &_pitch_u);
  _calo_tree->Branch("pitch_v", "std::vector<float>", &_pitch_v);
  _calo_tree->Branch("pitch_y", "std::vector<float>", &_pitch_y);

  _calo_tree->Branch("sx", "std::vector<float>", & _sx); 
  _calo_tree->Branch("sy", "std::vector<float>", & _sy); 
  _calo_tree->Branch("sz", "std::vector<float>", & _sz); 

  _calo_tree->Branch("w_u", "std::vector<float>", &_w_u);
  _calo_tree->Branch("w_v", "std::vector<float>", &_w_v);
  _calo_tree->Branch("w_y", "std::vector<float>", &_w_y);

  _calo_tree->Branch("t_u", "std::vector<float>", &_t_u);
  _calo_tree->Branch("t_v", "std::vector<float>", &_t_v);
  _calo_tree->Branch("t_y", "std::vector<float>", &_t_y);

  _calo_tree->Branch("x_u", "std::vector<float>", &_x_u);
  _calo_tree->Branch("x_v", "std::vector<float>", &_x_v);
  _calo_tree->Branch("x_y", "std::vector<float>", &_x_y);

  _calo_tree->Branch("y_u", "std::vector<float>", &_y_u);
  _calo_tree->Branch("y_v", "std::vector<float>", &_y_v);
  _calo_tree->Branch("y_y", "std::vector<float>", &_y_y);

  _calo_tree->Branch("z_u", "std::vector<float>", &_z_u);
  _calo_tree->Branch("z_v", "std::vector<float>", &_z_v);
  _calo_tree->Branch("z_y", "std::vector<float>", &_z_y);

  _calo_tree->Branch("dir_x_u", "std::vector<float>", &_dir_x_u);
  _calo_tree->Branch("dir_x_v", "std::vector<float>", &_dir_x_v);
  _calo_tree->Branch("dir_x_y", "std::vector<float>", &_dir_x_y);

  _calo_tree->Branch("dir_y_u", "std::vector<float>", &_dir_y_u);
  _calo_tree->Branch("dir_y_v", "std::vector<float>", &_dir_y_v);
  _calo_tree->Branch("dir_y_y", "std::vector<float>", &_dir_y_y);

  _calo_tree->Branch("dir_z_u", "std::vector<float>", &_dir_z_u);
  _calo_tree->Branch("dir_z_v", "std::vector<float>", &_dir_z_v);
  _calo_tree->Branch("dir_z_y", "std::vector<float>", &_dir_z_y);
}

void CalorimetryAnalysis::resetTTree(TTree *_tree)
{
}

std::vector<unsigned> FinddEdxHits(const anab::Calorimetry &calo, const std::vector<art::Ptr<recob::Hit>> &hits) {
  std::vector<unsigned> ret;

  // NOTE: TpIndices are __actually__ the art::Ptr keys of the Hit's
  // Don't ask me why
  for (size_t ind: calo.TpIndices()) {
    bool found = false;
    for (const art::Ptr<recob::Hit> &hit: hits) {
      if (hit.key() == ind) {
        ret.push_back(hit->Channel());
        found = true;
        break;
      }
    }
    if (!found) {
      throw cet::exception("dEdx w/out Hit!");
    }
  }
  return ret;
}

float calcPitch(const geo::GeometryCore *geo, const geo::PlaneID &plane, TVector3 location, TVector3 direction, const spacecharge::SpaceCharge *SCEService=NULL) {
  float angletovert = geo->WireAngleToVertical(geo->View(plane), plane) - 0.5*::util::pi<>();

  // apply space charge change to dir
  if (SCEService) {
    geo::Point_t location_p(location);
    geo::Vector_t offset_v = SCEService->GetPosOffsets(location_p);
    TVector3 offset = TVector3(offset_v.X(), offset_v.Y(), offset_v.Z());
    // fix the x-sign
    if (location.X() < 0.) offset.SetX(-offset.X());

    // location a ~wire's distance along the track 
    TVector3 location_dx = location + geo->WirePitch(plane) * direction;

    // Apply space charge to the offset location
    geo::Point_t location_dx_p(location_dx);
    geo::Vector_t offset_dx_v = SCEService->GetPosOffsets(location_dx_p);
    TVector3 offset_dx(offset_dx_v.X(), offset_dx_v.Y(), offset_dx_v.Z());
    // fix the x-sign
    if (location_dx.X() < 0.) offset_dx.SetX(-offset_dx.X());

    // this is the actual direction
    direction = (offset_dx + location_dx - offset - location).Unit();
  }

  // get the pitch
  float cosgamma = abs(cos(angletovert) * direction.Z() + sin(angletovert) * direction.Y());
  float pitch = geo->WirePitch(plane) / cosgamma;

  // map pitch back onto particle trajectory
  if (SCEService) {
    geo::Point_t location_p(location);
    geo::Vector_t offset_v = SCEService->GetPosOffsets(location_p);
    TVector3 offset(offset_v.X(), offset_v.Y(), offset_v.Z());
    // fix the x-sign
    if (location.X() < 0.) offset.SetX(-offset.X());

    TVector3 location_sce = location + offset;

    TVector3 location_sce_dx = location_sce + direction * pitch;

    // map this back onto the particle trajectory
    geo::Point_t location_sce_dx_p(location_sce_dx);
    geo::Vector_t back_v = SCEService->GetCalPosOffsets(location_sce_dx_p, plane.TPC);
    TVector3 back(back_v.X(), back_v.Y(), back_v.Z());
    TVector3 location_dx = location_sce_dx + back;

    pitch = (location - location_dx).Mag();
  }

  return pitch;
}


void CalorimetryAnalysis::FillCalorimetry(art::Event const &e,
            const std::vector<art::Ptr<recob::Wire>> &allWires,
            const std::vector<art::Ptr<recob::Hit>> &allHits,
            const art::Ptr<recob::PFParticle> &pfp,
            const art::Ptr<recob::Track> &trk,
            const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints,
            const std::vector<art::Ptr<anab::Calorimetry>> &calos,
            const std::vector<art::Ptr<anab::ParticleID>> &pids,
            const std::vector<art::Ptr<recob::Hit>> &areaHits,
            const std::vector<art::Ptr<recob::Hit>> &trkHits,
            const std::vector<art::Ptr<recob::Hit>> &caloHits,
            const recob::MCSFitResult &muon_mcs,
            const sbn::RangeP &muon_range,
            const sbn::RangeP &proton_range,
            const bool fData,
            const std::vector<std::pair<int, float>> &particleMatches,
            const art::Ptr<simb::MCParticle> &trueParticle,
            const std::array<std::map<unsigned, std::vector<const sim::IDE*>>, 3> particle_ide_map,
            const std::array<std::vector<const sim::IDE*>,3> &allIDEs,
            const std::vector<geo::BoxBoundedGeo> &activeVolumes)
            //const std::vector<searchingfornues::BtPart> btparts_v,
            //const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart)
{

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  const spacecharge::SpaceCharge *spacecharge = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const dprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  // TODO: variables to set
  _generation = -1;
  _shr_daughters = 0;
  _trk_daughters = 0;
  _trk_score = -1;
  _longest = -1;

  _daughters = pfp->NumDaughters();
  _pfpno = pfp->Self();
  _trkid = trk->ID();

  for (unsigned sp_i = 0; sp_i < spacepoints.size(); sp_i++) {
    _sx.push_back(spacepoints[sp_i]->XYZ()[0]);
    _sy.push_back(spacepoints[sp_i]->XYZ()[1]);
    _sz.push_back(spacepoints[sp_i]->XYZ()[2]);
  }

  for (const auto &wire: allWires) {
    unsigned plane_id = geometry->ChannelToWire(wire->Channel()).at(0).Plane;

    if (plane_id == 0) {
      _integrated_charge_u.push_back(0);
      for (auto const &range: wire->SignalROI().get_ranges()) {
        for (float val: range) _integrated_charge_u.back() += val;
      }
    }
    else if (plane_id == 1) {
      _integrated_charge_v.push_back(0);
      for (auto const &range: wire->SignalROI().get_ranges()) {
        for (float val: range) _integrated_charge_v.back() += val;
      }
    }
    else if (plane_id == 2) {
      _integrated_charge_y.push_back(0);
      for (auto const &range: wire->SignalROI().get_ranges()) {
        for (float val: range) _integrated_charge_y.back() += val;
      }
    }
  }

  FillHits(clock_data, allHits, 
    _allhit_charge_u,
    _allhit_charge_v,
    _allhit_charge_y,
    _allhit_wire_u,
    _allhit_wire_v,
    _allhit_wire_y,
    _allhit_channel_u,
    _allhit_channel_v,
    _allhit_channel_y,
    _allhit_multiplicity_u,
    _allhit_multiplicity_v,
    _allhit_multiplicity_y,
    _allhit_width_u,
    _allhit_width_v,
    _allhit_width_y,
    _allhit_time_u,
    _allhit_time_v,
    _allhit_time_y,
    _allhit_nhit_u,
    _allhit_nhit_v,
    _allhit_nhit_y
  );

  FillHits(clock_data, areaHits, 
    _areahit_charge_u,
    _areahit_charge_v,
    _areahit_charge_y,
    _areahit_wire_u,
    _areahit_wire_v,
    _areahit_wire_y,
    _areahit_channel_u,
    _areahit_channel_v,
    _areahit_channel_y,
    _areahit_multiplicity_u,
    _areahit_multiplicity_v,
    _areahit_multiplicity_y,
    _areahit_width_u,
    _areahit_width_v,
    _areahit_width_y,
    _areahit_time_u,
    _areahit_time_v,
    _areahit_time_y,
    _areahit_nhit_u,
    _areahit_nhit_v,
    _areahit_nhit_y
  );

  FillHits(clock_data, trkHits, 
    _trkhit_charge_u,
    _trkhit_charge_v,
    _trkhit_charge_y,
    _trkhit_wire_u,
    _trkhit_wire_v,
    _trkhit_wire_y,
    _trkhit_channel_u,
    _trkhit_channel_v,
    _trkhit_channel_y,
    _trkhit_multiplicity_u,
    _trkhit_multiplicity_v,
    _trkhit_multiplicity_y,
    _trkhit_width_u,
    _trkhit_width_v,
    _trkhit_width_y,
    _trkhit_time_u,
    _trkhit_time_v,
    _trkhit_time_y,
    _trkhit_nhit_u,
    _trkhit_nhit_v,
    _trkhit_nhit_y
  );
    
  FillHits(clock_data, caloHits, 
    _calohit_charge_u,
    _calohit_charge_v,
    _calohit_charge_y,
    _calohit_wire_u,
    _calohit_wire_v,
    _calohit_wire_y,
    _calohit_channel_u,
    _calohit_channel_v,
    _calohit_channel_y,
    _calohit_multiplicity_u,
    _calohit_multiplicity_v,
    _calohit_multiplicity_y,
    _calohit_width_u,
    _calohit_width_v,
    _calohit_width_y,
    _calohit_time_u,
    _calohit_time_v,
    _calohit_time_y,
    _calohit_nhit_u,
    _calohit_nhit_v,
    _calohit_nhit_y
  );

  FillHits(clock_data, caloHits, 
    _sumhit_charge_u,
    _sumhit_charge_v,
    _sumhit_charge_y,
    _sumhit_wire_u,
    _sumhit_wire_v,
    _sumhit_wire_y,
    _sumhit_channel_u,
    _sumhit_channel_v,
    _sumhit_channel_y,
    _sumhit_multiplicity_u,
    _sumhit_multiplicity_v,
    _sumhit_multiplicity_y,
    _sumhit_width_u,
    _sumhit_width_v,
    _sumhit_width_y,
    _sumhit_time_u,
    _sumhit_time_v,
    _sumhit_time_y,
    _sumhit_nhit_u,
    _sumhit_nhit_v,
    _sumhit_nhit_y,
    false);

  // store Backtracking
  if (!fData)
  {
    _backtracked_e = trueParticle->E();
    _backtracked_pdg = trueParticle->PdgCode();
    float trueTrackE = 0.;
    float matchE = -1;
    for (auto pair: particleMatches) {
      trueTrackE += pair.second;
      if (pair.first == trueParticle->TrackId()) matchE = pair.second;
    }
    assert(matchE > 0.);
    _backtracked_purity = matchE / trueTrackE;

    float trueParticleE = 0.;
    _backtracked_q_u = 0.;
    _backtracked_q_v = 0.;
    _backtracked_q_y = 0.;
    _backtracked_qtotal_u = 0.;
    _backtracked_qtotal_v = 0.;
    _backtracked_qtotal_y = 0.;
    _backtracked_deposited_e_u = 0.;
    _backtracked_deposited_e_v = 0.;
    _backtracked_deposited_e_y = 0.;

    for (auto const &ide_pair: particle_ide_map[0]) {
      _backtracked_c_u.push_back(ide_pair.first);
      _backtracked_pc_q_u.push_back(0);
      _backtracked_pc_e_u.push_back(0);
      for (const sim::IDE *ide: ide_pair.second) {
        trueParticleE += ide->energy;
        _backtracked_q_u += ide->numElectrons;
        _backtracked_deposited_e_u += ide->energy;
        _backtracked_pc_q_u.back() += ide->numElectrons;
        _backtracked_pc_e_u.back() += ide->energy;
      }
    }

    for (auto const &ide_pair: particle_ide_map[1]) {
      _backtracked_c_v.push_back(ide_pair.first);
      _backtracked_pc_q_v.push_back(0);
      _backtracked_pc_e_v.push_back(0);
      for (const sim::IDE *ide: ide_pair.second) {
        trueParticleE += ide->energy;
        _backtracked_q_v += ide->numElectrons;
        _backtracked_deposited_e_v += ide->energy;
        _backtracked_pc_q_v.back() += ide->numElectrons;
        _backtracked_pc_e_v.back() += ide->energy;
      }
    }

    for (auto const &ide_pair: particle_ide_map[2]) {
      _backtracked_c_y.push_back(ide_pair.first);
      _backtracked_pc_q_y.push_back(0);
      _backtracked_pc_e_y.push_back(0);
      for (const sim::IDE *ide: ide_pair.second) {
        trueParticleE += ide->energy;
        _backtracked_q_y += ide->numElectrons;
        _backtracked_deposited_e_y += ide->energy;
        _backtracked_pc_q_y.back() += ide->numElectrons;
        _backtracked_pc_e_y.back() += ide->energy;
      }
    }

    for (auto ide: allIDEs[0]) {
      _backtracked_qtotal_u += ide->numElectrons;
    }
    for (auto ide: allIDEs[1]) {
      _backtracked_qtotal_v += ide->numElectrons;
    }
    for (auto ide: allIDEs[2]) {
      _backtracked_qtotal_y += ide->numElectrons;
    }

    _backtracked_completeness = matchE / trueParticleE;

    // TODO: set overlay purity
    
   _backtracked_theta = trueParticle->Momentum().Theta();
   _backtracked_phi = trueParticle->Momentum().Phi();

    _backtracked_px = trueParticle->Momentum().Px();
    _backtracked_py = trueParticle->Momentum().Py();
    _backtracked_pz = trueParticle->Momentum().Pz();

    _backtracked_end_x = trueParticle->EndPosition().X();
    _backtracked_end_y = trueParticle->EndPosition().Y();
    _backtracked_end_z = trueParticle->EndPosition().Z();

    _backtracked_start_x = trueParticle->Position().X();
    _backtracked_start_y = trueParticle->Position().Y();
    _backtracked_start_z = trueParticle->Position().Z();
    _backtracked_start_t = trueParticle->Position().T();

    _backtracked_start_U = XYZtoPlanecoordinate(_backtracked_start_x, _backtracked_start_y, _backtracked_start_z, 0);
    _backtracked_start_V = XYZtoPlanecoordinate(_backtracked_start_x, _backtracked_start_y, _backtracked_start_z, 1);
    _backtracked_start_Y = XYZtoPlanecoordinate(_backtracked_start_x, _backtracked_start_y, _backtracked_start_z, 2);

    // save the trajectory
    for (unsigned i_traj = 0; i_traj < trueParticle->NumberTrajectoryPoints(); i_traj++) {
      _backtracked_x.push_back(trueParticle->Vx(i_traj));
      _backtracked_y.push_back(trueParticle->Vy(i_traj));
      _backtracked_z.push_back(trueParticle->Vz(i_traj));

     _backtracked_dirx.push_back(trueParticle->Px(i_traj) / trueParticle->P(i_traj));
     _backtracked_diry.push_back(trueParticle->Py(i_traj) / trueParticle->P(i_traj));
     _backtracked_dirz.push_back(trueParticle->Pz(i_traj) / trueParticle->P(i_traj));
    }

    // also map between wires and the trajectory
    //
    // First setup all the traj containers
    _backtracked_nloc_c_u.insert(_backtracked_nloc_c_u.end(), _backtracked_c_u.size(), 0);
    _backtracked_locx_c_u.insert(_backtracked_locx_c_u.end(), _backtracked_c_u.size(), 0.);
    _backtracked_locy_c_u.insert(_backtracked_locy_c_u.end(), _backtracked_c_u.size(), 0.);
    _backtracked_locz_c_u.insert(_backtracked_locz_c_u.end(), _backtracked_c_u.size(), 0.);
    _backtracked_dirx_c_u.insert(_backtracked_dirx_c_u.end(), _backtracked_c_u.size(), 0.);
    _backtracked_diry_c_u.insert(_backtracked_diry_c_u.end(), _backtracked_c_u.size(), 0.);
    _backtracked_dirz_c_u.insert(_backtracked_dirz_c_u.end(), _backtracked_c_u.size(), 0.);
    _backtracked_de_c_u.insert(_backtracked_de_c_u.end(),   _backtracked_c_u.size(), 0.);
    _backtracked_pitch_c_u.insert(_backtracked_pitch_c_u.end(), _backtracked_c_u.size(), 0.);
    _backtracked_scepitch_c_u.insert(_backtracked_scepitch_c_u.end(), _backtracked_c_u.size(), 0.);

    _backtracked_nloc_c_v.insert(_backtracked_nloc_c_v.end(), _backtracked_c_v.size(), 0);
    _backtracked_locx_c_v.insert(_backtracked_locx_c_v.end(), _backtracked_c_v.size(), 0.);
    _backtracked_locy_c_v.insert(_backtracked_locy_c_v.end(), _backtracked_c_v.size(), 0.);
    _backtracked_locz_c_v.insert(_backtracked_locz_c_v.end(), _backtracked_c_v.size(), 0.);
    _backtracked_dirx_c_v.insert(_backtracked_dirx_c_v.end(), _backtracked_c_v.size(), 0.);
    _backtracked_diry_c_v.insert(_backtracked_diry_c_v.end(), _backtracked_c_v.size(), 0.);
    _backtracked_dirz_c_v.insert(_backtracked_dirz_c_v.end(), _backtracked_c_v.size(), 0.);
    _backtracked_de_c_v.insert(_backtracked_de_c_v.end(),   _backtracked_c_v.size(), 0.);
    _backtracked_pitch_c_v.insert(_backtracked_pitch_c_v.end(), _backtracked_c_v.size(), 0.);
    _backtracked_scepitch_c_v.insert(_backtracked_scepitch_c_v.end(), _backtracked_c_v.size(), 0.);

    _backtracked_nloc_c_y.insert(_backtracked_nloc_c_y.end(), _backtracked_c_y.size(), 0);
    _backtracked_locx_c_y.insert(_backtracked_locx_c_y.end(), _backtracked_c_y.size(), 0.);
    _backtracked_locy_c_y.insert(_backtracked_locy_c_y.end(), _backtracked_c_y.size(), 0.);
    _backtracked_locz_c_y.insert(_backtracked_locz_c_y.end(), _backtracked_c_y.size(), 0.);
    _backtracked_dirx_c_y.insert(_backtracked_dirx_c_y.end(), _backtracked_c_y.size(), 0.);
    _backtracked_diry_c_y.insert(_backtracked_diry_c_y.end(), _backtracked_c_y.size(), 0.);
    _backtracked_dirz_c_y.insert(_backtracked_dirz_c_y.end(), _backtracked_c_y.size(), 0.);
    _backtracked_de_c_y.insert(_backtracked_de_c_y.end(),   _backtracked_c_y.size(), 0.);
    _backtracked_pitch_c_y.insert(_backtracked_pitch_c_y.end(), _backtracked_c_y.size(), 0.);
    _backtracked_scepitch_c_y.insert(_backtracked_scepitch_c_y.end(), _backtracked_c_y.size(), 0.);

    for (geo::TPCID tpc: geometry->IterateTPCIDs()) {
      for (geo::PlaneID planeID: geometry->IteratePlaneIDs(tpc)) {
        auto const &plane = planeID.Plane;
        if (plane > 2) continue; // protect against bad plane ID

        int lastWire = -1000000;
        int last_traj = 0;
        for (unsigned i_traj = 0; i_traj < trueParticle->NumberTrajectoryPoints()-1; i_traj++) {
          TVector3 thispoint = trueParticle->Position(i_traj).Vect();
          TVector3 nextpoint = trueParticle->Position(i_traj+1).Vect();

          geo::Point_t thispoint_p(thispoint.X(), thispoint.Y(), thispoint.Z());
          geo::Point_t nextpoint_p(nextpoint.X(), nextpoint.Y(), nextpoint.Z());

          double thiswirecoord = geometry->WireCoordinate(thispoint_p, planeID);
          double nextwirecoord = geometry->WireCoordinate(nextpoint_p, planeID);
          int wireStart = std::nearbyint((thiswirecoord >= nextwirecoord) ? std::floor(thiswirecoord) : std::ceil(thiswirecoord));
          int wireEnd   = std::nearbyint((thiswirecoord >= nextwirecoord) ? std::floor(nextwirecoord) : std::ceil(nextwirecoord));
          // if we're not crossing a wire continue
          if (wireStart == wireEnd) continue;

          // check the validity of the range of wires
          if (!(wireStart >= 0 && wireStart < (int)geometry->Plane(planeID).Nwires())) continue;
          if (!(wireEnd >= 0 && wireEnd <= (int)geometry->Plane(planeID).Nwires())) continue;

          // if this wire is the same as the last skip it
          if (wireStart == lastWire) {
            if (wireStart < wireEnd) wireStart ++;
            else wireStart --;
          }
          // again continue if this wouldn't cross a new wire
          if (wireStart == wireEnd) continue;

          // load the trajectory information at this point
          double locX = trueParticle->Position(i_traj).X();
          double locY = trueParticle->Position(i_traj).Y();
          double locZ = trueParticle->Position(i_traj).Z();
          double dirX = trueParticle->Momentum(i_traj).Vect().Unit().X();
          double dirY = trueParticle->Momentum(i_traj).Vect().Unit().Y();
          double dirZ = trueParticle->Momentum(i_traj).Vect().Unit().Z();
          double deltaE = (trueParticle->Momentum(last_traj).E() - trueParticle->Momentum(i_traj+1).E()) / abs(wireEnd - wireStart); // average dE over wires

          double pitch = calcPitch(geometry, planeID, trueParticle->Position(i_traj).Vect(), trueParticle->Momentum(i_traj).Vect().Unit());
          double pitchSCE = calcPitch(geometry, planeID, trueParticle->Position(i_traj).Vect(), trueParticle->Momentum(i_traj).Vect().Unit(), spacecharge);

          // save
          int incl = wireStart < wireEnd ? 1: -1;
          for (int wire = wireStart; wire != wireEnd; wire += incl) {
            lastWire = wire;
            unsigned channel = geometry->PlaneWireToChannel(planeID.Plane, wire, planeID.TPC, planeID.Cryostat);

            if (plane == 0) {
              int index = -1;
              for (unsigned i_test = 0; i_test < _backtracked_c_u.size(); i_test++) {
                if (_backtracked_c_u[i_test] == channel) {
                  index = i_test;
                  break;
                }
              }
              if (index < 0) continue;

              _backtracked_nloc_c_u[index] ++;
              _backtracked_locx_c_u[index] = locX;
              _backtracked_locy_c_u[index] = locY;
              _backtracked_locz_c_u[index] = locZ;
              _backtracked_dirx_c_u[index] = dirX;
              _backtracked_diry_c_u[index] = dirY;
              _backtracked_dirz_c_u[index] = dirZ;
              _backtracked_de_c_u[index] = deltaE;
              _backtracked_pitch_c_u[index] = pitch;
              _backtracked_scepitch_c_u[index] = pitchSCE;
            }
            else if (plane == 1) {
              int index = -1;
              for (unsigned i_test = 0; i_test < _backtracked_c_v.size(); i_test++) {
                if (_backtracked_c_v[i_test] == channel) {
                  index = i_test;
                  break;
                }
              }
              if (index < 0) continue;

              _backtracked_nloc_c_v[index] ++;
              _backtracked_locx_c_v[index] = locX;
              _backtracked_locy_c_v[index] = locY;
              _backtracked_locz_c_v[index] = locZ;
              _backtracked_dirx_c_v[index] = dirX;
              _backtracked_diry_c_v[index] = dirY;
              _backtracked_dirz_c_v[index] = dirZ;
              _backtracked_de_c_v[index] = deltaE;
              _backtracked_pitch_c_v[index] = pitch;
              _backtracked_scepitch_c_v[index] = pitchSCE;
            }
            else if (plane == 2) {
              int index = -1;
              for (unsigned i_test = 0; i_test < _backtracked_c_y.size(); i_test++) {
                if (_backtracked_c_y[i_test] == channel) {
                  index = i_test;
                  break;
                }
              }
              if (index < 0) continue;

              _backtracked_nloc_c_y[index] ++;
              _backtracked_locx_c_y[index] = locX;
              _backtracked_locy_c_y[index] = locY;
              _backtracked_locz_c_y[index] = locZ;
              _backtracked_dirx_c_y[index] = dirX;
              _backtracked_diry_c_y[index] = dirY;
              _backtracked_dirz_c_y[index] = dirZ;
              _backtracked_de_c_y[index] = deltaE;
              _backtracked_pitch_c_y[index] = pitch;
              _backtracked_scepitch_c_y[index] = pitchSCE;
            }

            
          }
          last_traj = i_traj;
        }
      }
    }

    // TODO: implement spacecharge
    float reco_st[3] = {_backtracked_start_x, _backtracked_start_y, _backtracked_start_z};
    // TODO: do something special for electrons and photons???
    // if (mcp.pdg == 11 || mcp.pdg == 22)
    // {
    //   reco_st[0] += searchingfornues::x_offset(mcp.start_t);
    // }
    True2RecoMappingXYZ(_backtracked_start_t, _backtracked_start_x, _backtracked_start_y, _backtracked_start_z, reco_st);
    _backtracked_sce_start_x = reco_st[0];
    _backtracked_sce_start_y = reco_st[1];
    _backtracked_sce_start_z = reco_st[2];
    
    _backtracked_sce_start_U = XYZtoPlanecoordinate(reco_st[0], reco_st[1], reco_st[2], 0);
    _backtracked_sce_start_V = XYZtoPlanecoordinate(reco_st[0], reco_st[1], reco_st[2], 1);
    _backtracked_sce_start_Y = XYZtoPlanecoordinate(reco_st[0], reco_st[1], reco_st[2], 2);
    _backtracked_process_is_stopping = (trueParticle->EndProcess() == "Decay" ||
                                        trueParticle->EndProcess() == "CoupledTransportation" ||
                                        trueParticle->EndProcess() == "FastScintillation" ||
                                        trueParticle->EndProcess() == "muMinusCaptureAtRest" ||
                                        trueParticle->EndProcess() == "LArVoxelReadoutScoringProcess");

    std::cout << "PID: " << _backtracked_pdg << " Process: " << trueParticle->EndProcess() << " is stopping: " << _backtracked_process_is_stopping << std::endl;

    std::vector<geoalgo::AABox> aa_volumes;
    _backtracked_end_in_tpc = false;
    for (const geo::BoxBoundedGeo &AV: activeVolumes) {
      if (AV.ContainsPosition(trueParticle->Position().Vect()) && AV.ContainsPosition(trueParticle->EndPosition().Vect())) {
        _backtracked_end_in_tpc = true;
        aa_volumes.emplace_back(AV.MinX(), AV.MinY(), AV.MinZ(), AV.MaxX(), AV.MaxY(), AV.MaxZ());
        break;
      }
    }

    _backtracked_length_start_to_end = (trueParticle->Position().Vect() - trueParticle->EndPosition().Vect()).Mag();
    _backtracked_length = 0.;
    for (unsigned i = 1; i < trueParticle->NumberTrajectoryPoints(); i++) {
      _backtracked_length += ContainedLength(trueParticle->Position(i-1).Vect(), trueParticle->Position(i).Vect(), aa_volumes);
    }
  }

  // Kinetic energy using tabulated stopping power (GeV)
  float mcs_momentum_muon = muon_range.range_p;
  float momentum_proton = proton_range.range_p;
  float energy_proton = std::sqrt(std::pow(momentum_proton, 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
  float energy_muon = std::sqrt(std::pow(mcs_momentum_muon, 2) + std::pow(muon->Mass(), 2)) - muon->Mass();

  _trk_mcs_muon_mom = mcs_momentum_muon;
  _trk_energy_proton = energy_proton;
  _trk_energy_muon = energy_muon;

  _trk_theta = trk->Theta();
  _trk_phi = trk->Phi();
  // TODO: SCE
  _trk_len = trk->Length(); 

  _trk_dir_x = trk->StartDirection().X();
  _trk_dir_y = trk->StartDirection().Y();
  _trk_dir_z = trk->StartDirection().Z();

  _trk_start_x = trk->Start().X();
  _trk_start_y = trk->Start().Y();
  _trk_start_z = trk->Start().Z();

  _trk_end_x = trk->End().X();
  _trk_end_y = trk->End().Y();
  _trk_end_z = trk->End().Z();

  // TODO: SCE
  float _trk_start_sce[3];
  (void) _trk_start_sce;
  _trk_sce_start_x = _trk_start_x;
  _trk_sce_start_y = _trk_start_y;
  _trk_sce_start_z = _trk_start_z;

  float _trk_end_sce[3];
  (void) _trk_end_sce;
  _trk_sce_end_x = _trk_end_x;
  _trk_sce_end_y = _trk_end_y;
  _trk_sce_end_z = _trk_end_z;

  // fill other particle ID
  for (const art::Ptr<anab::ParticleID> pid: pids) {
    auto const &plane = pid->PlaneID().Plane;
    if (plane > 2) continue;

    if (plane == 0) {
       // TODO: bragg variable??
      _trk_bragg_p_u = -1;
      _trk_bragg_mu_u = -1;
      _trk_bragg_mip_u = -1;

      _trk_pid_chipr_u = pid->Chi2Proton();
      _trk_pid_chika_u = pid->Chi2Kaon();
      _trk_pid_chipi_u = pid->Chi2Pion();
      _trk_pid_chimu_u = pid->Chi2Muon();
      _trk_pida_u = pid->PIDA();
    }
    else if (plane == 1) {
       // TODO: bragg variable??
      _trk_bragg_p_v = -1;
      _trk_bragg_mu_v = -1;
      _trk_bragg_mip_v = -1;

      _trk_pid_chipr_v = pid->Chi2Proton();
      _trk_pid_chika_v = pid->Chi2Kaon();
      _trk_pid_chipi_v = pid->Chi2Pion();
      _trk_pid_chimu_v = pid->Chi2Muon();
      _trk_pida_v = pid->PIDA();
    }
    else if (plane == 2) {
       // TODO: bragg variable??
      _trk_bragg_p_y = -1;
      _trk_bragg_mu_y = -1;
      _trk_bragg_mip_y = -1;

      _trk_pid_chipr_y = pid->Chi2Proton();
      _trk_pid_chika_y = pid->Chi2Kaon();
      _trk_pid_chipi_y = pid->Chi2Pion();
      _trk_pid_chimu_y = pid->Chi2Muon();
    }
  }

  // fill Calorimetry
  for (const art::Ptr<anab::Calorimetry> calo : calos)
  {
    // TODO: will this work for ICARUS?
    auto const& plane = calo->PlaneID().Plane;
    if (plane>2)
    {
      continue;
    }
    auto const& xyz_v = calo->XYZ();

    if (plane == 0)
    {
      _trk_calo_range_u = calo->Range();
      _dqdx_u = calo->dQdx();
      _rr_u = calo->ResidualRange();
      _pitch_u = calo->TrkPitchVec();
      for (auto xyz : xyz_v)
      {
        _w_u.push_back(geometry->WireCoordinate(xyz, calo->PlaneID()));
        _t_u.push_back(dprop.ConvertXToTicks(xyz.X(), calo->PlaneID()));
        _x_u.push_back(xyz.X());
        _y_u.push_back(xyz.Y());
        _z_u.push_back(xyz.Z());

        float _dir_u[3];
        TrkDirectionAtXYZ(*trk, xyz.X(), xyz.Y(), xyz.Z(), _dir_u);
        _dir_x_u.push_back(_dir_u[0]);
        _dir_y_u.push_back(_dir_u[1]);
        _dir_z_u.push_back(_dir_u[2]);
      }
      // TODO: ShrFit???
      _dedx_u = calo->dEdx();
      _dedx_channel_u = FinddEdxHits(*calo, trkHits);

    }
    else if (plane == 1)
    {
      _trk_calo_range_v = calo->Range();

      _dqdx_v = calo->dQdx();
      _rr_v = calo->ResidualRange();
      _pitch_v = calo->TrkPitchVec();
      for (auto xyz : xyz_v)
      {
        _w_v.push_back(geometry->WireCoordinate(xyz, calo->PlaneID()));
        _t_v.push_back(dprop.ConvertXToTicks(xyz.X(), calo->PlaneID()));
        _x_v.push_back(xyz.X());
        _y_v.push_back(xyz.Y());
        _z_v.push_back(xyz.Z());

        float _dir_v[3];
        TrkDirectionAtXYZ(*trk, xyz.X(), xyz.Y(), xyz.Z(), _dir_v);
        _dir_x_v.push_back(_dir_v[0]);
        _dir_y_v.push_back(_dir_v[1]);
        _dir_z_v.push_back(_dir_v[2]);
      }
      // TODO: ShrFit???
      _dedx_v = calo->dEdx();
      _dedx_channel_v = FinddEdxHits(*calo, trkHits);
    }
    else if (plane == 2) //collection
    {
      _trk_calo_range_y = calo->Range();

      _dqdx_y = calo->dQdx();
      _rr_y = calo->ResidualRange();
      _pitch_y = calo->TrkPitchVec();
      for (auto xyz : xyz_v)
      {
        _w_y.push_back(geometry->WireCoordinate(xyz, calo->PlaneID()));
        _t_y.push_back(dprop.ConvertXToTicks(xyz.X(), calo->PlaneID()));
        _x_y.push_back(xyz.X());
        _y_y.push_back(xyz.Y());
        _z_y.push_back(xyz.Z());

        float _dir_y[3];
        TrkDirectionAtXYZ(*trk, xyz.X(), xyz.Y(), xyz.Z(), _dir_y);
        _dir_x_y.push_back(_dir_y[0]);
        _dir_y_y.push_back(_dir_y[1]);
        _dir_z_y.push_back(_dir_y[2]);
      }
      // TODO: ShrFit???
      _dedx_y = calo->dEdx();
      _dedx_channel_y = FinddEdxHits(*calo, trkHits);
    }
  }

  _calo_tree->Fill();
}

// Define helper functions
  void TrkDirectionAtXYZ(const recob::Track trk, const double x, const double y, const double z, float out[3])
  {
    float min_dist = 1000;
    size_t i_min = -1;
    for(size_t i=0; i < trk.NumberTrajectoryPoints(); i++)
    {
      if (trk.HasValidPoint(i))
      { // check this point is valid
        auto point_i = trk.LocationAtPoint(i);
        float distance = distance3d((double)point_i.X(), (double)point_i.Y(), (double)point_i.Z(),
                  x, y, z);
        if (distance < min_dist)
        {
          min_dist = distance;
          i_min = i;
        }
      }// if point is valid
    }// for all track points

    if (i_min == (size_t)-1) return;

    auto direction = trk.DirectionAtPoint(i_min);
    out[0] = (float)direction.X();
    out[1] = (float)direction.Y();
    out[2] = (float)direction.Z();

    float norm;
    norm = out[0]*out[0] + out[1]*out[1] + out[2]*out[2];
    if (fabs(norm -1) > 0.001)
    {
      std::cout << "i_min = " << i_min << std::endl;
      std::cout << "minimum distance = " << min_dist << std::endl;
      std::cout << "out[0], out[1], out[2] = " << out[0] << " , " << out[1] << " , " << out[2] << std::endl;
      std::cout << "norm = " << norm << std::endl;
    }
  }

void FillHits(const detinfo::DetectorClocksData &dclock,
              const std::vector<art::Ptr<recob::Hit>> &hits,
              std::vector<float> &charge_u,
              std::vector<float> &charge_v,
              std::vector<float> &charge_y,
              std::vector<unsigned> &wire_u,
              std::vector<unsigned> &wire_v,
              std::vector<unsigned> &wire_y,
              std::vector<unsigned> &channel_u,
              std::vector<unsigned> &channel_v,
              std::vector<unsigned> &channel_y,
              std::vector<int> &multiplicity_u,
              std::vector<int> &multiplicity_v,
              std::vector<int> &multiplicity_y,
              std::vector<float> &width_u,
              std::vector<float> &width_v,
              std::vector<float> &width_y,
              std::vector<float> &time_u,
              std::vector<float> &time_v,
              std::vector<float> &time_y,
              std::vector<unsigned> &nhit_u,
              std::vector<unsigned> &nhit_v,
              std::vector<unsigned> &nhit_y,
              bool use_integral) {

  struct HitInfo {
    float charge;
    int multiplicity;
    float width;
    int start;
    int end;
    float time;
    int nhit;
    int wire;
    HitInfo(): charge(0), multiplicity(0), width(0), start(0), end(0), time(0.), nhit(0) {}
  };

  std::map<unsigned, HitInfo> hit_map_u;
  std::map<unsigned, HitInfo> hit_map_v;
  std::map<unsigned, HitInfo> hit_map_y;

  for (const art::Ptr<recob::Hit> hit: hits) {
    if (hit->WireID().Plane == 0) {
      // if using summed ADC method to count charge, don't double-count hits in the same snippet
      if (!use_integral && hit_map_u[hit->Channel()].start == hit->StartTick() && hit_map_u[hit->Channel()].end == hit->EndTick()) {
        continue;
      }
      hit_map_u[hit->Channel()].charge += (use_integral) ? hit->Integral() : hit->SummedADC();
      hit_map_u[hit->Channel()].multiplicity = std::max(hit_map_u[hit->Channel()].multiplicity, (int) hit->Multiplicity());
      hit_map_u[hit->Channel()].start = hit_map_u[hit->Channel()].nhit ? std::min(hit_map_u[hit->Channel()].start, hit->StartTick()) : hit->StartTick();
      hit_map_u[hit->Channel()].end = hit_map_u[hit->Channel()].nhit ? std::max(hit_map_u[hit->Channel()].end, hit->EndTick()) : hit->EndTick();
      hit_map_u[hit->Channel()].width = (hit_map_u[hit->Channel()].end - hit_map_u[hit->Channel()].start);
      hit_map_u[hit->Channel()].time = (hit->PeakTime() + hit_map_u[hit->Channel()].time * hit_map_u[hit->Channel()].nhit) / (hit_map_u[hit->Channel()].nhit + 1);
      hit_map_u[hit->Channel()].nhit ++;
      hit_map_u[hit->Channel()].wire = hit->WireID().Wire;
    }
    else if (hit->WireID().Plane == 1) {
      // if using summed ADC method to count charge, don't double-count hits in the same snippet
      if (!use_integral && hit_map_v[hit->Channel()].start == hit->StartTick() && hit_map_v[hit->Channel()].end == hit->EndTick()) {
        continue;
      }
      hit_map_v[hit->Channel()].charge += (use_integral) ? hit->Integral() : hit->SummedADC();
      hit_map_v[hit->Channel()].multiplicity = std::max(hit_map_v[hit->Channel()].multiplicity, (int) hit->Multiplicity());
      hit_map_v[hit->Channel()].start = hit_map_v[hit->Channel()].nhit ? std::min(hit_map_v[hit->Channel()].start, hit->StartTick()) : hit->StartTick();
      hit_map_v[hit->Channel()].end = hit_map_v[hit->Channel()].nhit ? std::max(hit_map_v[hit->Channel()].end, hit->EndTick()) : hit->EndTick();
      hit_map_v[hit->Channel()].width = (hit_map_v[hit->Channel()].end - hit_map_v[hit->Channel()].start);
      hit_map_v[hit->Channel()].time = (hit->PeakTime() + hit_map_v[hit->Channel()].time * hit_map_v[hit->Channel()].nhit) / (hit_map_v[hit->Channel()].nhit + 1);
      hit_map_v[hit->Channel()].nhit ++;
      hit_map_v[hit->Channel()].wire = hit->WireID().Wire;
    }
    else {
      // if using summed ADC method to count charge, don't double-count hits in the same snippet
      if (!use_integral && hit_map_y[hit->Channel()].start == hit->StartTick() && hit_map_y[hit->Channel()].end == hit->EndTick()) {
        continue;
      }
      hit_map_y[hit->Channel()].charge += (use_integral) ? hit->Integral() : hit->SummedADC();
      hit_map_y[hit->Channel()].multiplicity = std::max(hit_map_y[hit->Channel()].multiplicity, (int) hit->Multiplicity());
      hit_map_y[hit->Channel()].start = hit_map_y[hit->Channel()].nhit ? std::min(hit_map_y[hit->Channel()].start, hit->StartTick()) : hit->StartTick();
      hit_map_y[hit->Channel()].end = hit_map_y[hit->Channel()].nhit ? std::max(hit_map_y[hit->Channel()].end, hit->EndTick()) : hit->EndTick();
      hit_map_y[hit->Channel()].width = (hit_map_y[hit->Channel()].end - hit_map_y[hit->Channel()].start);
      hit_map_y[hit->Channel()].time = (hit->PeakTime() + hit_map_y[hit->Channel()].time * hit_map_y[hit->Channel()].nhit) / (hit_map_y[hit->Channel()].nhit + 1);
      hit_map_y[hit->Channel()].nhit ++;
      hit_map_y[hit->Channel()].wire = hit->WireID().Wire;
    }
  }

  for (auto const &pair: hit_map_u) {
    channel_u.push_back(pair.first);
    charge_u.push_back(pair.second.charge);
    multiplicity_u.push_back(pair.second.multiplicity);
    width_u.push_back(pair.second.width);
    time_u.push_back(pair.second.time);
    nhit_u.push_back(pair.second.nhit);
    wire_u.push_back(pair.second.wire);
  }
  for (auto const &pair: hit_map_v) {
    channel_v.push_back(pair.first);
    charge_v.push_back(pair.second.charge);
    multiplicity_v.push_back(pair.second.multiplicity);
    width_v.push_back(pair.second.width);
    time_v.push_back(pair.second.time);
    nhit_v.push_back(pair.second.nhit);
    wire_v.push_back(pair.second.wire);
  }
  for (auto const &pair: hit_map_y) {
    channel_y.push_back(pair.first);
    charge_y.push_back(pair.second.charge);
    multiplicity_y.push_back(pair.second.multiplicity);
    width_y.push_back(pair.second.width);
    time_y.push_back(pair.second.time);
    nhit_y.push_back(pair.second.nhit);
    wire_y.push_back(pair.second.wire);
  }

}

// TODO: include space charge
void True2RecoMappingXYZ(float t,float x, float y, float z, float out[3]) {
  (void) t;
  out[0] = x;
  out[1] = y;
  out[2] = z;
}

  float XYZtoPlanecoordinate(const float x, const float y, const float z, const int plane)
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    geo::Point_t p(x,y,z);
    geo::TPCID tpc = geom->FindTPCAtPosition(p);
    if (!tpc) return -1e6;
    geo::PlaneID planeID(tpc, plane);
    double _wire2cm = geom->WirePitch(planeID);

    return geom->WireCoordinate(p, planeID) * _wire2cm;
  }

  float distance3d(const float& x1, const float& y1, const float& z1,
                  const float& x2, const float& y2, const float& z2)
  {
    return sqrt((x1-x2)*(x1-x2) +
                (y1-y2)*(y1-y2) +
                (z1-z2)*(z1-z2));

  }

float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geoalgo::AABox> &boxes) {
  static const geoalgo::GeoAlgo algo;

  // if points are the same, return 0
  if ((v0 - v1).Mag() < 1e-6) return 0;

  // construct individual points
  geoalgo::Point_t p0(v0);
  geoalgo::Point_t p1(v1);

  // construct line segment
  geoalgo::LineSegment line(p0, p1);

  double length = 0;

  // total contained length is sum of lengths in all boxes
  // assuming they are non-overlapping
  for (auto const &box: boxes) {
    int n_contained = box.Contain(p0) + box.Contain(p1);
    // both points contained -- length is total length (also can break out of loop)
    if (n_contained == 2) {
      length = (v1 - v0).Mag();
      break;
    }
    // one contained -- have to find intersection point (which must exist)
    if (n_contained == 1) {
      auto intersections = algo.Intersection(line, box);
      // Because of floating point errors, it can sometimes happen
      // that there is 1 contained point but no "Intersections"
      // if one of the points is right on the edge
      if (intersections.size() == 0) {
        // determine which point is on the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        assert(p0_edge || p1_edge);
        // contained one is on edge -- can treat both as not contained
        //
        // In this case, no length
        if ((p0_edge && box.Contain(p0)) || (box.Contain(p1) && p1_edge))
          continue;
        // un-contaned one is on edge -- treat both as contained
        else if ((p0_edge && box.Contain(p1)) || (box.Contain(p0) && p1_edge)) {
	  length = (v1 - v0).Mag();
	  break;
        }
        else {
          assert(false); // bad
        }
      }
      // floating point errors can also falsely cause 2 intersection points
      //
      // in this case, one of the intersections must be very close to the 
      // "contained" point, so the total contained length will be about
      // the same as the distance between the two intersection points
      else if (intersections.size() == 2) {
        length += (intersections.at(0).ToTLorentzVector().Vect() - intersections.at(1).ToTLorentzVector().Vect()).Mag();
        continue;
      }
      // "Correct"/ideal case -- 1 intersection point
      else if (intersections.size() == 1) {
        // get TVector at intersection point
        TVector3 int_tv(intersections.at(0).ToTLorentzVector().Vect());
        length += ( box.Contain(p0) ? (v0 - int_tv).Mag() : (v1 - int_tv).Mag() ); 
      }
      else assert(false); // bad
    }
    // none contained -- either must have zero or two intersections
    if (n_contained == 0) {
      auto intersections = algo.Intersection(line, box);
      if (!(intersections.size() == 0 || intersections.size() == 2)) {
        // more floating point error fixes...
        //
        // figure out which points are near the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        // and which points are near the intersection
        TVector3 vint = intersections.at(0).ToTLorentzVector().Vect();

        bool p0_int = (v0 - vint).Mag() < tol;
        bool p1_int = (v1 - vint).Mag() < tol;
        // exactly one of them should produce the intersection
        assert((p0_int && p0_edge) != (p1_int && p1_edge));
        // void variables when assert-ions are turned off
        (void) p0_int; (void) p1_int;

        // both close to edge -- full length is contained
        if (p0_edge && p1_edge) {
          length += (v0 - v1).Mag();
        }
        // otherwise -- one of them is not on an edge, no length is contained
        else {}
      }
      // assert(intersections.size() == 0 || intersections.size() == 2);
      else if (intersections.size() == 2) {
        TVector3 start(intersections.at(0).ToTLorentzVector().Vect());
        TVector3 end(intersections.at(1).ToTLorentzVector().Vect());
        length += (start - end).Mag();
      }
    }
  }

  return length;
}//ContainedLength

} // namespace analysis

DEFINE_ART_MODULE(analysis::CalorimetryAnalysis)
#endif
