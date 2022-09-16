/**
 * @file        NuVertexChargeTree_module.cc
 * @author      Gray Putnam (grayputnam@uchicago.edu)
 */

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

#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "../LArRecoProducer/LArReco/TrackMomentumCalculator.h"

#include "sbncode/TPCReco/VertexStub/StubMergeAlgorithms.h"

#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbnobj/Common/Reco/Stub.h"

namespace sbn {
  class NuVertexChargeTree;
}

/**
 * @brief Analyzer module for use with sbn::Stub and sbn::VertexHit objects.
 * 
 * Intended for use with neutrino interactions and neutrino-like particlegun
 * events. Not intended for use on data or on MC with cosmics.
 */
class sbn::NuVertexChargeTree : public art::EDAnalyzer {
public:
  explicit NuVertexChargeTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuVertexChargeTree(NuVertexChargeTree const&) = delete;
  NuVertexChargeTree(NuVertexChargeTree&&) = delete;
  NuVertexChargeTree& operator=(NuVertexChargeTree const&) = delete;
  NuVertexChargeTree& operator=(NuVertexChargeTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // keep track of file index
  void respondToOpenInputFile(const art::FileBlock& fb) override {
    (void) fb;
    fIFile += 1;
  }

private:
  // helpers
  void Clear();
  bool InFV(TVector3 pos);
  void FillMeta(const art::Event &evt);
  void FillNeutrino(const simb::MCTruth &nu, 
   const recob::Vertex &vert,
   const std::vector<art::Ptr<recob::PFParticle>> &slice_pfps,
   const std::vector<std::vector<std::pair<int, float>>> &slice_pfp_matches,
   const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles, 
   const std::vector<const sim::GeneratedParticleInfo *> &infos,
   const geo::GeometryCore &geo,
   const detinfo::DetectorPropertiesData &dprop);

  void FillVertexHits(
    const std::vector<art::Ptr<sbn::VertexHit>> &hits, 
    const std::vector<art::Ptr<recob::Hit>> &hit_hits,
    const detinfo::DetectorPropertiesData &dprop,
    const detinfo::DetectorClocksData &clock_data,
    const std::vector<art::Ptr<simb::MCParticle>> &trueParticles);

  void FillStubs(
    const std::vector<sbn::StubInfo> &stubs,
    const std::vector<unsigned> &stub_hit_inds,
    const geo::GeometryCore *geo,
    const detinfo::DetectorPropertiesData &dprop,
    const detinfo::DetectorClocksData &clock_data,
    const std::vector<art::Ptr<simb::MCParticle>> &trueParticles);

  void FillStubPFPs(
    const std::vector<art::Ptr<recob::PFParticle>> &stub_pfps,
    const std::vector<std::vector<art::Ptr<recob::Hit>>> &stub_pfp_hits,
    const detinfo::DetectorClocksData &clock_data,
    const std::vector<art::Ptr<simb::MCParticle>> &trueParticles);

  void FillSlice(const std::vector<art::Ptr<recob::PFParticle>> &slice_particles);

  // config
  std::vector<std::string> fPFParticleTags;
  std::vector<std::string> fVertexHitTags;
  std::vector<std::string> fStubTags;
  std::array<float, 6> fFiducialInset;
  std::vector<geo::BoxBoundedGeo> fFiducialVolumes;
  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  trkf::TrackMomentumCalculator fRangeCalculator;
  bool fCorrectSCE;

  // output
  TTree *_tree;

  int fIEvt;
  int fIFile;
  int fEvt;

  float fNuE;
  int fNuMode;
  int fNuIntType;

  std::vector<float> fFSPPx;
  std::vector<float> fFSPPy;
  std::vector<float> fFSPPz;

  std::vector<float> fFSPE;
  std::vector<float> fFSPQ;
  std::vector<float> fFSPKE;
  std::vector<float> fFSPVisE;
  std::vector<float> fFSPLen;
  std::vector<int> fFSPPID;
  std::vector<bool> fFSPMatched;
  std::vector<int> fFSPMatchedPID;
  std::vector<bool> fFSPMatchedPFPPrimary;
  std::vector<bool> fFSPMtachedPrimary;
  std::vector<bool> fFSPIsStopping;

  std::vector<float> fFSPStartX;
  std::vector<float> fFSPStartY;
  std::vector<float> fFSPStartZ;
  std::vector<float> fFSPEndX;
  std::vector<float> fFSPEndY;
  std::vector<float> fFSPEndZ;

  std::vector<float> fFSPMaxQ0;
  std::vector<float> fFSPMaxQ1;
  std::vector<float> fFSPMaxQ2;
  std::vector<float> fFSPPitch0;
  std::vector<float> fFSPPitch1;
  std::vector<float> fFSPPitch2;

  float fRecoVertDist;

  std::vector<int> fVertexHitSPID;
  std::vector<float> fVertexHitPitch;
  std::vector<float> fVertexHitdQdx;
  std::vector<float> fVertexHitdEdx;
  std::vector<float> fVertexHitCharge;
  std::vector<float> fVertexHitDist;
  std::vector<float> fVertexHitVtxW;
  std::vector<int> fVertexHitPlane;
  std::vector<int> fVertexHitWire;
  std::vector<float> fVertexHitPeakX;
  std::vector<float> fVertexHitTrueEnergy;
  std::vector<float> fVertexHitTrueProtonEnergy;
  std::vector<float> fVertexHitTrueBraggProtonEnergy;
  std::vector<float> fVertexHitTrueX;
  std::vector<float> fVertexHitTrueY;
  std::vector<float> fVertexHitTrueZ;

  std::vector<float> fStubEndx;
  std::vector<float> fStubEndy;
  std::vector<float> fStubEndz;
  std::vector<float> fStubDirx;
  std::vector<float> fStubDiry;
  std::vector<float> fStubDirz;
  std::vector<float> fStubLength;
  std::vector<float> fStubCharge;
  std::vector<float> fStubPitch;
  std::vector<float> fStubTrkPitch;
  std::vector<float> fStubLenP;
  std::vector<int> fStubNWire;
  std::vector<int> fStubHitInd;
  std::vector<int> fStubPFPID;
  std::vector<int> fStubNPlane;
  std::vector<float> fStubVtxW;
  std::vector<int> fStubEndW;
  std::vector<float> fStubVtxEField;
  std::vector<float> fStubEndEField;

  std::vector<int> fStubQLength;
  std::vector<float> fStubQs;
  std::vector<int> fStubOnTracks;
  std::vector<int> fStubWires;

  std::vector<float> fStubTrueKE;
  std::vector<float> fStubTrueLength;
  std::vector<int> fStubTrueID;

  std::vector<int> fStubTruePdg;
  std::vector<float> fStubTruePx;
  std::vector<float> fStubTruePy;
  std::vector<float> fStubTruePz;

  std::vector<int> fStubPFPTruePdg;
  std::vector<float> fStubPFPTruePx;
  std::vector<float> fStubPFPTruePy;
  std::vector<float> fStubPFPTruePz;

  std::vector<int> fStubMatchLength;
  std::vector<int> fStubMatch;
  std::vector<int> fStubMatchOverlaps;
  std::vector<float> fStubMatchDot;

  std::vector<int> fStubXPlaneMatchLength;
  std::vector<int> fStubXPlaneMatch;
  std::vector<float> fStubXPlaneMatchTOff;
  std::vector<float> fStubXPlaneMatchQOff;
  std::vector<float> fStubXPlaneMatchdQdxOff;
  std::vector<float> fStubXPlaneMatchPeakQOff;

  int fNSliceParticles;
  int fNSlicePrimaryParticles;
  int fNSlicePrimaryTracks;
  int fNSlicePrimaryNus;
  int fNSlicePrimaryShowers;
};

sbn::NuVertexChargeTree::NuVertexChargeTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},  // ,
    fPFParticleTags(p.get<std::vector<std::string>>("ParticleTags")),
    fVertexHitTags(p.get<std::vector<std::string>>("VertexHitTags")),
    fStubTags(p.get<std::vector<std::string>>("StubTags")),
    fFiducialInset({p.get<float>("xmin"), p.get<float>("xmax"), p.get<float>("ymin"), p.get<float>("ymax"), p.get<float>("zmin"), p.get<float>("zmax")}), 
    fRangeCalculator(p.get<float>("MinTrackLength", 0.1)),
    fCorrectSCE(p.get<bool>("CorrectSCE"))
{
  fIEvt = 0;
  fIFile = 0;
  art::ServiceHandle<art::TFileService> tfs;

  _tree = tfs->make<TTree>("VertexEAnalyzer", "kink_tree");

  _tree->Branch("ifile", &fIFile, "ifile/i");
  _tree->Branch("ievt", &fIEvt, "ievt/i");
  _tree->Branch("evt", &fEvt, "evt/i");

  _tree->Branch("nue", &fNuE, "nue/F");
  _tree->Branch("numode", &fNuMode, "numode/i");
  _tree->Branch("nuinttype", &fNuIntType, "nuinttype/i");

  _tree->Branch("fsp_px", &fFSPPx);
  _tree->Branch("fsp_py", &fFSPPy);
  _tree->Branch("fsp_pz", &fFSPPz);

  _tree->Branch("fsp_start_x", &fFSPStartZ);
  _tree->Branch("fsp_start_y", &fFSPStartY);
  _tree->Branch("fsp_start_z", &fFSPStartZ);
  _tree->Branch("fsp_end_x", &fFSPEndZ);
  _tree->Branch("fsp_end_y", &fFSPEndY);
  _tree->Branch("fsp_end_z", &fFSPEndZ);

  _tree->Branch("fsp_e", &fFSPE);
  _tree->Branch("fsp_ke", &fFSPKE);
  _tree->Branch("fsp_vise", &fFSPVisE);
  _tree->Branch("fsp_q", &fFSPQ);
  _tree->Branch("fsp_len", &fFSPLen);
  _tree->Branch("fsp_pid", &fFSPPID);
  _tree->Branch("fsp_matched", &fFSPMatched);
  _tree->Branch("fsp_matched_pfp_primary", &fFSPMatchedPFPPrimary);
  _tree->Branch("fsp_matched_pid", &fFSPMatchedPID);
  _tree->Branch("fsp_matched_primary", &fFSPMtachedPrimary);
  _tree->Branch("fsp_process_is_topping", &fFSPIsStopping);

  _tree->Branch("fsp_maxq0", &fFSPMaxQ0);
  _tree->Branch("fsp_maxq1", &fFSPMaxQ1);
  _tree->Branch("fsp_maxq2", &fFSPMaxQ2);
  _tree->Branch("fsp_pitch0", &fFSPPitch0);
  _tree->Branch("fsp_pitch1", &fFSPPitch1);
  _tree->Branch("fsp_pitch2", &fFSPPitch2);

  _tree->Branch("reco_vert_dist", &fRecoVertDist, "reco_vert_dist/F");

  _tree->Branch("vertex_hit_spid", &fVertexHitSPID);
  _tree->Branch("vertex_hit_pitch", &fVertexHitPitch);
  _tree->Branch("vertex_hit_dqdx", &fVertexHitdQdx);
  _tree->Branch("vertex_hit_dedx", &fVertexHitdEdx);
  _tree->Branch("vertex_hit_charge", &fVertexHitCharge);
  _tree->Branch("vertex_hit_dist", &fVertexHitDist);
  _tree->Branch("vertex_hit_vtxw", &fVertexHitVtxW);
  _tree->Branch("vertex_hit_plane", &fVertexHitPlane);
  _tree->Branch("vertex_hit_wire", &fVertexHitWire);
  _tree->Branch("vertex_hit_peak_x", &fVertexHitPeakX);
  _tree->Branch("vertex_hit_true_e", &fVertexHitTrueEnergy);
  _tree->Branch("vertex_hit_true_proton_e", &fVertexHitTrueProtonEnergy);
  _tree->Branch("vertex_hit_true_bragg_proton_e", &fVertexHitTrueBraggProtonEnergy);
  _tree->Branch("vertex_hit_true_x", &fVertexHitTrueX);
  _tree->Branch("vertex_hit_true_y", &fVertexHitTrueY);
  _tree->Branch("vertex_hit_true_z", &fVertexHitTrueZ);

  _tree->Branch("stub_endx", &fStubEndx);
  _tree->Branch("stub_endy", &fStubEndy);
  _tree->Branch("stub_endz", &fStubEndz);

  _tree->Branch("stub_dirx", &fStubDirx);
  _tree->Branch("stub_diry", &fStubDiry);
  _tree->Branch("stub_dirz", &fStubDirz);
  _tree->Branch("stub_length", &fStubLength);
  _tree->Branch("stub_lenp", &fStubLenP);
  _tree->Branch("stub_charge", &fStubCharge);
  _tree->Branch("stub_trk_pitch", &fStubTrkPitch);
  _tree->Branch("stub_pitch", &fStubPitch);
  _tree->Branch("stub_nwire", &fStubNWire);
  _tree->Branch("stub_hit_ind", &fStubHitInd);
  _tree->Branch("stub_pfpid", &fStubPFPID);
  _tree->Branch("stub_nplane", &fStubNPlane);
  _tree->Branch("stub_vtx_w", &fStubVtxW);
  _tree->Branch("stub_end_w", &fStubEndW);
  _tree->Branch("stub_vtx_efield", &fStubVtxEField);
  _tree->Branch("stub_end_efield", &fStubEndEField);

  _tree->Branch("stub_qlength", &fStubQLength);
  _tree->Branch("stub_qs", &fStubQs);
  _tree->Branch("stub_ontracks", &fStubOnTracks);
  _tree->Branch("stub_wires", &fStubWires);

  _tree->Branch("stub_true_ke", &fStubTrueKE);
  _tree->Branch("stub_true_length", &fStubTrueLength);
  _tree->Branch("stub_true_id", &fStubTrueID);

  _tree->Branch("stub_true_pdg", &fStubTruePdg);
  _tree->Branch("stub_true_px", &fStubTruePx);
  _tree->Branch("stub_true_py", &fStubTruePy);
  _tree->Branch("stub_true_pz", &fStubTruePz);

  _tree->Branch("stub_pfp_true_pdg", &fStubPFPTruePdg);
  _tree->Branch("stub_pfp_true_px", &fStubPFPTruePx);
  _tree->Branch("stub_pfp_true_py", &fStubPFPTruePy);
  _tree->Branch("stub_pfp_true_pz", &fStubPFPTruePz);

  _tree->Branch("stub_match_length", &fStubMatchLength);
  _tree->Branch("stub_match", &fStubMatch);
  _tree->Branch("stub_match_overlaps", &fStubMatchOverlaps);
  _tree->Branch("stub_match_dot", &fStubMatchDot);

  _tree->Branch("stub_xplane_match_length", &fStubXPlaneMatchLength);
  _tree->Branch("stub_xplane_match", &fStubXPlaneMatch);
  _tree->Branch("stub_xplane_match_toff", &fStubXPlaneMatchTOff);
  _tree->Branch("stub_xplane_match_qoff", &fStubXPlaneMatchQOff);
  _tree->Branch("stub_xplane_match_dqdxoff", &fStubXPlaneMatchdQdxOff);
  _tree->Branch("stub_xplane_match_peakqoff", &fStubXPlaneMatchPeakQOff);

  _tree->Branch("n_slice_particles", &fNSliceParticles, "n_slice_particles/i");
  _tree->Branch("n_slice_primary_particles", &fNSlicePrimaryParticles, "n_slice_primary_particles/i");
  _tree->Branch("n_slice_primary_tracks", &fNSlicePrimaryTracks, "n_slice_primary_tracks/i");
  _tree->Branch("n_slice_primary_nus", &fNSlicePrimaryNus, "n_slice_primary_nus/i");
  _tree->Branch("n_slice_primary_showers", &fNSlicePrimaryShowers, "n_slice_primary_showers/i");

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

    // fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
    fFiducialVolumes.emplace_back(XMin + fFiducialInset[0], XMax - fFiducialInset[1], 
                                  YMin + fFiducialInset[2], YMax - fFiducialInset[3], 
                                  ZMin + fFiducialInset[4], ZMax - fFiducialInset[5]);
  }
}

bool sbn::NuVertexChargeTree::InFV(TVector3 pos) {
  for (auto const &FV: fFiducialVolumes) {
    if (FV.ContainsPosition(pos)) return true;
  }
  return false;
}

const simb::MCParticle *Genie2G4MCParticle(
  const simb::MCParticle &genie_part,
  const simb::MCTruth &mctruth,
  const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles,
  const std::vector<const sim::GeneratedParticleInfo *> infos) {

  const simb::MCParticle *ret = nullptr;
  for (int iparticle = 0; iparticle < (int)g4_mcparticles.size(); iparticle++) {
    if (infos[iparticle]->hasGeneratedParticleIndex() &&
        (int)infos[iparticle]->generatedParticleIndex() < mctruth.NParticles() && // TODO: why is this number sometimes bigger than the number of particles?
        mctruth.GetParticle(infos[iparticle]->generatedParticleIndex()).TrackId() == genie_part.TrackId() &&
        g4_mcparticles[iparticle]->Process() == "primary" /* TODO: will have to remove this restriction to include secondary particles*/) {

      // if a genie particle re-scatters in g4 and makes more particles, then multiple g4 particles can match to a 
      // genie particle. Thus, we also check that the start location of the associated genie particle matches the g4
      // and that the pdgid matches (to be on the safe side)
      //
      // Note that this should be accounted for by requiring the process to be primary. This is a bit of redundancy.
      const simb::MCParticle& matched_genie_particle = mctruth.GetParticle(infos[iparticle]->generatedParticleIndex());
      if ((matched_genie_particle.Position().Vect() - g4_mcparticles[iparticle]->Position().Vect()).Mag() < 1e-4 &&
        matched_genie_particle.PdgCode() == g4_mcparticles[iparticle]->PdgCode()) {

        // this should only be true for one particle
        assert(ret == nullptr);
        ret = g4_mcparticles[iparticle].get();
      }
    }
  }
  return ret;
}

void sbn::NuVertexChargeTree::Clear() {
  fEvt = 0;

  fNuE = 0;
  fNuMode = -1;
  fNuIntType = -1;

  fFSPPx.clear();
  fFSPPy.clear();
  fFSPPz.clear();

  fFSPStartX.clear();
  fFSPStartY.clear();
  fFSPStartZ.clear();
  fFSPEndX.clear();
  fFSPEndY.clear();
  fFSPEndZ.clear();

  fFSPQ.clear();
  fFSPE.clear();
  fFSPKE.clear();
  fFSPVisE.clear();
  fFSPLen.clear();
  fFSPPID.clear();
  fFSPMatched.clear();
  fFSPMatchedPFPPrimary.clear();
  fFSPMatchedPID.clear();
  fFSPMtachedPrimary.clear();
  fFSPIsStopping.clear();

  fFSPMaxQ0.clear();
  fFSPMaxQ1.clear();
  fFSPMaxQ2.clear();
  fFSPPitch0.clear();
  fFSPPitch1.clear();
  fFSPPitch2.clear();

  fRecoVertDist = 0.;
  fVertexHitSPID.clear();
  fVertexHitPitch.clear();
  fVertexHitdQdx.clear();
  fVertexHitdEdx.clear();
  fVertexHitCharge.clear();
  fVertexHitDist.clear();
  fVertexHitVtxW.clear();
  fVertexHitPlane.clear();
  fVertexHitWire.clear();
  fVertexHitPeakX.clear();
  fVertexHitTrueEnergy.clear();
  fVertexHitTrueProtonEnergy.clear();
  fVertexHitTrueBraggProtonEnergy.clear();
  fVertexHitTrueX.clear();
  fVertexHitTrueY.clear();
  fVertexHitTrueZ.clear();

  fStubEndx.clear();
  fStubEndy.clear();
  fStubEndz.clear();
  fStubDirx.clear();
  fStubDiry.clear();
  fStubDirz.clear();
  fStubLength.clear();
  fStubLenP.clear();
  fStubCharge.clear();
  fStubTrkPitch.clear();
  fStubPitch.clear();
  fStubNWire.clear();
  fStubHitInd.clear();
  fStubPFPID.clear();
  fStubNPlane.clear();
  fStubVtxW.clear();
  fStubEndW.clear();
  fStubVtxEField.clear();
  fStubEndEField.clear();

  fStubQLength.clear();
  fStubQs.clear();
  fStubOnTracks.clear();
  fStubWires.clear();

  fStubTrueKE.clear();
  fStubTrueLength.clear();
  fStubTrueID.clear();
  fStubTruePdg.clear();
  fStubTruePx.clear();
  fStubTruePy.clear();
  fStubTruePz.clear();
  fStubPFPTruePdg.clear();
  fStubPFPTruePx.clear();
  fStubPFPTruePy.clear();
  fStubPFPTruePz.clear();

  fStubMatchLength.clear();
  fStubMatch.clear();
  fStubMatchOverlaps.clear();
  fStubMatchDot.clear();

  fStubXPlaneMatchLength.clear();
  fStubXPlaneMatch.clear();
  fStubXPlaneMatchTOff.clear();
  fStubXPlaneMatchQOff.clear();
  fStubXPlaneMatchdQdxOff.clear();
  fStubXPlaneMatchPeakQOff.clear();

  fNSliceParticles = 0;
  fNSlicePrimaryParticles = 0;
  fNSlicePrimaryTracks = 0;
  fNSlicePrimaryNus = 0;
  fNSlicePrimaryShowers = 0;
}

void sbn::NuVertexChargeTree::FillMeta(const art::Event &evt) {
  fEvt = evt.event();
  fIEvt += 1;
}

void sbn::NuVertexChargeTree::FillSlice(const std::vector<art::Ptr<recob::PFParticle>> &slice_particles) {
  fNSliceParticles = slice_particles.size();
  for (unsigned i = 0; i < slice_particles.size(); i++) {
    bool is_primary = slice_particles[i]->IsPrimary();
    if (!is_primary) {
      int parent = slice_particles[i]->Parent();
      for (unsigned j = 0; j < slice_particles.size(); j++) {
        if ((int)slice_particles[j]->Self() == parent) {
          is_primary = slice_particles[j]->IsPrimary() && lar_pandora::LArPandoraHelper::IsNeutrino(slice_particles[j]);
          break;
        }
      }
    }
    if (is_primary) {
      fNSlicePrimaryParticles ++;
    }

    if (is_primary && lar_pandora::LArPandoraHelper::IsTrack(slice_particles[i])) {
      fNSlicePrimaryTracks ++;
    }
    else if (is_primary && lar_pandora::LArPandoraHelper::IsShower(slice_particles[i])) {
      fNSlicePrimaryShowers ++;
    }
    else if (is_primary && lar_pandora::LArPandoraHelper::IsNeutrino(slice_particles[i])) {
      fNSlicePrimaryNus ++;
    }
   
  }
}

void sbn::NuVertexChargeTree::FillStubPFPs(
    const std::vector<art::Ptr<recob::PFParticle>> &stub_pfps,
    const std::vector<std::vector<art::Ptr<recob::Hit>>> &stub_pfp_hits,
    const detinfo::DetectorClocksData &clock_data,
    const std::vector<art::Ptr<simb::MCParticle>> &trueParticles) {

  for (unsigned i_stub = 0; i_stub < stub_pfps.size(); i_stub++) {
    // fill PFP info
    if (stub_pfps[i_stub]) {
      fStubPFPID.push_back(stub_pfps[i_stub]->Self());
    }
    else {
      fStubPFPID.push_back(-1);
    }

    // truth-matching
    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, stub_pfp_hits[i_stub]);
    // sort
    std::sort(matches.begin(), matches.end(), 
      [](auto const &lhs, auto const &rhs) {
        return lhs.second > rhs.second;
    });

    // get the matched particle
    art::Ptr<simb::MCParticle> G4;
    int matched_track_id = matches.size() ? abs(matches[0].first) : -1;
    for (unsigned i_g4 = 0; i_g4 < trueParticles.size(); i_g4++) {
      if (matched_track_id == trueParticles[i_g4]->TrackId()) {
        G4 = trueParticles[i_g4];
        break;
      }
    }

    if (G4) {
      fStubPFPTruePdg.push_back(G4->PdgCode());
      fStubPFPTruePx.push_back(G4->Momentum().Px());
      fStubPFPTruePy.push_back(G4->Momentum().Py());
      fStubPFPTruePz.push_back(G4->Momentum().Pz());
    }
    else {
      fStubPFPTruePdg.push_back(-1);
      fStubPFPTruePx.push_back(-1);
      fStubPFPTruePy.push_back(-1);
      fStubPFPTruePz.push_back(-1);
    }

  } // end iterate over Stub-PFParticle's
}

void sbn::NuVertexChargeTree::FillStubs(
    const std::vector<sbn::StubInfo> &stubs,
    const std::vector<unsigned> &stub_hit_inds,
    const geo::GeometryCore *geo,
    const detinfo::DetectorPropertiesData &dprop,
    const detinfo::DetectorClocksData &clock_data,
    const std::vector<art::Ptr<simb::MCParticle>> &trueParticles) {

  // TODO: fix -- for now, use a null space-charge service
  const spacecharge::SpaceCharge *sce = lar::providerFrom<spacecharge::SpaceChargeService>();
      
  for (unsigned i_stub = 0; i_stub < stubs.size(); i_stub++) {
    const sbn::Stub &stub = stubs[i_stub].stub;

    // fill reco info
    fStubLength.push_back((stub.vtx - stub.end).r());
    fStubEndx.push_back(stub.end.X());
    fStubEndy.push_back(stub.end.Y());
    fStubEndz.push_back(stub.end.Z());

    fStubDirx.push_back((stub.end - stub.vtx).Unit().X());
    fStubDiry.push_back((stub.end - stub.vtx).Unit().Y());
    fStubDirz.push_back((stub.end - stub.vtx).Unit().Z());
    fStubNWire.push_back(stub.CoreNHit());
    fStubTrkPitch.push_back(stub.trkpitch.front());
    fStubPitch.push_back(stub.pitch.front());
    fStubLenP.push_back(fRangeCalculator.GetTrackMomentum((stub.vtx - stub.end).r(), 2212));

    // Only count charge on the main stub
    fStubCharge.push_back(stub.CoreCharge());

    // Compute per-charge stuff
    // Save per-charge stuff
    fStubQLength.push_back(stub.hits.front().size());
    for (const sbn::StubHit &h: stub.hits.front()) {
      fStubQs.push_back(h.charge);
      fStubOnTracks.push_back(h.ontrack);
      fStubWires.push_back(h.wire);
    }

    fStubHitInd.push_back(stub_hit_inds[i_stub]);
    fStubNPlane.push_back(stub.plane.size());
    fStubVtxW.push_back(stub.vtx_w.front());
    fStubEndW.push_back(stub.hit_w.front());

    fStubVtxEField.push_back(stub.efield_vtx);
    fStubEndEField.push_back(stub.efield_end);

    // truth-matching
    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, stubs[i_stub].hits);
    // sort
    std::sort(matches.begin(), matches.end(), 
      [](auto const &lhs, auto const &rhs) {
        return lhs.second > rhs.second;
    });

    // get the matched particle
    art::Ptr<simb::MCParticle> G4;
    int matched_track_id = matches.size() ? abs(matches[0].first) : -1;
    for (unsigned i_g4 = 0; i_g4 < trueParticles.size(); i_g4++) {
      if (matched_track_id == trueParticles[i_g4]->TrackId()) {
        G4 = trueParticles[i_g4];
        break;
      }
    }

    // Fill in truth information
    if (G4) {
      fStubTruePdg.push_back(G4->PdgCode());
      fStubTrueID.push_back(G4->TrackId());
      fStubTrueKE.push_back(G4->Momentum().E() - G4->Momentum().M());
      fStubTrueLength.push_back((G4->Position().Vect() - G4->EndPosition().Vect()).Mag());

      fStubTruePx.push_back(G4->Momentum().Px());
      fStubTruePy.push_back(G4->Momentum().Py());
      fStubTruePz.push_back(G4->Momentum().Pz());
    }
    else {
      fStubTruePdg.push_back(-1);
      fStubTrueID.push_back(-1);
      fStubTrueKE.push_back(-1);
      fStubTrueLength.push_back(-1);

      fStubTruePx.push_back(-1);
      fStubTruePy.push_back(-1);
      fStubTruePz.push_back(-1);
    }

    // intra plane matching
    {
      unsigned start_size = fStubMatch.size();
      for (unsigned j_stub = 0; j_stub < stubs.size(); j_stub++) {
        if (i_stub == j_stub) continue;
        if (stubs[i_stub].stub.plane.size() != 1 || stubs[j_stub].stub.plane.size() != 1) continue; // Ignore multi-plane stubs
        if (stubs[i_stub].stub.plane.front() == stubs[j_stub].stub.plane.front()) {
          fStubMatch.push_back(j_stub);
          fStubMatchOverlaps.push_back(sbn::StubContains(stubs[i_stub], stubs[j_stub]));
          fStubMatchDot.push_back(sbn::StubDirectionDot(stubs[i_stub], stubs[j_stub], geo, dprop)); 
        }
      }
      fStubMatchLength.push_back(fStubMatch.size() - start_size);
    }

    // inter plane matching 
    {
      unsigned start_size = fStubXPlaneMatch.size();
      for (unsigned j_stub = 0; j_stub < stubs.size(); j_stub++) {
        if (i_stub == j_stub) continue;
        if (stubs[i_stub].stub.plane.size() != 1 || stubs[j_stub].stub.plane.size() != 1) continue; // Ignore multi-plane stubs
        if (stubs[i_stub].stub.plane.front().Plane != stubs[j_stub].stub.plane.front().Plane) {
          fStubXPlaneMatch.push_back(j_stub);
          fStubXPlaneMatchTOff.push_back(sbn::StubTimeOffset(stubs[i_stub], stubs[j_stub], clock_data, dprop));
          fStubXPlaneMatchQOff.push_back(sbn::StubChargeOffset(stubs[i_stub], stubs[j_stub]));
          fStubXPlaneMatchPeakQOff.push_back(sbn::StubPeakChargeOffset(stubs[i_stub], stubs[j_stub]));
          fStubXPlaneMatchdQdxOff.push_back((stubs[i_stub].stub.plane.front().TPC == stubs[j_stub].stub.plane.front().TPC) ?
            sbn::StubPeakdQdxOffset(stubs[i_stub], stubs[j_stub], geo, sce, dprop) : std::numeric_limits<float>::signaling_NaN());
        }
      }
      fStubXPlaneMatchLength.push_back(fStubXPlaneMatch.size() - start_size);
    }

  } // end iterate over vertex-stubs
}
void sbn::NuVertexChargeTree::FillVertexHits(
    const std::vector<art::Ptr<sbn::VertexHit>> &hits, 
    const std::vector<art::Ptr<recob::Hit>> &hit_hits,
    const detinfo::DetectorPropertiesData &dprop,
    const detinfo::DetectorClocksData &clock_data,
    const std::vector<art::Ptr<simb::MCParticle>> &trueParticles) {

  for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
    const sbn::VertexHit &hit = *hits[i_hit];
    unsigned plane = hit.wire.Plane;
    fVertexHitPlane.push_back(plane);
    fVertexHitWire.push_back(hit.wire.Wire);
    fVertexHitPeakX.push_back(dprop.ConvertTicksToX(hit_hits[i_hit]->PeakTime(), hit_hits[i_hit]->WireID()));
    fVertexHitCharge.push_back(hit.charge);
    fVertexHitDist.push_back(hit.proj_dist_to_vertex);
    fVertexHitVtxW.push_back(hit.vtxw);
    fVertexHitdQdx.push_back(hit.dqdx);
    fVertexHitSPID.push_back(hit.spID);
    fVertexHitPitch.push_back(hit.pitch);
    fVertexHitdEdx.push_back(hit.dedx);

    fVertexHitTrueEnergy.push_back(0.);
    fVertexHitTrueProtonEnergy.push_back(0.);
    fVertexHitTrueBraggProtonEnergy.push_back(0.);

    fVertexHitTrueX.push_back(0.);
    fVertexHitTrueY.push_back(0.);
    fVertexHitTrueZ.push_back(0.);

    // Get the matching energy of the hit 
    art::ServiceHandle<cheat::BackTrackerService> bt;
    std::vector<sim::IDE> hit_ides;
    try {
      std::vector<sim::IDE> try_get_hit_ides = bt->HitToAvgSimIDEs(clock_data, *hit_hits[i_hit]);
      hit_ides = try_get_hit_ides;
    }
    catch(...) {}

    for (unsigned i_ide = 0; i_ide < hit_ides.size(); i_ide++) {
      fVertexHitTrueEnergy.back() += hit_ides[i_ide].energy;
      fVertexHitTrueX.back() += hit_ides[i_ide].energy * hit_ides[i_ide].x;
      fVertexHitTrueY.back() += hit_ides[i_ide].energy * hit_ides[i_ide].y;
      fVertexHitTrueZ.back() += hit_ides[i_ide].energy * hit_ides[i_ide].z;

      for (unsigned i_g4 = 0; i_g4 < trueParticles.size(); i_g4++) {
        if (abs(hit_ides[i_ide].trackID) == trueParticles[i_g4]->TrackId()) {
          if (abs(trueParticles[i_g4]->PdgCode()) == 2212) {
            fVertexHitTrueProtonEnergy.back() += hit_ides[i_ide].energy;
            TVector3 ide_pos(hit_ides[i_ide].x, hit_ides[i_ide].y, hit_ides[i_ide].z);
            if ((trueParticles[i_g4]->EndPosition().Vect() - ide_pos).Mag() < 0.5) {
              fVertexHitTrueBraggProtonEnergy.back() += hit_ides[i_ide].energy;
            }
          }
          break;
        }
      }
    }

    if (hit_ides.size()) {
      fVertexHitTrueX.back() /= fVertexHitTrueEnergy.back();
      fVertexHitTrueY.back() /= fVertexHitTrueEnergy.back();
      fVertexHitTrueZ.back() /= fVertexHitTrueEnergy.back();
    }

  } // end iterate over hits
}

void sbn::NuVertexChargeTree::FillNeutrino(const simb::MCTruth &nu, 
   const recob::Vertex &vert,
   const std::vector<art::Ptr<recob::PFParticle>> &slice_pfps,
   const std::vector<std::vector<std::pair<int, float>>> &slice_pfp_matches,
   const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles, 
   const std::vector<const sim::GeneratedParticleInfo *> &infos,
   const geo::GeometryCore &geo,
   const detinfo::DetectorPropertiesData &dprop) {

  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  // TODO: fix -- for now, use a null space-charge service
  const spacecharge::SpaceCharge *sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  // If we are not correcting for SCE, then blank out the service
  if (!fCorrectSCE) sce = nullptr;

  TVector3 vpos(vert.position().X(), vert.position().Y(), vert.position().Z());

  if (nu.NeutrinoSet()) {
    fRecoVertDist = (nu.GetNeutrino().Nu().Position().Vect() - vpos).Mag();

    fNuE = nu.GetNeutrino().Nu().E();
    fNuMode = nu.GetNeutrino().Mode();
    fNuIntType = nu.GetNeutrino().InteractionType();
  }
  else {
    fRecoVertDist = (nu.GetParticle(0).Position().Vect() - vpos).Mag();

    fNuE = -1;
    fNuMode = -1;
    fNuIntType = -1;
  }

  // save each stable final-state particle
  for (int i_part = 0; i_part < nu.NParticles(); i_part++) {
    const simb::MCParticle &particle = nu.GetParticle(i_part);
    // is stable?
    if (particle.StatusCode() != 1) continue;

    // std::cout << "Stable FSP -- PdgID: " << particle.PdgCode() << " E: " << (particle.Momentum().E() - particle.Momentum().M()) << std::endl;

    // Find the corresponding G4 particle
    const simb::MCParticle *G4 = Genie2G4MCParticle(particle, nu, g4_mcparticles, infos);

    if (!G4) continue;
    // std::cout << "Tracked by G4!\n";

    // lookup its hits

    // save some crap
    fFSPPx.push_back(particle.Momentum().Px());
    fFSPPy.push_back(particle.Momentum().Py());
    fFSPPz.push_back(particle.Momentum().Pz());

    fFSPStartX.push_back(G4->Position().X());
    fFSPStartY.push_back(G4->Position().Y());
    fFSPStartZ.push_back(G4->Position().Z());
    fFSPEndX.push_back(G4->EndPosition().X());
    fFSPEndY.push_back(G4->EndPosition().Y());
    fFSPEndZ.push_back(G4->EndPosition().Z());

    fFSPE.push_back(particle.Momentum().E());
    fFSPKE.push_back(particle.Momentum().E() - particle.Momentum().M());
    fFSPPID.push_back(particle.PdgCode());

    fFSPLen.push_back((G4->Position().Vect() - G4->EndPosition().Vect()).Mag());

    bool is_stopping =  (G4->EndProcess() == "Decay" ||
                                        G4->EndProcess() == "CoupledTransportation" ||
                                        G4->EndProcess() == "FastScintillation" ||
                                        G4->EndProcess() == "muMinusCaptureAtRest" ||
                                        G4->EndProcess() == "LArVoxelReadoutScoringProcess");
    fFSPIsStopping.push_back(is_stopping);

    // see if there is a match
    //
    // Total up the energy fromthis particle
    float thisVisE = 0.;
    float thisQ = 0.;

    std::map<geo::WireID, float> chargemap;
    for (geo::View_t view: geo.Views()) {
      std::vector<const sim::IDE*> particle_ides(backtracker->TrackIdToSimIDEs_Ps(G4->TrackId(), view));
      for (const sim::IDE *ide: particle_ides) {
        thisVisE += ide->energy;

        // Correct for electron lifetime 
        geo::Point_t p {ide->x, ide->y, ide->z};
        geo::TPCGeo const* tpc = geo.PositionToTPCptr(p);
        if (tpc) {
          float driftT = tpc->FirstPlane().DistanceFromPlane(p) / dprop.DriftVelocity();
          thisQ += ide->numElectrons * exp(driftT / dprop.ElectronLifetime());

          geo::PlaneGeo const& plane = tpc->Plane(view);
          try {
            geo::WireID w = plane.NearestWireID(p);
            chargemap[w] += ide->numElectrons * exp(driftT / dprop.ElectronLifetime());
          }
          catch(geo::InvalidWireError const&) {} 
        }
      }
    }

    fFSPVisE.push_back(thisVisE / 3. /* average over each plane */);
    fFSPQ.push_back(thisQ / 3.);

    for (unsigned p = 0; p < 3; p++) {
      float maxQ = -1;
      for (auto const &pair: chargemap) {
        if (pair.first.Plane == p) {
          if (pair.second > maxQ) maxQ = pair.second;
        }
      }

      if (p == 0) fFSPMaxQ0.push_back(maxQ);
      if (p == 1) fFSPMaxQ1.push_back(maxQ);
      if (p == 2) fFSPMaxQ2.push_back(maxQ);
    }

    geo::Point_t start_loc(particle.Position().Vect());
    geo::Vector_t start_dir(particle.Momentum().Vect().Unit());
    geo::TPCID firstTPC = geo.FindTPCAtPosition(start_loc);
    if (firstTPC) {
      fFSPPitch0.push_back(sbn::GetPitch(&geo, sce, start_loc, start_dir, geo.View(geo::PlaneID(firstTPC, 0)), firstTPC, true, true));
      fFSPPitch1.push_back(sbn::GetPitch(&geo, sce, start_loc, start_dir, geo.View(geo::PlaneID(firstTPC, 1)), firstTPC, true, true));
      fFSPPitch2.push_back(sbn::GetPitch(&geo, sce, start_loc, start_dir, geo.View(geo::PlaneID(firstTPC, 2)), firstTPC, true, true));
    }
    else {
      fFSPPitch0.push_back(-1);
      fFSPPitch1.push_back(-1);
      fFSPPitch2.push_back(-1);
    }

    bool found_match = false;
    bool primary_match = false;
    int pfp_pdgid = -1;
    bool is_primary = false;
    for (unsigned i_pfp = 0; i_pfp < slice_pfp_matches.size(); i_pfp++) {
      for (unsigned i_match = 0; i_match < slice_pfp_matches[i_pfp].size(); i_match++) {

        if (slice_pfp_matches[i_pfp][i_match].first == G4->TrackId() &&
            (slice_pfp_matches[i_pfp][i_match].second / thisVisE > 0.5)) {
          found_match = true;     
          primary_match = i_match == 0;
          pfp_pdgid = slice_pfps[i_pfp]->PdgCode(); 

          is_primary = slice_pfps[i_pfp]->IsPrimary();
          if (!is_primary) {
            int parent = slice_pfps[i_pfp]->Parent();
            for (unsigned j_pfp = 0; j_pfp < slice_pfps.size(); j_pfp++) {
              if ((int)slice_pfps[j_pfp]->Self() == parent) {
                is_primary = slice_pfps[j_pfp]->IsPrimary() && lar_pandora::LArPandoraHelper::IsNeutrino(slice_pfps[j_pfp]);
                break;
              }
            }
          }

          // std::cout << "Found match! For particle ID: " << particle.PdgCode() << " gen ID: " << particle.TrackId() << std::endl;
          // std::cout << "G4 particleID: " << G4->PdgCode() << " track ID: " <<  G4->TrackId() << " vis E: " << thisVisE << std::endl;
          // std::cout << "Match info\n";
          // for (unsigned i_match = 0; i_match < slice_pfp_matches[i_pfp].size(); i_match++) {
          //   std::cout << "Match ID: " << slice_pfp_matches[i_pfp][i_match].first << " E: " << slice_pfp_matches[i_pfp][i_match].second << std::endl; 
          // }

          break;
        }
      }
      if (found_match) break;
    }

    fFSPMatched.push_back(found_match);
    fFSPMatchedPFPPrimary.push_back(is_primary);
    fFSPMatchedPID.push_back(pfp_pdgid);
    fFSPMtachedPrimary.push_back(primary_match);

  }
}

void sbn::NuVertexChargeTree::analyze(art::Event const& evt) {
  // services
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> inventory_service;
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clock_data);
  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();

  // input data
  std::vector<std::vector<art::Ptr<recob::Slice>>> sliceList;
  std::vector<art::FindManyP<recob::Hit>> sliceHitList;
  std::vector<art::FindManyP<recob::PFParticle>> slicePFParticleList;
  std::vector<art::FindManyP<sbn::VertexHit>> sliceVertexHitList;

  for (unsigned i = 0; i < fPFParticleTags.size(); i++) {
    art::Handle<std::vector<recob::Slice>> slice_handle;
    evt.getByLabel(fPFParticleTags.at(i), slice_handle);
    sliceList.emplace_back();
    art::fill_ptr_vector(sliceList.back(), slice_handle);
    sliceHitList.emplace_back(sliceList.back(), evt, fPFParticleTags.at(i));
    slicePFParticleList.emplace_back(sliceList.back(), evt, fPFParticleTags.at(i));
    sliceVertexHitList.emplace_back(sliceList.back(), evt, fVertexHitTags.at(i));
  }

  // flatten
  std::vector<art::Ptr<recob::Slice>> slices;
  std::vector<std::vector<art::Ptr<recob::Hit>>> sliceHits;

  std::vector<std::vector<art::Ptr<recob::PFParticle>>> sliceParticles;
  std::vector<std::vector<std::vector<art::Ptr<recob::Hit>>>> sliceParticleHits;
  std::vector<std::vector<std::vector<art::Ptr<recob::Vertex>>>> sliceParticleVertexs;
  
  std::vector<std::vector<art::Ptr<sbn::VertexHit>>> sliceVertexHits;
  std::vector<std::vector<art::Ptr<recob::Vertex>>> sliceVertexHitVtxs;
  std::vector<std::vector<art::Ptr<recob::Hit>>> sliceVertexHitHits;
  std::vector<std::vector<art::Ptr<recob::SpacePoint>>> sliceVertexHitSPs;
  std::vector<std::vector<art::Ptr<sbn::Stub>>> sliceStubs;
  std::vector<std::vector<unsigned>> sliceStubHitInds;
  std::vector<std::vector<std::vector<art::Ptr<recob::Hit>>>> sliceStubHits;
  std::vector<std::vector<art::Ptr<recob::PFParticle>>> sliceStubPFPs;
  std::vector<std::vector<std::vector<art::Ptr<recob::Hit>>>> sliceStubPFPHits;

  for (unsigned i = 0; i < sliceList.size(); i++) {
    slices.insert(slices.end(), sliceList[i].begin(), sliceList[i].end());
    for (unsigned j = 0; j < sliceList[i].size(); j++) {
      sliceHits.push_back(sliceHitList[i].at(j));
      sliceParticles.push_back(slicePFParticleList[i].at(j));

      sliceParticleHits.emplace_back();
      art::FindManyP<recob::Cluster> sliceParticleClusters(sliceParticles.back(), evt, fPFParticleTags.at(i));
      for (unsigned i_pfp = 0; i_pfp < sliceParticles.back().size(); i_pfp++) {
        sliceParticleHits.back().emplace_back();
        const std::vector<art::Ptr<recob::Cluster>> &thisParticleClusters = sliceParticleClusters.at(i_pfp);
        art::FindManyP<recob::Hit> clusterHits(thisParticleClusters, evt, fPFParticleTags.at(i));
        for (unsigned i_clus = 0; i_clus < thisParticleClusters.size(); i_clus++) {
          sliceParticleHits.back().back().insert(sliceParticleHits.back().back().end(), clusterHits.at(i_clus).begin(), clusterHits.at(i_clus).end());
        }
      }

      sliceParticleVertexs.emplace_back();
      art::FindManyP<recob::Vertex> findSliceParticleVertexs(sliceParticles.back(), evt, fPFParticleTags.at(i));
      for (unsigned i_pfp = 0; i_pfp < sliceParticles.back().size(); i_pfp++) {
        sliceParticleVertexs.back().push_back(findSliceParticleVertexs.at(i_pfp));
      }

      sliceVertexHits.push_back(sliceVertexHitList[i].at(j));

      art::FindManyP<recob::Vertex> vertexHitVtxs(sliceVertexHits.back(), evt, fVertexHitTags.at(i));
      sliceVertexHitVtxs.emplace_back();
      for (unsigned i_vhit = 0; i_vhit < sliceVertexHits.back().size(); i_vhit++) {
        sliceVertexHitVtxs.back().push_back(vertexHitVtxs.at(i_vhit).at(0));
      }

      art::FindManyP<recob::Hit> vertexHitHits(sliceVertexHits.back(), evt, fVertexHitTags.at(i));
      sliceVertexHitHits.emplace_back();
      for (unsigned i_vhit = 0; i_vhit < sliceVertexHits.back().size(); i_vhit++) {
        sliceVertexHitHits.back().push_back(vertexHitHits.at(i_vhit).at(0));
      }

      art::Ptr<recob::SpacePoint> nullSP;
      art::FindManyP<recob::SpacePoint> vertexHitSPs(sliceVertexHitHits.back(), evt, fPFParticleTags.at(i));
      sliceVertexHitSPs.emplace_back();
      for (unsigned i_sp = 0; i_sp < sliceVertexHitHits.back().size(); i_sp++) {
        const std::vector<art::Ptr<recob::SpacePoint>> &this_sp = vertexHitSPs.at(i_sp);
        if (this_sp.size()) sliceVertexHitSPs.back().push_back(this_sp.at(0));
        else sliceVertexHitSPs.back().push_back(nullSP);
      }

      // art::Ptr<sbn::Stub> nullVH;
      art::FindManyP<sbn::Stub> vertexStubs(sliceVertexHits.back(), evt, fStubTags.at(i));
      sliceStubs.emplace_back();
      sliceStubHitInds.emplace_back();
      for (unsigned i_stub = 0; i_stub < sliceVertexHits.back().size(); i_stub++) {
        const std::vector<art::Ptr<sbn::Stub>> &this_vh = vertexStubs.at(i_stub);
        for (unsigned i_vh_stub = 0; i_vh_stub < this_vh.size(); i_vh_stub++) {
          sliceStubs.back().push_back(this_vh.at(i_vh_stub));
          sliceStubHitInds.back().push_back(i_stub);
        }
      }

      art::FindManyP<recob::Hit> vertexStubHits(sliceStubs.back(), evt, fStubTags.at(i));
      sliceStubHits.emplace_back();
      for (unsigned i_stub = 0; i_stub < sliceStubs.back().size(); i_stub++) {
        sliceStubHits.back().push_back(vertexStubHits.at(i_stub));
      }

      art::FindManyP<recob::PFParticle> vertexStubPFPs(sliceStubs.back(), evt, fStubTags.at(i));
      art::Ptr<recob::PFParticle> nullPFP;
      sliceStubPFPs.emplace_back();
      for (unsigned i_stub = 0; i_stub < sliceStubs.back().size(); i_stub++) {
        const std::vector<art::Ptr<recob::PFParticle>> &thisStubPFP = vertexStubPFPs.at(i_stub);
        if (thisStubPFP.size()) sliceStubPFPs.back().push_back(thisStubPFP.at(0));
        else sliceStubPFPs.back().push_back(nullPFP);
      }

      sliceStubPFPHits.emplace_back();
      for (unsigned i_stub = 0; i_stub < sliceStubPFPs.back().size(); i_stub++) {
        bool found = false;
        for (unsigned i_pfp = 0; i_pfp < sliceParticles.back().size(); i_pfp++) {
          if (sliceStubPFPs.back()[i_stub] == sliceParticles.back()[i_pfp]) {
            sliceStubPFPHits.back().push_back(sliceParticleHits.back()[i_pfp]);
            found = true;
            break;
          }
        }
        if (!found) {
          sliceStubPFPHits.back().emplace_back();
        }
      }

    } // end iterate over slices
  } // end iterate over pfp-tags

  // also get truth stuff
  art::Handle<std::vector<simb::MCTruth>> trueNuHandle;
  evt.getByLabel("generator", trueNuHandle);
  std::vector<art::Ptr<simb::MCTruth>> trueNus;
  art::fill_ptr_vector(trueNus, trueNuHandle);

  art::Handle<std::vector<simb::MCParticle>> trueParticleHandle;
  evt.getByLabel("largeant", trueParticleHandle);
  std::vector<art::Ptr<simb::MCParticle>> trueParticles;
  art::fill_ptr_vector(trueParticles, trueParticleHandle);

  // get associations from MCTruth information to truth 
  art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> truth_to_particles(trueNus, evt, "largeant");

  // for (unsigned i_g4 = 0; i_g4 < trueParticles.size(); i_g4++) {
  //   const simb::MCParticle &particle = *trueParticles[i_g4];
  //   std::cout << "particle ID: " << particle.TrackId() << " pdg: " << particle.PdgCode() << " start E: " << (particle.E() - particle.Momentum().M()) << " start X: " << particle.Vx() << std::endl;
  // }

  // iterate over each neutrino
  for (unsigned i_nu = 0; i_nu < trueNus.size(); i_nu++) {
    const simb::MCTruth &nu = *trueNus[i_nu];
    // require inside fiducial volume
    TVector3 true_vertex = nu.NeutrinoSet() ? nu.GetNeutrino().Nu().Position().Vect() : nu.GetParticle(0).Position().Vect();
    if (!InFV(true_vertex)) continue;

    // try to do reco-matching
    //
    // Total up the visible energy of this neutrino
    float visE = 0.;
    for (unsigned i_g4 = 0; i_g4 < trueParticles.size(); i_g4++) {
      const simb::MCParticle &particle = *trueParticles[i_g4];
      art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(particle.TrackId());
      if (truth == trueNus[i_nu]) {
        std::vector<const sim::IDE*> particle_ides(backtracker->TrackIdToSimIDEs_Ps(particle.TrackId()));
        for (const sim::IDE *ide: particle_ides) visE += ide->energy;
      }
    }

    int slice_match = -1;
    // see if the visible energy from any slice matches
    for (unsigned i_slice = 0; i_slice < slices.size(); i_slice++) {
      const recob::Slice &slice = *slices[i_slice];
      (void) slice;
      float matchingE = 0.;
      float totalE = 0.;
      std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, sliceHits.at(i_slice));
      for (unsigned i_match = 0; i_match < matches.size(); i_match++) {
        totalE += matches[i_match].second;
      }

      for (unsigned i_g4 = 0; i_g4 < trueParticles.size(); i_g4++) {
        for (unsigned i_match = 0; i_match < matches.size(); i_match++) {
          if (matches[i_match].first == trueParticles[i_g4]->TrackId()) {
            art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(trueParticles[i_g4]->TrackId());
            if (truth == trueNus[i_nu]) {
              matchingE += matches[i_match].second;
            }
            break;
          }
        }
      }

      // std::cout << "Slice: " << i_slice << " match E: " << matchingE << " vis E: " << visE << " total E: " << totalE << std::endl;

      if (matchingE / visE > 0.5) {
        if (matchingE / totalE > 0.5) slice_match = i_slice;
        break;
      }
    }

    if (slice_match < 0) continue;
    // We found a matching slice! 

    // Do reco stuff
    // Collect reco stuff
    const std::vector<art::Ptr<recob::PFParticle>> thisSliceParticles = sliceParticles.at(slice_match);
    const std::vector<std::vector<art::Ptr<recob::Hit>>> thisSliceParticleHits = sliceParticleHits.at(slice_match);
    const std::vector<std::vector<art::Ptr<recob::Vertex>>> thisSliceParticleVertexs = sliceParticleVertexs.at(slice_match);

    const std::vector<art::Ptr<sbn::VertexHit>> &thisVertexHits = sliceVertexHits.at(slice_match);
    const std::vector<art::Ptr<recob::Hit>> &thisVertexHitHits = sliceVertexHitHits.at(slice_match);
    const std::vector<art::Ptr<recob::SpacePoint>> &thisVertexHitSPs = sliceVertexHitSPs.at(slice_match);
    const std::vector<art::Ptr<recob::Vertex>> &thisVertexHitVtxs = sliceVertexHitVtxs.at(slice_match);

    const std::vector<art::Ptr<sbn::Stub>> &thisStubs = sliceStubs.at(slice_match);
    const std::vector<std::vector<art::Ptr<recob::Hit>>> &thisStubHits = sliceStubHits.at(slice_match);
    const std::vector<unsigned> &thisStubHitInds = sliceStubHitInds.at(slice_match);
    const std::vector<art::Ptr<recob::PFParticle>> &thisStubPFPs = sliceStubPFPs.at(slice_match);
    const std::vector<std::vector<art::Ptr<recob::Hit>>> &thisStubPFPHits = sliceStubPFPHits.at(slice_match);

    // Build the stub infos
    std::vector<sbn::StubInfo> thisStubInfos;
    for (unsigned i = 0; i < thisStubs.size(); i++) {
      thisStubInfos.emplace_back();

      sbn::StubInfo& newStub = thisStubInfos.back();

      newStub.stub = *thisStubs[i];
      newStub.pfp = thisStubPFPs.at(i);
      newStub.hits = thisStubHits.at(i);
      newStub.vhit = thisVertexHits.at(thisStubHitInds[i]);
      newStub.vhit_hit = thisVertexHitHits.at(thisStubHitInds[i]); 
    }

    (void) thisVertexHitSPs;
    (void) thisVertexHitVtxs;
    //for (unsigned i_hit = 0; i_hit < thisVertexHitHits.size(); i_hit++) {
    //  std::cout << "Hit: " << i_hit << std::endl;
    //  art::Ptr<recob::SpacePoint> sp = thisVertexHitSPs.at(i_hit);
    //  if (sp) {
    //    std::cout << "SP! " << sp->ID() << " X: " << sp->XYZ()[0] << " Y: " << sp->XYZ()[1] << " Z: " << sp->XYZ()[2] << std::endl;
    //  }
    //}

    // Get the primary vertex
    art::Ptr<recob::Vertex> vert;
    for (unsigned i_pfp = 0; i_pfp < thisSliceParticles.size(); i_pfp++) {
      if (thisSliceParticles[i_pfp]->IsPrimary() &&
          /* TODO: remove. Needed because of issue in SCE correction module */ thisSliceParticleVertexs.at(i_pfp).size()) {
        vert = thisSliceParticleVertexs.at(i_pfp).at(0);
        break;
      }
    }

    double badpos[3] {-999., -999., -999.};
    recob::Vertex badvert(badpos);
    const recob::Vertex &nu_vert = (vert) ? *vert : badvert;

    // Lookup truth-matching
    std::vector<std::vector<std::pair<int, float>>> thisSlicePFPMatches; 
    for (unsigned i_pfp = 0; i_pfp < thisSliceParticles.size(); i_pfp++) {
      std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(clock_data, thisSliceParticleHits.at(i_pfp));
      // sort
      std::sort(matches.begin(), matches.end(), 
        [](auto const &lhs, auto const &rhs) {
          return lhs.second > rhs.second;
      });
      thisSlicePFPMatches.push_back(matches);
    }

    Clear();
    FillMeta(evt);

    // save the true-neutrino info
    FillNeutrino(nu, nu_vert, thisSliceParticles, thisSlicePFPMatches, truth_to_particles.at(i_nu), truth_to_particles.data(i_nu), *geo, dprop);

    // The vertex hits!
    FillVertexHits(thisVertexHits, thisVertexHitHits, dprop, clock_data, trueParticles);

    FillStubs(thisStubInfos, thisStubHitInds, geo, dprop, clock_data, trueParticles);

    FillStubPFPs(thisStubPFPs, thisStubPFPHits, clock_data, trueParticles);

    FillSlice(thisSliceParticles);

    _tree->Fill();
  }
}

DEFINE_ART_MODULE(sbn::NuVertexChargeTree)
