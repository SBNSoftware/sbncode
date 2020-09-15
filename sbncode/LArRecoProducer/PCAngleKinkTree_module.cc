////////////////////////////////////////////////////////////////////////
// Class:       PCAngleKinkTree
// Plugin Type: analyzer (art v3_05_01)
// File:        PCAngleKinkTree_module.cc
//
// Generated at Tue Sep  8 10:02:52 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_10_00.
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

#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "Products/PCAnglePlane.h"
#include "Products/PCAngleKink.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"



namespace sbn {
  class PCAngleKinkTree;
}


class sbn::PCAngleKinkTree : public art::EDAnalyzer {
public:
  explicit PCAngleKinkTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PCAngleKinkTree(PCAngleKinkTree const&) = delete;
  PCAngleKinkTree(PCAngleKinkTree&&) = delete;
  PCAngleKinkTree& operator=(PCAngleKinkTree const&) = delete;
  PCAngleKinkTree& operator=(PCAngleKinkTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // keep track of file index
  void respondToOpenInputFile(const art::FileBlock& fb) {
    (void) fb;
    fIFile += 1;
  }


private:
  // helpers
  void Clear();
  void FillTruth(const recob::PFParticle &particle, const simb::MCParticle &trueParticle);
  void FillAngles(const recob::PFParticle &particle, const std::vector<art::Ptr<sbn::PCAnglePlane>> &angles);
  void FillMeta(const art::Event &evt, unsigned i_part);
  void FillKinks(const recob::PFParticle &particle, const std::vector<art::Ptr<sbn::PCAngleKink>> &kinks);
  void MatchKinks();
  void MatchPlaneKinks(unsigned Plane, 
      std::vector<int> &scatter_match_plane,
      std::vector<float> &scatter_dist_plane,
      const std::vector<float> &scatterW,
      const std::vector<float> &reco_kinkT, 
      const std::vector<float> &reco_kinkW);

  // config
  std::vector<std::string> fPFParticleTags;
  std::vector<std::string> fAngleTags;
  std::vector<std::string> fKinkTags;
  float fAngleCut;
  bool fOnlyPrimary;

  // output
  TTree *_tree;

  // List of angles per-plane
  std::vector<sbn::PCAngle> fPCAngleU;
  std::vector<sbn::PCAngle> fPCAngleV;
  std::vector<sbn::PCAngle> fPCAngleY;

  std::vector<int> fPCAngleUGen;
  std::vector<int> fPCAngleVGen;
  std::vector<int> fPCAngleYGen;

  std::vector<int> fPCAngleUID;
  std::vector<int> fPCAngleVID;
  std::vector<int> fPCAngleYID;

  std::vector<float> fTrajX;
  std::vector<float> fTrajY;
  std::vector<float> fTrajZ;

  std::vector<float> fTrajPU;
  std::vector<float> fTrajPV;
  std::vector<float> fTrajPY;
  std::vector<float> fTrajPT;

  std::vector<float> fScatterX;
  std::vector<float> fScatterY;
  std::vector<float> fScatterZ;

  std::vector<float> fScatterPU;
  std::vector<float> fScatterPV;
  std::vector<float> fScatterPY;
  std::vector<float> fScatterPT;

  std::vector<float> fScatterMag;
  std::vector<float> fScatterMagU;
  std::vector<float> fScatterMagV;
  std::vector<float> fScatterMagY;

  // matching
  std::vector<int> fScatterMatchU;
  std::vector<float> fScatterMatchDistU;
  std::vector<int> fScatterMatchV;
  std::vector<float> fScatterMatchDistV;
  std::vector<int> fScatterMatchY;
  std::vector<float> fScatterMatchDistY;

  // hit stuff
  std::vector<float> fKinkTimeMaxU;
  std::vector<float> fKinkWireMaxU;
  std::vector<float> fKinkTimeLoU;
  std::vector<float> fKinkWireLoU;
  std::vector<float> fKinkTimeHiU;
  std::vector<float> fKinkWireHiU;
  std::vector<float> fKinkEstAngleU;
  std::vector<float> fKinkMaxAngleU;
  std::vector<float> fKinkLoHiAngleU;
  std::vector<float> fKinkFitAngleU;
  std::vector<float> fKinkFitPitchU;
  std::vector<float> fKinkFitChi2U;

  std::vector<float> fKinkTimeMaxV;
  std::vector<float> fKinkWireMaxV;
  std::vector<float> fKinkTimeLoV;
  std::vector<float> fKinkWireLoV;
  std::vector<float> fKinkTimeHiV;
  std::vector<float> fKinkWireHiV;
  std::vector<float> fKinkEstAngleV;
  std::vector<float> fKinkMaxAngleV;
  std::vector<float> fKinkLoHiAngleV;
  std::vector<float> fKinkFitAngleV;
  std::vector<float> fKinkFitPitchV;
  std::vector<float> fKinkFitChi2V;

  std::vector<float> fKinkTimeMaxY;
  std::vector<float> fKinkWireMaxY;
  std::vector<float> fKinkTimeLoY;
  std::vector<float> fKinkWireLoY;
  std::vector<float> fKinkTimeHiY;
  std::vector<float> fKinkWireHiY;
  std::vector<float> fKinkEstAngleY;
  std::vector<float> fKinkMaxAngleY;
  std::vector<float> fKinkLoHiAngleY;
  std::vector<float> fKinkFitAngleY;
  std::vector<float> fKinkFitPitchY;
  std::vector<float> fKinkFitChi2Y;

  int fTruePDG;
  float fTrueE;
  float fTrueP;

  int fIEvt;
  int fIFile;
  int fEvt;
  int fNo;
};

// static helper functions
float VecAngle(std::array<float, 2> A, std::array<float, 2> B) {
  float costh = (A[0] * B[0] + A[1] * B[1]) \
    / (sqrt(A[0]*A[0] + A[1] * A[1]) * sqrt(B[0]*B[0] + B[1] * B[1]));

  return acos(costh);
}

sbn::PCAngleKinkTree::PCAngleKinkTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},  // ,
    fPFParticleTags(p.get<std::vector<std::string>>("ParticleTags")),
    fAngleTags(p.get<std::vector<std::string>>("AngleTags")),
    fKinkTags(p.get<std::vector<std::string>>("KinkTags")),
    fAngleCut(p.get<float>("AngleCut", 5.)),
    fOnlyPrimary(p.get<bool>("OnlyPrimary", true))
  // More initializers here.
{
  fIEvt = 0;
  fIFile = 0;
  art::ServiceHandle<art::TFileService> tfs;

  _tree = tfs->make<TTree>("PCAngleKinkAnalyzer", "kink_tree");

  _tree->Branch("pcangle_u", &fPCAngleU);
  _tree->Branch("pcangle_v", &fPCAngleV);
  _tree->Branch("pcangle_y", &fPCAngleY);

  _tree->Branch("pcangle_u_id", &fPCAngleUID);
  _tree->Branch("pcangle_v_id", &fPCAngleVID);
  _tree->Branch("pcangle_y_id", &fPCAngleYID);

  _tree->Branch("pcangle_u_gen", &fPCAngleUGen);
  _tree->Branch("pcangle_v_gen", &fPCAngleVGen);
  _tree->Branch("pcangle_y_gen", &fPCAngleYGen);

  _tree->Branch("kink_time_max_u", &fKinkTimeMaxU);
  _tree->Branch("kink_wire_max_u", &fKinkWireMaxU);
  _tree->Branch("kink_time_lo_u", &fKinkTimeLoU);
  _tree->Branch("kink_wire_lo_u", &fKinkWireLoU);
  _tree->Branch("kink_time_hi_u", &fKinkTimeHiU);
  _tree->Branch("kink_wire_hi_u", &fKinkWireHiU);
  _tree->Branch("kink_est_angle_u", &fKinkEstAngleU);
  _tree->Branch("kink_max_angle_u", &fKinkMaxAngleU);
  _tree->Branch("kink_lohi_angle_u", &fKinkLoHiAngleU);
  _tree->Branch("kink_fit_angle_u", &fKinkFitAngleU);
  _tree->Branch("kink_fit_pitch_u", &fKinkFitPitchU);
  _tree->Branch("kink_fit_chi2_u", &fKinkFitChi2U);

  _tree->Branch("kink_time_max_v", &fKinkTimeMaxV);
  _tree->Branch("kink_wire_max_v", &fKinkWireMaxV);
  _tree->Branch("kink_time_lo_v", &fKinkTimeLoV);
  _tree->Branch("kink_wire_lo_v", &fKinkWireLoV);
  _tree->Branch("kink_time_hi_v", &fKinkTimeHiV);
  _tree->Branch("kink_wire_hi_v", &fKinkWireHiV);
  _tree->Branch("kink_est_angle_v", &fKinkEstAngleV);
  _tree->Branch("kink_max_angle_v", &fKinkMaxAngleV);
  _tree->Branch("kink_lohi_angle_v", &fKinkLoHiAngleV);
  _tree->Branch("kink_fit_angle_v", &fKinkFitAngleV);
  _tree->Branch("kink_fit_pitch_v", &fKinkFitPitchV);
  _tree->Branch("kink_fit_chi2_v", &fKinkFitChi2V);

  _tree->Branch("kink_time_max_y", &fKinkTimeMaxY);
  _tree->Branch("kink_wire_max_y", &fKinkWireMaxY);
  _tree->Branch("kink_time_lo_y", &fKinkTimeLoY);
  _tree->Branch("kink_wire_lo_y", &fKinkWireLoY);
  _tree->Branch("kink_time_hi_y", &fKinkTimeHiY);
  _tree->Branch("kink_wire_hi_y", &fKinkWireHiY);
  _tree->Branch("kink_est_angle_y", &fKinkEstAngleY);
  _tree->Branch("kink_max_angle_y", &fKinkMaxAngleY);
  _tree->Branch("kink_lohi_angle_y", &fKinkLoHiAngleY);
  _tree->Branch("kink_fit_angle_y", &fKinkFitAngleY);
  _tree->Branch("kink_fit_pitch_y", &fKinkFitPitchY);
  _tree->Branch("kink_fit_chi2_y", &fKinkFitChi2Y);
  
  _tree->Branch("traj_x", &fTrajX);
  _tree->Branch("traj_y", &fTrajY);
  _tree->Branch("traj_z", &fTrajZ);
  _tree->Branch("traj_pu", &fTrajPU);
  _tree->Branch("traj_pv", &fTrajPV);
  _tree->Branch("traj_py", &fTrajPY);
  _tree->Branch("traj_pt", &fTrajPT);

  _tree->Branch("scatter_x", &fScatterX);
  _tree->Branch("scatter_y", &fScatterY);
  _tree->Branch("scatter_z", &fScatterZ);
  _tree->Branch("scatter_pu", &fScatterPU);
  _tree->Branch("scatter_pv", &fScatterPV);
  _tree->Branch("scatter_py", &fScatterPY);
  _tree->Branch("scatter_pt", &fScatterPT);

  _tree->Branch("scatter_mag", &fScatterMag);
  _tree->Branch("scatter_mag_u", &fScatterMagU);
  _tree->Branch("scatter_mag_v", &fScatterMagV);
  _tree->Branch("scatter_mag_y", &fScatterMagY);

  _tree->Branch("scatter_match_u", &fScatterMatchU);
  _tree->Branch("scatter_match_dist_u", &fScatterMatchDistU);
  _tree->Branch("scatter_match_v", &fScatterMatchV);
  _tree->Branch("scatter_match_dist_v", &fScatterMatchDistV);
  _tree->Branch("scatter_match_y", &fScatterMatchY);
  _tree->Branch("scatter_match_dist_y", &fScatterMatchDistY);

  _tree->Branch("true_pdg", &fTruePDG, "true_pdg/i");
  _tree->Branch("true_p", &fTrueP, "true_pdg/F");
  _tree->Branch("true_e", &fTrueE, "true_pdg/F");

  _tree->Branch("ifile", &fIFile, "ifile/i");
  _tree->Branch("ievt", &fIEvt, "ievt/i");
  _tree->Branch("evt", &fEvt, "evt/i");
  _tree->Branch("no", &fNo, "no/i");
}

void sbn::PCAngleKinkTree::Clear() {
  fPCAngleU.clear();
  fPCAngleV.clear();
  fPCAngleY.clear();

  fPCAngleUID.clear();
  fPCAngleVID.clear();
  fPCAngleYID.clear();

  fPCAngleUGen.clear();
  fPCAngleVGen.clear();
  fPCAngleYGen.clear();

  fKinkTimeMaxU.clear();
  fKinkWireMaxU.clear();
  fKinkTimeLoU.clear();
  fKinkWireLoU.clear();
  fKinkTimeHiU.clear();
  fKinkWireHiU.clear();
  fKinkEstAngleU.clear();
  fKinkMaxAngleU.clear();
  fKinkLoHiAngleU.clear();
  fKinkFitAngleU.clear();
  fKinkFitPitchU.clear();
  fKinkFitChi2U.clear();

  fKinkTimeMaxV.clear();
  fKinkWireMaxV.clear();
  fKinkTimeLoV.clear();
  fKinkWireLoV.clear();
  fKinkTimeHiV.clear();
  fKinkWireHiV.clear();
  fKinkEstAngleV.clear();
  fKinkMaxAngleV.clear();
  fKinkLoHiAngleV.clear();
  fKinkFitAngleV.clear();
  fKinkFitPitchV.clear();
  fKinkFitChi2V.clear();

  fKinkTimeMaxY.clear();
  fKinkWireMaxY.clear();
  fKinkTimeLoY.clear();
  fKinkWireLoY.clear();
  fKinkTimeHiY.clear();
  fKinkWireHiY.clear();
  fKinkEstAngleY.clear();
  fKinkMaxAngleY.clear();
  fKinkLoHiAngleY.clear();
  fKinkFitAngleY.clear();
  fKinkFitPitchY.clear();
  fKinkFitChi2Y.clear();

  fTrajX.clear();
  fTrajY.clear();
  fTrajZ.clear();

  fTrajPU.clear();
  fTrajPV.clear();
  fTrajPY.clear();
  fTrajPT.clear();

  fScatterX.clear();
  fScatterY.clear();
  fScatterZ.clear();

  fScatterMag.clear();
  fScatterMagU.clear();
  fScatterMagV.clear();
  fScatterMagY.clear();

  fScatterMatchU.clear();
  fScatterMatchDistU.clear();
  fScatterMatchV.clear();
  fScatterMatchDistV.clear();
  fScatterMatchY.clear();
  fScatterMatchDistY.clear();

  fScatterPU.clear();
  fScatterPV.clear();
  fScatterPY.clear();
  fScatterPT.clear();

  fTruePDG = 0;
  fTrueE = 0;
  fTrueP = 0;

  fEvt = 0;
  fNo = 0;
}

void sbn::PCAngleKinkTree::FillMeta(const art::Event &evt, unsigned i_part) {
  fNo = i_part;
  fEvt = evt.event();
  fIEvt += 1;
}

void sbn::PCAngleKinkTree::FillTruth(const recob::PFParticle &particle, const simb::MCParticle &trueParticle) {
  (void) particle;

  const geo::GeometryCore *geo = lar::providerFrom<geo::Geometry>();

  // save simple stuff
  fTruePDG = trueParticle.PdgCode();
  fTrueE = trueParticle.Momentum().E();
  fTrueP = trueParticle.Momentum().P();

  // save the particle trajectory
  for (unsigned i_traj = 0; i_traj < trueParticle.NumberTrajectoryPoints(); i_traj ++) {
    TLorentzVector pos = trueParticle.Position(i_traj);

    fTrajX.push_back(pos.X());
    fTrajY.push_back(pos.Y());
    fTrajZ.push_back(pos.Z());

    fTrajPU.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 0, 0, 0) * geo->WirePitch());
    fTrajPV.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 1, 0, 0) * geo->WirePitch());
    fTrajPY.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 2, 0, 0) * geo->WirePitch());
    fTrajPT.push_back(pos.X());

    // also find any elastic scatters
    if (i_traj > 1) {
      // make sure contained in a TPC
      double posarr[3];
      pos.Vect().GetXYZ(posarr);
      geo::TPCID tpc = geo->FindTPCAtPosition(posarr);
      if (tpc.TPC == geo::TPCID::InvalidID) continue;

      // angle in degree
      float angle = trueParticle.Momentum(i_traj).Vect().Angle(trueParticle.Momentum(i_traj-1).Vect());
      if (angle*(180. / M_PI) > fAngleCut) {
	fScatterX.push_back(pos.X());
	fScatterY.push_back(pos.Y());
	fScatterZ.push_back(pos.Z());
	
	fScatterPU.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 0, 0, 0) * geo->WirePitch());
	fScatterPV.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 1, 0, 0) * geo->WirePitch());
	fScatterPY.push_back(geo->WireCoordinate(pos.Y(), pos.Z(), 2, 0, 0) * geo->WirePitch());
	fScatterPT.push_back(pos.X());

        fScatterMag.push_back(angle);

        // project the angle onto each wire-plane
        double scatter_wireplane[3];
        for (unsigned i_plane = 0; i_plane < 3; i_plane++) {
          geo::PlaneID plane(tpc, i_plane);
          double wire_angle_tovert = geo->WireAngleToVertical(geo->View(plane), plane);
          TVector3 v0(trueParticle.Momentum(i_traj-1).Px(), 
            cos(wire_angle_tovert) * trueParticle.Momentum(i_traj-1).Pz() + sin(wire_angle_tovert) * trueParticle.Momentum(i_traj-1).Py(), 0.);
          TVector3 v1(trueParticle.Momentum(i_traj).Px(), 
            cos(wire_angle_tovert) * trueParticle.Momentum(i_traj).Pz() + sin(wire_angle_tovert) * trueParticle.Momentum(i_traj).Py(), 0.);
          if (v0.Mag() < 1e-6 || v1.Mag() < 1e-6) scatter_wireplane[i_plane] = 0.;
          else scatter_wireplane[i_plane] = v1.Unit().Angle(v0.Unit());
        }
        fScatterMagU.push_back(scatter_wireplane[0]);
        fScatterMagV.push_back(scatter_wireplane[1]);
        fScatterMagY.push_back(scatter_wireplane[2]);
      }
    }
  }
}

void sbn::PCAngleKinkTree::FillKinks(const recob::PFParticle &particle, const std::vector<art::Ptr<sbn::PCAngleKink>> &kinks) {
  for (unsigned i_kink = 0; i_kink < kinks.size(); i_kink++) {
    const sbn::PCAngleKink &kink = *kinks[i_kink];
    unsigned plane_no = kink.maxWire.Plane;
    if (plane_no == 0) {
      fKinkTimeMaxU.push_back(kink.position_max[0]);
      fKinkWireMaxU.push_back(kink.position_max[1]);
      fKinkTimeLoU.push_back(kink.position_lo[0]);
      fKinkWireLoU.push_back(kink.position_lo[1]);
      fKinkTimeHiU.push_back(kink.position_hi[0]);
      fKinkWireHiU.push_back(kink.position_hi[1]);
      fKinkEstAngleU.push_back(kink.est_angle);
      fKinkMaxAngleU.push_back(kink.max_angle);
      fKinkLoHiAngleU.push_back(M_PI - VecAngle(kink.vec_lo_at_halfmax_lo, kink.vec_hi_at_halfmax_hi));
      fKinkFitAngleU.push_back(kink.fit_angle);
      fKinkFitPitchU.push_back(kink.fit_pitch);
      fKinkFitChi2U.push_back(kink.fit_chi2);
    }
    else if (plane_no == 1) {
      fKinkTimeMaxV.push_back(kink.position_max[0]);
      fKinkWireMaxV.push_back(kink.position_max[1]);
      fKinkTimeLoV.push_back(kink.position_lo[0]);
      fKinkWireLoV.push_back(kink.position_lo[1]);
      fKinkTimeHiV.push_back(kink.position_hi[0]);
      fKinkWireHiV.push_back(kink.position_hi[1]);
      fKinkEstAngleV.push_back(kink.est_angle);
      fKinkMaxAngleV.push_back(kink.max_angle);
      fKinkLoHiAngleV.push_back(M_PI - VecAngle(kink.vec_lo_at_halfmax_lo, kink.vec_hi_at_halfmax_hi));
      fKinkFitAngleV.push_back(kink.fit_angle);
      fKinkFitPitchV.push_back(kink.fit_pitch);
      fKinkFitChi2V.push_back(kink.fit_chi2);
    }
    else if (plane_no == 2) {
      fKinkTimeMaxY.push_back(kink.position_max[0]);
      fKinkWireMaxY.push_back(kink.position_max[1]);
      fKinkTimeLoY.push_back(kink.position_lo[0]);
      fKinkWireLoY.push_back(kink.position_lo[1]);
      fKinkTimeHiY.push_back(kink.position_hi[0]);
      fKinkWireHiY.push_back(kink.position_hi[1]);
      fKinkEstAngleY.push_back(kink.est_angle);
      fKinkMaxAngleY.push_back(kink.max_angle);
      fKinkLoHiAngleY.push_back(M_PI - VecAngle(kink.vec_lo_at_halfmax_lo, kink.vec_hi_at_halfmax_hi));
      fKinkFitAngleY.push_back(kink.fit_angle);
      fKinkFitPitchY.push_back(kink.fit_pitch);
      fKinkFitChi2Y.push_back(kink.fit_chi2);
    }
  }
}

void sbn::PCAngleKinkTree::FillAngles(const recob::PFParticle &particle, const std::vector<art::Ptr<sbn::PCAnglePlane>> &angles) {
  for (unsigned i_angle = 0; i_angle < angles.size(); i_angle++) {
    const sbn::PCAnglePlane &angle = *angles[i_angle];
    if (angle.plane.Plane == 0) {
      for (unsigned i = 0; i < angle.angles.size(); i++) {
        for (const sbn::PCAngle &a: angle.angles[i]) {
          fPCAngleU.push_back(a);
          fPCAngleUID.push_back(angle.branchIDs[i]);
          fPCAngleUGen.push_back(angle.generations[i]);
        }
      }
    }
    else if (angle.plane.Plane == 1) {
      for (unsigned i = 0; i < angle.angles.size(); i++) {
        for (const sbn::PCAngle &a: angle.angles[i]) {
          fPCAngleV.push_back(a);
          fPCAngleVID.push_back(angle.branchIDs[i]);
          fPCAngleVGen.push_back(angle.generations[i]);
        }
      }
    }
    else if (angle.plane.Plane == 2) {
      for (unsigned i = 0; i < angle.angles.size(); i++) {
        for (const sbn::PCAngle &a: angle.angles[i]) {
          fPCAngleY.push_back(a);
          fPCAngleYID.push_back(angle.branchIDs[i]);
          fPCAngleYGen.push_back(angle.generations[i]);
        }
      }
    }
  }
}

// static helper
void sbn::PCAngleKinkTree::MatchPlaneKinks(unsigned Plane, 
    std::vector<int> &scatter_match_plane,
    std::vector<float> &scatter_dist_plane,
    const std::vector<float> &scatterW,
    const std::vector<float> &reco_kinkT, 
    const std::vector<float> &reco_kinkW) {

  for (unsigned i = 0; i < scatterW.size(); i++) {
    TVector3 true_scatter(fScatterPT[i], scatterW[i], 0.);
    float min_dist = 100000000.;
    int min_ind = -1;
    for (unsigned j = 0; j < reco_kinkT.size(); j++) {
      TVector3 reco_kink(reco_kinkT[j], reco_kinkW[j], 0.);
      if ((true_scatter - reco_kink).Mag() < min_dist) {
        min_dist = (true_scatter - reco_kink).Mag();
        min_ind = j;
      }
    }
    scatter_match_plane.push_back(min_ind);
    scatter_dist_plane.push_back(min_dist);
  }
}

// iterate over all the kinks on each plane and find the closest spatial match
void sbn::PCAngleKinkTree::MatchKinks() {
  MatchPlaneKinks(0, fScatterMatchU, fScatterMatchDistU, fScatterPU, fKinkTimeMaxU, fKinkWireMaxU);
  MatchPlaneKinks(1, fScatterMatchV, fScatterMatchDistV, fScatterPV, fKinkTimeMaxV, fKinkWireMaxV);
  MatchPlaneKinks(2, fScatterMatchY, fScatterMatchDistY, fScatterPY, fKinkTimeMaxY, fKinkWireMaxY);
}

void sbn::PCAngleKinkTree::analyze(art::Event const& evt)
{
  // input data
  std::vector<std::vector<art::Ptr<recob::PFParticle>>> particleList;
  std::vector<art::FindManyP<sbn::PCAnglePlane>> particleAngleList;
  std::vector<art::FindManyP<sbn::PCAngleKink>> particleKinkList;
  std::vector<art::FindManyP<recob::Cluster>> particleClusterList;
  
  for (unsigned i = 0; i < fPFParticleTags.size(); i++) {
    art::Handle<std::vector<recob::PFParticle>> pfparticle_handle;
    evt.getByLabel(fPFParticleTags.at(i), pfparticle_handle);
    particleList.emplace_back();
    art::fill_ptr_vector(particleList.back(), pfparticle_handle);
    particleAngleList.emplace_back(particleList.back(), evt, fAngleTags.at(i));
    particleKinkList.emplace_back(particleList.back(), evt, fKinkTags.at(i));
    particleClusterList.emplace_back(particleList.back(), evt, fPFParticleTags.at(i));
  }

  // flatten
  std::vector<art::Ptr<recob::PFParticle>> particles;
  std::vector<std::vector<art::Ptr<sbn::PCAnglePlane>>> particleAngles;
  std::vector<std::vector<art::Ptr<sbn::PCAngleKink>>> particleKinks;
  std::vector<std::vector<art::Ptr<recob::Cluster>>> particleClusters;
  std::vector<std::vector<art::Ptr<recob::Hit>>> particleHits;

  for (unsigned i = 0; i < particleList.size(); i++) {
    particles.insert(particles.end(), particleList[i].begin(), particleList[i].end());

    for (unsigned j = 0; j < particleList[i].size(); j++) {
      particleAngles.push_back(particleAngleList[i].at(j));
      particleKinks.push_back(particleKinkList[i].at(j));
      particleClusters.push_back(particleClusterList[i].at(j));

      particleHits.emplace_back();
      art::FindManyP<recob::Hit> clusterHits(particleClusters.back(), evt, fPFParticleTags.at(i));
      for (unsigned i_clus = 0; i_clus < particleClusters.back().size(); i_clus++) {
        particleHits.back().insert(particleHits.back().end(), clusterHits.at(i_clus).begin(), clusterHits.at(i_clus).end());
      }
    }
  }

  // also get truth stuff
  art::Handle<std::vector<simb::MCParticle>> trueParticleHandle;
  evt.getByLabel("largeant", trueParticleHandle);
  std::vector<art::Ptr<simb::MCParticle>> trueParticles;
  art::fill_ptr_vector(trueParticles, trueParticleHandle);

  // process
  for (unsigned i_part = 0; i_part < particles.size(); i_part++) {
    const recob::PFParticle &particle = *particles[i_part];
    const std::vector<art::Ptr<sbn::PCAnglePlane>> &angle = particleAngles[i_part];
    const std::vector<art::Ptr<sbn::PCAngleKink>> &kinks = particleKinks[i_part];
    const std::vector<art::Ptr<recob::Hit>> &hits = particleHits[i_part];

    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(hits, true);
    int match_id = matches.size() ? std::max_element(matches.begin(), matches.end(),
                                                     [](const auto &a, const auto &b) { return a.second < b.second; })->first : -1;
    // float matchE = matches.size() ? std::max_element(matches.begin(), matches.end(),
    //                                                  [](const auto &a, const auto &b) { return a.second < b.second; })->second : -1;
    art::Ptr<simb::MCParticle> matched_particle;
    for (art::Ptr<simb::MCParticle> p: trueParticles) {
      if (p->TrackId() == match_id) {
        matched_particle = p;
        break;
      }
    }

    // ignore cases with no match
    if (!matched_particle) continue;

    if (fOnlyPrimary && matched_particle->Process() != "primary") continue;

    // fill stuff
    Clear();
    FillTruth(particle, *matched_particle);
    FillAngles(particle, angle);
    FillKinks(particle, kinks);
    MatchKinks();
    FillMeta(evt, i_part);
    _tree->Fill();
  }
}

DEFINE_ART_MODULE(sbn::PCAngleKinkTree)
