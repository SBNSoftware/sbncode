////////////////////////////////////////////////////////////////////////
///// \file  WireModMuon_module.cc
///// \brief Generates one muon hitting a bounding box, for use by WireMod
/////
///// \author  gputnam@fnal.gov
//////////////////////////////////////////////////////////////////////////

#include "sbncode/EventGenerator/WireModGen/CORSIKAGenSBN.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace sbn {

class WireModMuon: public evgen::CORSIKAGenSBN {
  public:
    void produce(art::Event& evt) override;
    void endRun(art::Run& run) override;
    explicit WireModMuon(fhicl::ParameterSet const& pset);

  private:
    bool GetMuon(const art::Event &e, const simb::MCTruth &intruth, simb::MCParticle &outparticle);
    // Determine if particle ray should be selected
    bool SelectRay(const art::Event &e, const TLorentzVector &pos, const TLorentzVector &mom);
    // Whether a point is in the TPC readout time
    bool PointInTime(const art::Event &e, const geo::TPCID &tpcid, const TVector3 &p, float t);

    struct Plane {
      std::array<geo::Point_t, 4> pts; 

      bool Intersects(const TVector3 &pos, const TVector3 &dir) const;
      bool Intersects(const geo::Point_t &pos, const geo::Vector_t &dir) const;
    };

    // Digested geometry information
    std::vector<Plane> crt_planes;
    std::vector<geo::BoxBoundedGeo> tpc_bounds;

    // Enable for debugging output
    bool fVerbose;

    // Selection configuration
    bool fSelectCXR;
    bool fSelectCRT;
    unsigned fPrescaleCXR;
    unsigned fPrescaleCRT;
    unsigned fNMuonPerEvent;

    // Counters
    unsigned fNCXR;
    unsigned fNCRT;

    unsigned fNSelCXR;
    unsigned fNSelCRT;

    unsigned fNCORSIKA;
};

}

sbn::WireModMuon::WireModMuon(fhicl::ParameterSet const& pset):
  CORSIKAGenSBN(pset)
{

  // Load the CRT planes
  std::vector<double> crt_x0s = pset.get<std::vector<double>>("CRTX0", {});
  std::vector<double> crt_x1s = pset.get<std::vector<double>>("CRTX1", {});
  std::vector<double> crt_y0s = pset.get<std::vector<double>>("CRTY0", {});
  std::vector<double> crt_y1s = pset.get<std::vector<double>>("CRTY1", {});
  std::vector<double> crt_z0s = pset.get<std::vector<double>>("CRTZ0", {});
  std::vector<double> crt_z1s = pset.get<std::vector<double>>("CRTZ1", {});

  for (unsigned i_crt = 0; i_crt < crt_x0s.size(); i_crt++) {
    std::array<geo::Point_t, 4> ps;

    if (crt_x0s[i_crt] == crt_x1s[i_crt] || crt_y0s[i_crt] == crt_y1s[i_crt]) {
      geo::Point_t p1 {crt_x0s[i_crt], crt_y0s[i_crt], crt_z0s[i_crt]};
      geo::Point_t p2 {crt_x1s[i_crt], crt_y1s[i_crt], crt_z0s[i_crt]};
      geo::Point_t p3 {crt_x0s[i_crt], crt_y0s[i_crt], crt_z1s[i_crt]};
      geo::Point_t p4 {crt_x1s[i_crt], crt_y1s[i_crt], crt_z1s[i_crt]};
      ps = {p1, p2, p3, p4};
    }
    else {
      geo::Point_t p1 {crt_x0s[i_crt], crt_y0s[i_crt], crt_z0s[i_crt]};
      geo::Point_t p2 {crt_x0s[i_crt], crt_y1s[i_crt], crt_z0s[i_crt]};
      geo::Point_t p3 {crt_x1s[i_crt], crt_y0s[i_crt], crt_z1s[i_crt]};
      geo::Point_t p4 {crt_x1s[i_crt], crt_y1s[i_crt], crt_z1s[i_crt]};
      ps = {p1, p2, p3, p4};
    }

    sbn::WireModMuon::Plane crt_plane {ps};
    crt_planes.push_back(crt_plane);
  }

  // Load selection config
  fSelectCXR = pset.get<bool>("SelectCXR", true);
  fSelectCRT = pset.get<bool>("SelectCRT", true);

  fPrescaleCXR = pset.get<unsigned>("PrescaleCXR", 1);
  fPrescaleCRT = pset.get<unsigned>("PrescaleCRT", 1);

  fNMuonPerEvent = pset.get<unsigned>("NMuonPerEvent", 1);

  fNCXR = 0;
  fNCRT = 0;

  fNSelCXR = 0;
  fNSelCRT = 0;

  fNCORSIKA = 0;

  // debug output config
  fVerbose = pset.get<bool>("Verbose", true);

}

void sbn::WireModMuon::endRun(art::Run& run) {
  (void) run;

  mf::LogInfo("WireModMuon") << "Selected " << fNSelCXR << " of " << fNCXR << " cathode crossers and " << fNSelCRT << " of " << fNCRT << " CRT plane crossers. Called CORSIKA " << fNCORSIKA << " times.\n";
}

void sbn::WireModMuon::produce(art::Event& evt) {
  // output
  std::unique_ptr<std::vector<simb::MCTruth>> truthcol(new std::vector<simb::MCTruth>);

  simb::MCTruth truth;
  truth.SetOrigin(simb::kCosmicRay);

  for (unsigned i_muon = 0; i_muon < fNMuonPerEvent; i_muon++) {
    // found a muon
    bool foundmuon = false;
    simb::MCParticle particle;

    while (!foundmuon) {
      simb::MCTruth pretruth;
      GetSample(pretruth);
      fNCORSIKA ++;

      foundmuon = GetMuon(evt, pretruth, particle);
    }

    mf::LogInfo("WireModMuon") << "Selected muon at pos(" << particle.Vx() << ", " << particle.Vy() << ", " << particle.Vz() << ") [cm], dir(" 
      << particle.Momentum().Vect().Unit().X() << ", " << particle.Momentum().Vect().Unit().Y() << ", " << particle.Momentum().Vect().Unit().Z() 
      << "), with momentum " << particle.Momentum().Vect().Mag() << " [GeV/c].\n";

    truth.Add(particle);
  }

  // save output
  truthcol->push_back(truth);
  evt.put(std::move(truthcol));

  return;
}

bool sbn::WireModMuon::PointInTime(const art::Event &e, const geo::TPCID &tpcid, const TVector3 &p, float t) {
  auto const &clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const &dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  if (fVerbose) std::cout << "Point: " << p.X() << " " << p.Y() << " " << p.Z() << " is in TPC: " << tpcid << std::endl;

  if (!tpcid) return false;

  geo::PlaneID plane(tpcid, 0); // check against front plane

  // get the tick of the point. Subtract out the trigger time, which gets put into t_tick below
  double tick = dprop.ConvertXToTicks(p.X(), plane);

  // convert the gen time to a TDC
  double t_tdc = clock_data.TPCTick2TDC(clock_data.Time2Tick(t/1e3));

  // add in the offset from the drift
  double x_tick = t_tdc + tick;

  if (fVerbose) std::cout << "Gen time: " << t << " is at TPC tdc: " << t_tdc << ". After drift from X=" << p.X() << " arrival tick is: " << x_tick << ". TPC window is 0 to " << dprop.NumberTimeSamples() << std::endl;

  // this needs to be inside the tpc window
  return (x_tick > 0) && (x_tick < dprop.NumberTimeSamples());
}

bool sbn::WireModMuon::SelectRay(const art::Event &e, const TLorentzVector &pos, const TLorentzVector &mom) {
  // load services
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  bool CXR = false;
  bool CRT = false;

  TVector3 dir = mom.Vect().Unit();
  TVector3 pos3 = pos.Vect();
  const double c_cm_per_ns = 29.9792;
  double vel = mom.Beta()*c_cm_per_ns;
  const double eps = 1e-1; // wiggle room to ensure inside volume

  for (auto const &cryo: geometry->Iterate<geo::CryostatGeo>()) {
    for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryo.ID())) {
      std::vector<TVector3> intersections = TPC.ActiveBoundingBox().GetIntersections(pos3, dir);
      const double check_eps = 1e-1;
      float cathode = (TPC.DriftDir().X() > 0) ? TPC.ActiveBoundingBox().MinX() : TPC.ActiveBoundingBox().MaxX();
      CXR = intersections.size() == 2 && ((abs(intersections[0].X() - cathode) < check_eps) || (abs(intersections[1].X() - cathode) < check_eps));

      // check if cathode point is in-time
      if (CXR) {
        TVector3 cathode_point = (abs(intersections[0].X() - cathode) < abs(intersections[1].X() - cathode)) ? (intersections[0] + dir*eps) : (intersections[1] - dir*eps);

        if (fVerbose) std::cout << "Ray with initial time: " << pos.T() << " and arrival time: " << (pos.T() + cathode_point.Mag() / vel) << " with vel: " << vel << std::endl;
        if (fVerbose) std::cout << "Checking if cathode intersection point: " << cathode_point.X() << " " << cathode_point.Y() << " " << cathode_point.Z() << " is in time." << std::endl;
        CXR = PointInTime(e, TPC.ID(), cathode_point, pos.T() + (cathode_point - pos3).Mag() / vel);
      }
      if (CXR) break;
    }
    if (CXR) break;
  }
  if (fVerbose) std::cout << "CXR: " << CXR << std::endl;

  for (const sbn::WireModMuon::Plane &crt: crt_planes) {
    CRT = crt.Intersects(pos3, dir);
    if (CRT) {
      // require intersect active volume
      //
      // and require either entrance or exit point to be in-time
      for (auto const &cryo: geometry->Iterate<geo::CryostatGeo>()) {
        for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryo.ID())) {
          std::vector<TVector3> intersections = TPC.ActiveBoundingBox().GetIntersections(pos3, dir);
          CRT = intersections.size() == 2;
          if (CRT) {
            if (fVerbose) std::cout << "Ray with initial time: " << pos.T() << " and arrival time: " << (pos.T() + (intersections[0] - pos3).Mag() / vel) << " with vel: " << vel << std::endl;
            CRT = PointInTime(e, TPC.ID(), intersections[0] + dir*eps, pos.T() + (intersections[0] - pos3).Mag() / vel) ||
                  PointInTime(e, TPC.ID(), intersections[1] - dir*eps, pos.T() + (intersections[1] - pos3).Mag() / vel);
          }
          if (CRT) break;
        }
        if (CRT) break;
      }
    }
    if (CRT) break;
  }
  if (fVerbose) std::cout << "CRT: " << CRT << std::endl;

  if (CXR) fNCXR ++;
  if (CRT) fNCRT ++;

  if (CXR && fSelectCXR && (fNCXR % fPrescaleCXR == 0)) {
    fNSelCXR ++;
    return true;
  }

  if (CRT && fSelectCRT && (fNCRT % fPrescaleCRT == 0)) {
    fNSelCRT ++;
    return true;
  }

  return false;

}

bool sbn::WireModMuon::GetMuon(const art::Event &e, const simb::MCTruth &intruth, simb::MCParticle &outparticle) {
  for (int i = 0; i < intruth.NParticles(); ++i) {
    const simb::MCParticle &particle = intruth.GetParticle(i);
    if (abs(particle.PdgCode()) != 13) continue; // only muons

    if (fVerbose) std::cout << "Found muon at pos(" << particle.Vx() << ", " << particle.Vy() << ", " << particle.Vz() << ") [cm], dir(" 
      << particle.Momentum().Vect().Unit().X() << ", " << particle.Momentum().Vect().Unit().Y() << ", " << particle.Momentum().Vect().Unit().Z() 
      << "), with momentum " << particle.Momentum().Vect().Mag() << " [GeV/c].\n";

    if (SelectRay(e, particle.Position(), particle.Momentum())) {
      outparticle = particle;
      return true;
    }
  }

  return false;
}

bool sbn::WireModMuon::Plane::Intersects(const TVector3 &pos, const TVector3 &dir) const {
  geo::Point_t pos_t(pos.X(), pos.Y(), pos.Z());
  geo::Vector_t dir_t(dir.X(), dir.Y(), dir.Z());
  return Intersects(pos_t, dir_t);
}

bool sbn::WireModMuon::Plane::Intersects(const geo::Point_t &pos, const geo::Vector_t &dir) const {
  // Get the u, v vectors
  geo::Vector_t u = (geo::Vector_t(pts[1]) - geo::Vector_t(pts[0])).Unit();
  geo::Vector_t v = (geo::Vector_t(pts[2]) - geo::Vector_t(pts[0])).Unit();
  double ulen = (geo::Vector_t(pts[1]) - geo::Vector_t(pts[0])).Dot(u);
  double vlen = (geo::Vector_t(pts[2]) - geo::Vector_t(pts[0])).Dot(v);

  // calculate the normal of the plane
  geo::Vector_t normal = u.Cross(v);

  double denom = dir.Dot(normal);

  if (abs(denom) < 1e-5) { // parallel to plane
    return false;
  }

  // point along ray where it intersects plane
  double t = (pts[3] - pos).Dot(normal) / denom;

  // intersection happens behind the ray
  if (t < 0) {
    return false;
  }

  geo::Vector_t intersection = pos + dir*t - pts[0];

  double int_u = intersection.Dot(u);
  double int_v = intersection.Dot(v);

  // if (fVerbose) std::cout << "Plane check. ulen: " << ulen << " vlen: " << vlen << ". u: " << int_u << " v: " << int_v << " t: " << t << std::endl;

  return (int_u > 0) && (int_u < ulen) && (int_v > 0) && (int_v < vlen);
}

namespace sbn {
  DEFINE_ART_MODULE(WireModMuon)
}
