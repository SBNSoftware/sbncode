#include <cmath>
#include <TVector3.h>
#include "Interaction.hh"

namespace util {

double ECCQE(const TVector3& l_momentum, double l_energy) {
  // Based on D. Kaleko, LowEnergyExcess LArLite module ECCQECalculator
  double M_n = 0.939565; // GeV/c^2
  double M_p = 0.938272; // GeV/c^2
  double M_e = 0.000511; // GeV/c^2
  double bindingE = 0.0300; // GeV

  double mp2 = M_p * M_p;
  double me2 = M_e * M_e;
  double mnb = M_n - bindingE;

  TVector3 v(l_momentum);
  double l_mom = sqrt(l_energy * l_energy - me2);
  double l_theta = \
    acos(v.Pz() / sqrt(v.Px()*v.Px() + v.Py()*v.Py() + v.Pz()*v.Pz()));
  double enu_top = mp2 - mnb*mnb - me2 + 2.0 * mnb * l_energy;
  double enu_bot = 2.0 * (mnb - l_energy + l_mom * cos(l_theta));

  return enu_top / enu_bot;
}

double ContainedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geo::BoxBoundedGeo> &boxes) {
  // if points are the same, return 0
  if ((v0 - v1).Mag() < 1e-6) return 0;

  double length = 0;

  // total contained length is sum of lengths in all boxes
  // assuming they are non-overlapping
  for (auto const &box: boxes) {
    int n_contained = box.ContainsPosition(v0) + box.ContainsPosition(v1);
    // both points contained -- length is total length (also can break out of loop)
    if (n_contained == 2) {
      length = (v1 - v0).Mag();
      break;
    }
    // one contained -- have to find intersection point (which must exist)
    else if (n_contained == 1) {
      std::vector<TVector3> intersections = box.GetIntersections(v0, v1); 
      assert(intersections.size() == 1); // should only have one intersection point
      // get the length
      length += ( box.ContainsPosition(v0) ? (v0 - intersections[0]).Mag() : (v1 -  intersections[0]).Mag() );
    }
    else if (n_contained == 0) {
      std::vector<TVector3> intersections = box.GetIntersections(v0, v1); 
      // should have 0 or 2
      if (intersections.size() == 2) {
        length += (intersections[0] - intersections[1]).Mag();
      }
      else assert(intersections.size() == 0);
    }
    else assert(false); // bad
  }

  return length;
}

const simb::MCParticle* MatchMCParticleID(int track_id, const std::vector<simb::MCParticle> &mcparticle_list) {
  for (int i = 0; i < mcparticle_list.size(); i++) {
    if (track_id == mcparticle_list[i].TrackId()) {
      return &mcparticle_list[i];
    }
  }
  return NULL;
}

double MCParticleLength(const simb::MCParticle &particle) {
  double length = 0;
  TVector3 last = particle.Position().Vect();
  for (int i = 1; i < particle.NumberTrajectoryPoints(); i++) {
    TVector3 current = particle.Position(i).Vect();
    length += (current - last).Mag();
    last = current;
  }

  return length;
}

double MCParticleContainedLength(const simb::MCParticle &particle, const std::vector<geo::BoxBoundedGeo> &active_volumes) {
  double contained_length = 0;
  bool outside_AV = true;
  TVector3 start = particle.Position().Vect();
  
  // keep track of the containment position
  int volume_index = -1;
  // now check the active volumes
  for (int i = 0; i < active_volumes.size(); i++) {
    if (outside_AV && active_volumes[i].ContainsPosition(start)) {
      outside_AV = false;
      volume_index = i;
      break;
    }
  }
  
  // if outside AV, set length to 0
  if (outside_AV) return 0.;

  // keep track of the containment volume
  std::vector<geo::BoxBoundedGeo> volume { active_volumes.at(volume_index) };

  TVector3 last = start;
  for (int i = 1; i < particle.NumberTrajectoryPoints(); i++) {
    TVector3 current = particle.Position(i).Vect();
    if (!volume[0].ContainsPosition(current)) {
      break;
    }
    contained_length += ContainedLength(last, current, volume);   
    last = current;
  }

  return contained_length;
}



}  // namespace util

