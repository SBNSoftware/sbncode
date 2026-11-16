////////////////////////////////////////////////////////////////////////
/// \file   PhotonTruthTools.h
/// \brief  Helpers to decay generator-level pi0s to two photons and to
///         delete one final-state photon from an MCTruth.
///
/// Ported from the MicroBooNE implementation by L. Hagaman / E. Yandel
/// (uboone/ubsim PR #17, ubsim/NGEM), with these changes:
///  - randomness comes from a caller-supplied CLHEP engine (reproducible,
///    seeded through NuRandomService) instead of TRandom3(0);
///  - photons from each decayed pi0 get unique track IDs;
///  - the decayed pi0 keeps its trajectory, Gvtx, rescatter code, weight
///    and polarization, and records the photons as daughters.
////////////////////////////////////////////////////////////////////////

#ifndef SBNCODE_EVENTGENERATOR_SINGLEPHOTON_PHOTONTRUTHTOOLS_H
#define SBNCODE_EVENTGENERATOR_SINGLEPHOTON_PHOTONTRUTHTOOLS_H

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandomEngine.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <algorithm>
#include <utility>
#include <cmath>
#include <vector>

namespace sbn::singlephoton {

  // GENIE status codes (genie::GHepStatus_t)
  constexpr int kStableFinalState = 1;
  constexpr int kDecayedState = 3;

  constexpr int kPdgPhoton = 22;
  constexpr int kPdgPi0 = 111;

  /// Copy everything in `in` except the particle list (origin, generator
  /// info, neutrino) into `out`
  /// (which must be freshly constructed, with its particles already added).
  inline void CopyTruthHeader(simb::MCTruth const& in, simb::MCTruth& out)
  {
    out.SetOrigin(in.Origin());
    simb::MCGeneratorInfo const& gen = in.GeneratorInfo();
    out.SetGeneratorInfo(gen.generator, gen.generatorVersion, gen.generatorConfig);
    if (!in.NeutrinoSet()) return;
    simb::MCNeutrino const& nu = in.GetNeutrino();
    // SetNeutrino rebuilds the MCNeutrino from particle 0 (incoming neutrino)
    // and the outgoing lepton, neither of which is touched by these tools.
    out.SetNeutrino(nu.CCNC(),
                    nu.Mode(),
                    nu.InteractionType(),
                    nu.Target(),
                    nu.HitNuc(),
                    nu.HitQuark(),
                    nu.W(),
                    nu.X(),
                    nu.Y(),
                    nu.QSqr());
  }

  inline bool IsFinalStatePhoton(simb::MCParticle const& p)
  {
    return p.PdgCode() == kPdgPhoton && p.StatusCode() == kStableFinalState;
  }

  /// Whether `p`'s mother (looked up by track ID in `truth`) is a pi0.
  inline bool IsPi0Daughter(simb::MCTruth const& truth, simb::MCParticle const& p)
  {
    for (int j = 0; j < truth.NParticles(); ++j) {
      simb::MCParticle const& m = truth.GetParticle(j);
      if (m.TrackId() == p.Mother()) return m.PdgCode() == kPdgPi0;
    }
    return false;
  }

  /// Copy of `p` with a different status code and extra daughters.
  /// (MCParticle has no status-code setter.)
  inline simb::MCParticle WithStatus(simb::MCParticle const& p,
                                     int status,
                                     std::vector<int> const& extraDaughters = {})
  {
    simb::MCParticle out(p.TrackId(), p.PdgCode(), p.Process(), p.Mother(), p.Mass(), status);
    for (unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i)
      out.AddTrajectoryPoint(p.Position(i), p.Momentum(i));
    for (int i = 0; i < p.NumberDaughters(); ++i)
      out.AddDaughter(p.Daughter(i));
    for (int d : extraDaughters)
      out.AddDaughter(d);
    out.SetGvtx(p.GetGvtx());
    out.SetRescatter(p.Rescatter());
    out.SetWeight(p.Weight());
    out.SetPolarization(p.Polarization());
    out.SetEndProcess(p.EndProcess());
    return out;
  }

  /// Decay every final-state (status 1) pi0 in `in` isotropically to two
  /// photons in its rest frame. The pi0 is kept with status 3 so that LArG4
  /// does not track it; the photons are status 1 with the pi0 as mother.
  /// The Dalitz mode (~1.2%) is not simulated.
  /// \return number of pi0s decayed
  inline unsigned int DecayPi0sToTwoPhotons(simb::MCTruth const& in,
                                            simb::MCTruth& out,
                                            CLHEP::HepRandomEngine& engine)
  {
    CLHEP::RandFlat flat(engine);

    int nextTrackId = 0;
    for (int i = 0; i < in.NParticles(); ++i)
      nextTrackId = std::max(nextTrackId, in.GetParticle(i).TrackId() + 1);

    unsigned int nDecayed = 0;
    for (int i = 0; i < in.NParticles(); ++i) {
      simb::MCParticle const& p = in.GetParticle(i);

      if (p.PdgCode() != kPdgPi0 || p.StatusCode() != kStableFinalState) {
        out.Add(p);
        continue;
      }

      int const id1 = nextTrackId++;
      int const id2 = nextTrackId++;

      TLorentzVector const& x = p.Position();
      TLorentzVector const& P = p.Momentum();
      double const m = P.M(); // use the 4-vector mass so the photons sum exactly to P

      double const cosTheta = flat.fire(-1., 1.);
      double const sinTheta = std::sqrt(std::max(0., 1. - cosTheta * cosTheta));
      double const phi = flat.fire(0., 2. * M_PI);
      double const e = 0.5 * m;
      TVector3 const dir(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);

      TLorentzVector p1(e * dir, e);
      TLorentzVector p2(-e * dir, e);
      TVector3 const boost = P.BoostVector();
      p1.Boost(boost);
      p2.Boost(boost);

      out.Add(WithStatus(p, kDecayedState, {id1, id2}));

      for (auto const& [id, mom] : {std::pair{id1, p1}, std::pair{id2, p2}}) {
        simb::MCParticle gamma(id, kPdgPhoton, "primary", p.TrackId(), 0.0, kStableFinalState);
        gamma.AddTrajectoryPoint(x, mom);
        gamma.SetGvtx(p.GetGvtx());
        out.Add(gamma);
      }
      ++nDecayed;
    }

    CopyTruthHeader(in, out);
    return nDecayed;
  }

  /// Remove one final-state photon, chosen uniformly at random, from `in`.
  /// If `onlyFromPi0` is set, only photons whose mother is a pi0 in the same
  /// record are candidates. If there are no candidates `out` is a copy of `in`.
  /// \return index (in `in`) of the deleted particle, or -1 if none
  inline int DeleteOneRandomPhoton(simb::MCTruth const& in,
                                   simb::MCTruth& out,
                                   CLHEP::HepRandomEngine& engine,
                                   bool onlyFromPi0 = false)
  {
    std::vector<int> candidates;
    for (int i = 0; i < in.NParticles(); ++i) {
      simb::MCParticle const& p = in.GetParticle(i);
      if (!IsFinalStatePhoton(p)) continue;
      if (onlyFromPi0 && !IsPi0Daughter(in, p)) continue;
      candidates.push_back(i);
    }

    int deleted = -1;
    if (!candidates.empty()) {
      CLHEP::RandFlat flat(engine);
      deleted = candidates[flat.fireInt(static_cast<long>(candidates.size()))];
    }

    for (int i = 0; i < in.NParticles(); ++i)
      if (i != deleted) out.Add(in.GetParticle(i));

    CopyTruthHeader(in, out);
    return deleted;
  }

  /// Number of final-state photons; with `onlyFromPi0`, only those whose
  /// mother is a pi0 (so e.g. nuclear de-excitation photons are not counted).
  inline unsigned int CountFinalStatePhotons(simb::MCTruth const& truth, bool onlyFromPi0 = false)
  {
    unsigned int n = 0;
    for (int i = 0; i < truth.NParticles(); ++i) {
      simb::MCParticle const& p = truth.GetParticle(i);
      if (IsFinalStatePhoton(p) && (!onlyFromPi0 || IsPi0Daughter(truth, p))) ++n;
    }
    return n;
  }

} // namespace sbn::singlephoton

#endif // SBNCODE_EVENTGENERATOR_SINGLEPHOTON_PHOTONTRUTHTOOLS_H
