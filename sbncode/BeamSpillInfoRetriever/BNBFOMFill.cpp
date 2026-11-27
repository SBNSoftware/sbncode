/**
 * @file   sbncode/BeamSpillInfoRetriever/BNBFOMFill.cpp
 * @brief  BNB figure of merit for spills with missing inputs (see BNBFOMFill.h).
 */
#include "sbncode/BeamSpillInfoRetriever/BNBFOMFill.h"
#include "sbncode/BeamSpillInfoRetriever/getFOM.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace sbn {

  namespace {

    constexpr double MissingValue = -999.;
    bool isValid(double v) { return std::isfinite(v) && v != MissingValue; }

    double spillTime(BNBSpillInfo const& s) { return s.spill_time_s + 1e-9 * s.spill_time_ns; }

    bool measuredWidth(BNBBeamState const& s) { return s.widthSource != BNBBeamState::Nominal; }

    /// Sets the FOM fields of a spill that was filled or got a new width.
    void storeFOM(BNBSpillInfo& spill, BNBBeamState const& state, double fom) {
      if (measuredWidth(state)) {
        spill.FOM = fom;
        spill.PreFitFOM = MissingValue;
        spill.NoMultiWireFOM = computeFOM(state, true);
      }
      else {
        spill.FOM = spill.PreFitFOM = MissingValue;
        spill.NoMultiWireFOM = fom;
      }
    }

  } // local namespace


  std::vector<unsigned int> improveBNBqualityFOMs
    (std::vector<BNBSpillInfo>& spills, BNBFOMFillConfig const& cfg)
  {
    using namespace fomstatus;
    std::size_t const n = spills.size();

    // spills in time order
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
      [&spills](std::size_t a, std::size_t b){ return spillTime(spills[a]) < spillTime(spills[b]); });

    std::vector<BNBBeamState> S(n);
    std::vector<double> fom(n, MissingValue), t(n);
    for (std::size_t j = 0; j < n; ++j) {
      BNBSpillInfo const& spill = spills[order[j]];
      S[j] = getBNBBeamState(spill);
      if (S[j].hasFOM()) fom[j] = computeFOM(S[j]);
      t[j] = spillTime(spill);
    }
    auto const valid = [&](std::size_t j){ return S[j].hasFOM(); };
    // the spills just before and after j, both within nbMaxGap
    auto const adjacent = [&](std::size_t j){
      return j > 0 && j + 1 < n && t[j] - t[j-1] <= cfg.nbMaxGap && t[j+1] - t[j] <= cfg.nbMaxGap;
    };
    std::vector<bool> changed(n, false);

    // ---- 1. width of an empty multiwire from the adjacent spills
    if (cfg.neighbourWidth) {
      auto const goodWidth = [&](std::size_t j)
        { return valid(j) && measuredWidth(S[j]) && !(S[j].status & WidthFromNeighbors); };
      for (std::size_t j = 0; j < n; ++j) {
        if (!valid(j) || !(S[j].status & MWEmpty) || !adjacent(j)) continue;
        BNBBeamState const &p = S[j-1], &q = S[j+1];
        if (!goodWidth(j-1) || !goodWidth(j+1)) continue;
        if (std::abs(p.sx - q.sx) > cfg.nbMaxDSig || std::abs(p.sy - q.sy) > cfg.nbMaxDSig) continue;
        S[j].sx = 0.5 * (p.sx + q.sx);
        S[j].sy = 0.5 * (p.sy + q.sy);
        S[j].widthSource = p.widthSource;
        S[j].status = (S[j].status & ~NoMWWidth) | WidthFromNeighbors;
        fom[j] = computeFOM(S[j]);
        changed[j] = true;
      }
    }

    // ---- 2. neighbour fill: a BPM reading is missing, the adjacent spills are stable
    if (cfg.neighbourFill) {
      std::vector<BNBBeamState> const S0 = S;    // neighbours as they were
      std::vector<double> const fom0 = fom;
      std::vector<bool> filled(n, false);
      for (std::size_t j = 0; j < n; ++j) {
        BNBBeamState& s = S[j];
        if (valid(j) || (s.status & NoTOR) || s.tor <= cfg.minTor || !adjacent(j)) continue;
        if (!(s.status & (NoHBPM | NoVBPM))) continue;
        BNBBeamState const &p = S0[j-1], &q = S0[j+1];
        if (!p.hasFOM() || !q.hasFOM() || filled[j-1]) continue;
        double const tm = 0.5 * (p.tor + q.tor);
        if (std::abs(s.tor - tm) / tm > cfg.nbMaxDTor) continue;
        if (std::max(std::abs(p.hpos - q.hpos), std::abs(p.vpos - q.vpos)) > cfg.nbMaxDPos) continue;
        if (std::max(std::abs(p.hang - q.hang), std::abs(p.vang - q.vang)) > cfg.nbMaxDAng) continue;
        if (measuredWidth(p) && measuredWidth(q)    // a nominal width on either side: no requirement
          && std::max(std::abs(p.sx - q.sx), std::abs(p.sy - q.sy)) > cfg.nbMaxDSig) continue;
        if (std::abs(fom0[j-1] - fom0[j+1]) > cfg.nbMaxDFOM) continue;
        if (cfg.nbMinFOM >= 0. && (fom0[j-1] <= cfg.nbMinFOM || fom0[j+1] <= cfg.nbMinFOM)) continue;
        s.hpos = 0.5 * (p.hpos + q.hpos); s.hang = 0.5 * (p.hang + q.hang);
        s.vpos = 0.5 * (p.vpos + q.vpos); s.vang = 0.5 * (p.vang + q.vang);
        if (!measuredWidth(s)) {                    // no width of its own
          if (measuredWidth(p) && measuredWidth(q)) {
            s.sx = 0.5 * (p.sx + q.sx); s.sy = 0.5 * (p.sy + q.sy);
            s.widthSource = p.widthSource;
            s.status = (s.status & ~NoMWWidth) | WidthFromNeighbors;
          }
        }
        // the "missing" bits stay set: they record why the spill was filled
        s.status |= NeighborFilled;
        BNBBeamState computed = s;
        computed.status &= ~(NoHBPM | NoVBPM);
        fom[j] = computeFOM(computed);
        filled[j] = changed[j] = true;
      }
    }

    // ---- 3. 875-station drop-outs: own target BPMs + the measured spills on both sides
    if (cfg.burstFill) {
      auto const isGood = [&](std::size_t j){
        BNBSpillInfo const& sp = spills[order[j]];
        return valid(j) && !(S[j].status & NeighborFilled)
          && isValid(sp.HP875) && isValid(sp.VP875) && isValid(sp.VP873)
          && isValid(sp.HPTG1) && isValid(sp.VPTG2);
      };
      std::vector<std::size_t> good;
      for (std::size_t j = 0; j < n; ++j) if (isGood(j)) good.push_back(j);
      for (std::size_t j = 0; j < n && !good.empty(); ++j) {
        BNBBeamState& s = S[j];
        BNBSpillInfo& sp = spills[order[j]];
        if (valid(j) || (s.status & (NoTOR | NeighborFilled)) || s.tor <= cfg.minTor) continue;
        if (isValid(sp.HP875) || isValid(sp.VP875) || !isValid(sp.HPTG1) || !isValid(sp.VPTG2)) continue;
        auto const it = std::lower_bound(good.begin(), good.end(), j);   // first good after j
        if (it == good.begin() || it == good.end()) continue;
        std::size_t const a = *(it - 1), b = *it;
        if (t[j] - t[a] > cfg.burstMaxDt || t[b] - t[j] > cfg.burstMaxDt) continue;
        if (fom[a] <= cfg.burstMinFOM || fom[b] <= cfg.burstMinFOM) continue;
        BNBBeamState const &A = S[a], &B = S[b];
        if (std::abs(A.hang - B.hang) >= cfg.burstMaxDAng || std::abs(A.vang - B.vang) >= cfg.burstMaxDAng) continue;
        BNBSpillInfo const &spA = spills[order[a]], &spB = spills[order[b]];
        if (std::abs(sp.HPTG1 - 0.5 * (spA.HPTG1 + spB.HPTG1)) >= cfg.burstMaxDTgt) continue;
        if (std::abs(sp.VPTG2 - 0.5 * (spA.VPTG2 + spB.VPTG2)) >= cfg.burstMaxDTgt) continue;
        BNBBeamState e = s;
        e.hpos = 0.5 * ((A.hpos - spA.HPTG1) + (B.hpos - spB.HPTG1)) + sp.HPTG1;
        e.hang = 0.5 * (A.hang + B.hang);
        e.vpos = 0.5 * ((A.vpos - spA.VPTG2) + (B.vpos - spB.VPTG2)) + sp.VPTG2;
        e.vang = 0.5 * (A.vang + B.vang);
        // own width unless the M876 database width says the chamber was empty
        bool const m876empty = (isValid(sp.M876HS) && sp.M876HS > 4.0) || (isValid(sp.M876VS) && sp.M876VS > 4.0);
        if (m876empty || !measuredWidth(e)) {
          e.sx = e.sy = MissingValue;
          e.widthSource = BNBBeamState::Nominal;
          e.status |= NoMWWidth;
        }
        e.status &= ~(NoHBPM | NoVBPM);
        double const f = computeFOM(e);
        if (!(f > cfg.burstAccept)) continue;
        unsigned int const missing = s.status & (NoHBPM | NoVBPM);   // kept: why the spill was filled
        s = e;
        s.status |= BurstFill | missing;
        fom[j] = f;
        changed[j] = true;
      }
    }

    std::vector<unsigned int> status(n);
    for (std::size_t j = 0; j < n; ++j) {
      status[order[j]] = S[j].status;
      if (!changed[j]) continue;
      BNBBeamState computed = S[j];
      computed.status &= ~(NoHBPM | NoVBPM);
      storeFOM(spills[order[j]], computed, fom[j]);
    }
    return status;
  }

} // namespace sbn
