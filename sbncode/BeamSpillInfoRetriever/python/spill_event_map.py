#!/usr/bin/env python3
"""
spill_event_map.py -- for every BNB spill, the run/subrun/event it belongs to.

  ./spill_event_map.py beamqual_run4.root [--fom beamqual_run4_fom.root] -o spill_event_map_run4.root

The retriever (ICARUSBNBRetriever / SBNDBNBRetriever) stores, for each triggered
event, every spill since the previous trigger, tagged with that event's number
(rec.hdr.bnbinfo.event). So each spill belongs to the event whose window
(previous trigger, this trigger] contains it; the last spill of the window is the
spill that triggered the event (rec.hdr.spillbnbinfo of that CAF record).

Output TTree "map", one row per spill (rec.hdr.bnbinfo):
  run, subrun, event          event that collected the spill
  spill_time_sec/nsec         spill time
  is_trigger_spill            1 if this is the spill the event triggered on
  spill_index, n_spills       position of the spill in its event's window (time order), window size
  trigger_time_sec/nsec       time of the event's trigger spill (0 if the CAF record has no spill info)
  event_in_caf                1 if the event exists as a CAF record with spill info
  duplicate                   1 for the 2nd+ copy of a spill that appears in more than one CAF file
  file_idx                    index into the dump's file_list
  FOM, FOM_best, fom_source, fom_status, TOR860   (FOM_* only with --fom)
"""
import argparse
import sys

import numpy as np
import pandas as pd
import uproot

UI32 = np.iinfo(np.uint32).max


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("beamqual", help="dump_bnb_quality.py output")
    ap.add_argument("--fom", help="matching bnb_fom.py output (adds FOM_best etc.)")
    ap.add_argument("-o", "--output", required=True, help=".root or .csv")
    a = ap.parse_args(argv)

    F = uproot.open(a.beamqual)
    b = pd.DataFrame(F["bnbinfo"].arrays(["run", "subrun", "event", "spill_time_sec", "spill_time_nsec",
                                          "file_idx", "TOR860", "FOM"], library="np"))
    s = pd.DataFrame(F["spill"].arrays(["run", "subrun", "evt", "event", "spill_time_sec", "spill_time_nsec"],
                                       library="np"))
    if a.fom:
        G = uproot.open(a.fom)["bnbinfo"]
        f = G.arrays(["spill_time_sec", "spill_time_nsec", "FOM_best", "fom_source", "fom_status"], library="np")
        if not (np.array_equal(f["spill_time_sec"], b.spill_time_sec.values)
                and np.array_equal(f["spill_time_nsec"], b.spill_time_nsec.values)):
            sys.exit("ERROR: --fom file is not row-aligned with the beamqual file")
        for k in ("FOM_best", "fom_source", "fom_status"):
            b[k] = f[k]

    b["t"] = b.spill_time_sec.astype(np.int64) * 10**9 + b.spill_time_nsec.astype(np.int64)
    key = ["run", "subrun", "event", "t"]
    b["duplicate"] = b.duplicated(key).astype(np.int8)

    ok = (s.event != UI32) & ~s.duplicated(["run", "subrun", "evt"])     # records with spill info
    trig = s.loc[ok, ["run", "subrun", "evt", "spill_time_sec", "spill_time_nsec"]].rename(
        columns={"evt": "event", "spill_time_sec": "trigger_time_sec", "spill_time_nsec": "trigger_time_nsec"})
    b["_row"] = np.arange(len(b))
    m = b.merge(trig, on=["run", "subrun", "event"], how="left").sort_values("_row")
    m["event_in_caf"] = m.trigger_time_sec.notna().astype(np.int8)
    for c in ("trigger_time_sec", "trigger_time_nsec"):
        m[c] = m[c].fillna(0).astype(np.uint64)
    tt = m.trigger_time_sec.astype(np.int64) * 10**9 + m.trigger_time_nsec.astype(np.int64)
    m["is_trigger_spill"] = ((m.event_in_caf == 1) & (m.t == tt)).astype(np.int8)

    u = m.duplicate == 0                                     # index/size over unique spills
    g = m[u].sort_values(["run", "subrun", "event", "t"]).groupby(["run", "subrun", "event"])
    idx = g.cumcount()
    m["spill_index"] = -1
    m.loc[idx.index, "spill_index"] = idx.values
    m["n_spills"] = m.groupby(["run", "subrun", "event"])["duplicate"].transform(lambda d: int((d == 0).sum()))
    m["spill_index"] = m.spill_index.astype(np.int32)
    m["n_spills"] = m.n_spills.astype(np.int32)

    cols = ["run", "subrun", "event", "spill_time_sec", "spill_time_nsec", "is_trigger_spill", "spill_index",
            "n_spills", "trigger_time_sec", "trigger_time_nsec", "event_in_caf", "duplicate", "file_idx",
            "TOR860", "FOM"] + (["FOM_best", "fom_source", "fom_status"] if a.fom else [])
    out = {c: m[c].to_numpy() for c in cols}
    if a.output.endswith(".root"):
        with uproot.recreate(a.output) as w:
            tree = w.mktree("map", {c: v.dtype for c, v in out.items()})
            for i in range(0, len(m), 200000):
                tree.extend({c: v[i:i + 200000] for c, v in out.items()})
    else:
        pd.DataFrame(out).to_csv(a.output, index=False)

    n = len(m)
    print(f"{a.output}: {n:,} spills -> {m[u].groupby(['run','subrun','event']).ngroups:,} events; "
          f"event in CAF {m.event_in_caf.mean():.4%}; trigger spills {int(m.is_trigger_spill.sum()):,}; "
          f"duplicate copies {int(m.duplicate.sum()):,}")


if __name__ == "__main__":
    main()
