#!/usr/bin/env python3
"""
bunch_rotation.py -- Booster bunch-rotation settings for a time range, for update_flatcafs.py.

  B_BRRMPS[8]  scale factor of the bunch-rotation signal (0: bunch rotation off)   -> BRRMPS8
  B_BRRMPQ[1]  bunch-rotation frequency [kHz]                                     -> BRRMPQ1
  B_BRRMPL[9]  start time of the bunch-rotation curve [us]                        -> BRRMPL9

The settings are logged every 15 minutes. IFBeam has them from 2025-10-09 (one request per
day; longer windows fail). Earlier values come from the ACNet data logger (D44, node
"Backup", read-only, one request per week, >= 1.2 s apart); B_BRRMPL[9] is not logged there.
A spill gets the last value at or before its time, if that value is less than MAX_AGE old;
otherwise NaN.

  ./bunch_rotation.py 2025-04-01 2025-04-12      # print the settings changes in a range
"""
import csv
import datetime as dt
import io
import os
import sys
import time
import urllib.parse
import urllib.request

import numpy as np

DEVICES = [("BRRMPS8", "B_BRRMPS[8]"), ("BRRMPQ1", "B_BRRMPQ[1]"), ("BRRMPL9", "B_BRRMPL[9]")]
IFBEAM_URL = "https://dbweb1.fnal.gov:8443/ifbeam/data/data"
IFBEAM_START = dt.datetime(2025, 10, 9, tzinfo=dt.timezone.utc).timestamp()  # first IFBeam record
D44_URL = "https://www-bd.fnal.gov/D44/dataLogger"
D44_NODE, D44_EVENT = 96, "p,300000,true"    # Backup logger
MAX_AGE = 3600.0                             # [s] older settings are not used for a spill


def _get(url, timeout=120, tries=3):
    for k in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=timeout).read().decode()
        except Exception as e:               # noqa: BLE001 - network: retry, then give up
            print(f"[bunch_rotation] retry {k + 1}: {url[:90]}... ({e})", file=sys.stderr)
            time.sleep(5)
    return ""


def fetch_ifbeam(device, t0, t1, step=86400):
    ts, vs = [], []
    for a in np.arange(t0, t1, step):
        q = dict(v=device, e="p,900000", t0=int(a), t1=int(min(a + step, t1)), f="csv")
        body = _get(IFBEAM_URL + "?" + urllib.parse.urlencode(q))
        for row in csv.reader(io.StringIO(body)):
            if len(row) >= 5 and row[1] == device:
                try:
                    ts.append(float(row[2]) / 1000.); vs.append(float(row[-1]))
                except ValueError:
                    pass
    return ts, vs


def fetch_d44(device, t0, t1, step=7 * 86400):
    ts, vs = [], []
    for a in np.arange(t0, t1, step):
        b = min(a + step, t1)
        q = dict(t0=int(a * 1000), t1=int(b * 1000), how2plot="all", how_much="all", x="Time", y0=device,
                 y0_node=D44_NODE, y0_color="black", Y0_min="auto", Y0_max="auto", Y0_event=D44_EVENT, txt="true")
        body = _get(D44_URL + "?" + urllib.parse.urlencode(q, safe=":,[]"))
        for line in body.splitlines():
            p = line.split("\t")
            if len(p) == 2:
                try:
                    ts.append(a + float(p[0])); vs.append(float(p[1]))
                except ValueError:
                    pass
        time.sleep(1.2)
    return ts, vs


def fetch(t0, t1, path):
    """Settings for [t0 - MAX_AGE, t1] (UTC seconds), saved to and reused from `path` (.npz)."""
    if os.path.exists(path):
        return path
    t0 -= MAX_AGE
    out = {}
    for name, device in DEVICES:
        ts, vs = [], []
        if t0 < IFBEAM_START:
            a, b = fetch_d44(device, t0, min(t1, IFBEAM_START)); ts += a; vs += b
        if t1 > IFBEAM_START:
            a, b = fetch_ifbeam(device, max(t0, IFBEAM_START), t1); ts += a; vs += b
        o = np.argsort(ts, kind="stable")
        out[name + "_t"], out[name + "_v"] = np.asarray(ts)[o], np.asarray(vs)[o]
        print(f"[bunch_rotation] {device}: {len(ts)} values", file=sys.stderr)
    np.savez(path, **out)
    return path


def load(path):
    d = np.load(path)
    return {name: (d[name + "_t"], d[name + "_v"]) for name, _ in DEVICES}


def values_at(series, times):
    """Last value at or before each time (UTC s), NaN if none within MAX_AGE."""
    t, v = series
    times = np.asarray(times, float)
    if len(t) == 0:
        return np.full(len(times), np.nan, np.float32)
    i = np.searchsorted(t, times, side="right") - 1
    ok = (i >= 0) & (times - t[np.clip(i, 0, None)] <= MAX_AGE)
    return np.where(ok, v[np.clip(i, 0, None)], np.nan).astype(np.float32)


if __name__ == "__main__":
    a, b = (dt.datetime.fromisoformat(x).replace(tzinfo=dt.timezone.utc).timestamp() for x in sys.argv[1:3])
    path = fetch(a, b, "bunch_rotation_%d_%d.npz" % (a, b))
    for name, (t, v) in load(path).items():
        ch = ([0] if len(v) else []) + [k for k in range(1, len(v)) if v[k] != v[k - 1]]
        print(name, [(dt.datetime.fromtimestamp(t[k], dt.timezone.utc).strftime("%Y-%m-%d %H:%M"), float(v[k])) for k in ch][:20])
