#!/usr/bin/env python3
"""Read-only search of ACNet data loggers (via the D44 web endpoint) for
BNB BPM readings missing from IFBeam in ICARUS Run 2 (2023).

Every reply is cached under raw/d44/, and uncached requests are spaced >= 1.1 s apart.
Usage: python3 d44_search.py  -> prints a TSV summary, also written to results.tsv
"""
import hashlib, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

BASE = "https://www-bd.fnal.gov/D44"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "raw", "d44")
MAX_NEW = int(os.environ.get("MAX_NEW", "300"))   # ask before exceeding

WINDOWS = {
    1: dict(t0=1686355200000, t1=1686355800000, run=9941,
            missing=["E:HP875", "E:VP875", "E:VP873"], control="E:HPTG1"),
    2: dict(t0=1678924800000, t1=1678925400000, run=9649,
            missing=["E:HPTG1", "E:HPTG2", "E:VPTG2"], control="E:HP875"),
}
DEVICES = ["E:HP875", "E:VP875", "E:VP873", "E:HPTG1", "E:HPTG2", "E:VPTG2"]

_new = 0
_last = 0.0


def fetch(url, tag):
    """GET url with on-disk cache and rate limiting."""
    global _new, _last
    os.makedirs(CACHE, exist_ok=True)
    h = hashlib.sha1(url.encode()).hexdigest()[:10]
    path = os.path.join(CACHE, f"{re.sub(r'[^A-Za-z0-9_.-]', '_', tag)}__{h}.txt")
    if os.path.exists(path):
        return open(path).read()
    if _new >= MAX_NEW:
        sys.exit(f"Reached MAX_NEW={MAX_NEW} new queries; stopping.")
    wait = 1.1 - (time.time() - _last)
    if wait > 0:
        time.sleep(wait)
    _last = time.time()
    _new += 1
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            body = r.read().decode("utf-8", "replace")
    except Exception as e:  # cache errors too, so they are visible in raw/
        body = f"ERROR: {e!r}"
    with open(path, "w") as f:
        f.write(f"# URL: {url}\n{body}")
    return f"# URL: {url}\n{body}"


def node_list(device):
    """(source_id, logger_name, event) tuples the D44 page offers for a device."""
    txt = fetch(f"{BASE}/DLPNodes?devInfo=true&device={device}", f"nodes_{device}")
    root = ET.fromstring(txt.split("\n", 1)[1])
    return [(n.get("source_id"), n.get("logger_source"), n.get("event"))
            for n in root.iter("dlnode")]


def data_url(device, node, event, t0, t1, how2plot="all", how_much="all"):
    q = dict(t0=t0, t1=t1, how2plot=how2plot, how_much=how_much, x="Time",
             y0=device, y0_node=node, y0_color="black",
             Y0_min="auto", Y0_max="auto", Y0_event=event, txt="true")
    return f"{BASE}/dataLogger?" + urllib.parse.urlencode(q, safe=":,")


def parse(body):
    """Return (status, n_readings, example) from a dataLogger txt reply."""
    b = body.split("\n", 1)[1] if body.startswith("# URL:") else body
    if b.startswith("ERROR"):
        return "error", 0, b[:80]
    if "no data" in b.lower():
        return "No data", 0, ""
    rows = [l for l in b.splitlines()
            if re.search(r"\d", l) and re.search(r"-?\d+\.\d+\s*$", l.strip())]
    if not rows:
        return "error", 0, b.strip()[:80].replace("\n", " | ")
    return "data", len(rows), rows[0].strip()


def query(device, node, event, w):
    W = WINDOWS[w]
    body = fetch(data_url(device, node, event, W["t0"], W["t1"]),
                 f"w{w}_{device}_n{node}_{event}")
    return parse(body)


def main():
    out = []
    # (node, event) combos from the per-device node lists
    combos = {}
    for d in DEVICES:
        for sid, name, ev in node_list(d):
            combos.setdefault((sid, name, ev), set()).add(d)
    for (sid, name, ev), devs in sorted(combos.items()):
        for w in (1, 2):
            ctrl = WINDOWS[w]["control"]
            if ctrl not in devs:
                continue
            st = query(ctrl, sid, ev, w)
            out.append((sid, name, ev, ctrl, w, "control") + st)
            if st[0] != "data":
                continue
            for d in WINDOWS[w]["missing"]:
                if d in devs:
                    out.append((sid, name, ev, d, w, "missing") + query(d, sid, ev, w))
    hdr = "node\tlogger\tevent\tdevice\twindow\trole\tstatus\tn\texample"
    lines = [hdr] + ["\t".join(map(str, r)) for r in out]
    open(os.path.join(HERE, "results.tsv"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"# new network queries this run: {_new}", file=sys.stderr)


if __name__ == "__main__":
    main()
