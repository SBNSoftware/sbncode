import sys
import datetime as dt
from lib.glob import NTupleGlob
from lib import branches
import numpy as np

# load constants
from lib.constants import *

# EXTERNAL INPUT: The drift window in each TPC
tcathode_EE = 3198.5279397664003
tcathode_EW = 3207.147982327826
tcathode_WE = 3200.883742841676
tcathode_WW = 3199.9763136348492

# Get external information on run timing
run_times = {}
with open("/icarus/app/users/gputnam/calib/rundata") as f:
    for line in f:
        dat = line.split(" ")
        run_times[int(dat[0])] = dt.datetime.strptime(dat[1].rstrip("\n"), "%Y-%m-%dT%H:%M:%S").date()

plane2branches = [
    "h.p.x", "h.p.y", "h.p.z", "h.time", "h.tpc", "dqdx", "dir.x", "dir.y", "dir.z",
]
plane2branches = ["hits2.%s" % s for s in plane2branches]

def isTPCE(df):
    return df.tpc <= 1

def reduce_df(df):
    # use the external input to build the t0
    ccross_t0_E = df.hit_max_time_p2_tpcE - tcathode_EE
    ccross_t0_E[df.cryostat==1] = df.hit_max_time_p2_tpcE - tcathode_WE

    ccross_t0_W = df.hit_max_time_p2_tpcW - tcathode_EW
    ccross_t0_W[df.cryostat==1] = df.hit_max_time_p2_tpcW - tcathode_WW

    # Select anode + cathode crossing tracks
    select_track = df.selected == 1

    df["ccross_t0"] = ((ccross_t0_E + ccross_t0_W) / 2.) * tick_period
    
    df = df[(df.hits2.dqdx > 0) & select_track].copy()
    df["chunk"] = df.index.get_level_values(1) // 5
    df["tpcE"] = isTPCE(df.hits2.h)
    outdf = df.groupby(["entry", "chunk"])[('hits2', 'dqdx', '', ''),].median()
    outdf = outdf.join(df.groupby(["entry", "chunk"])[[('hits2', 'h', 'p', 'x'),
                                                      ('hits2', 'h', 'p', 'y'),
                                                      ('hits2', 'h', 'p', 'z'),
                                                      ('hits2', 'h', 'time', ''),
                                                      ('hits2', 'dir', 'x', ''),
                                                      ('hits2', 'dir', 'y', ''),
                                                      ('hits2', 'dir', 'z', ''),]].mean())
    
    outdf.columns = ["dqdx", "x", "y", "z", "time", "dirx", "diry", "dirz"]
    
    # fix the direction normalization
    norm = np.sqrt(outdf.dirx**2 + outdf.diry**2 + outdf.dirz**2)
    outdf.dirx = outdf.dirx / norm
    outdf.diry = outdf.diry / norm
    outdf.dirz = outdf.dirz / norm

    # TPC/Cryo info
    outdf["tpcE"] = df.groupby(["entry", "chunk"]).tpcE.all()
    outdf["tpcW"] = (~outdf.tpcE) & (df.groupby(["entry", "chunk"]).tpcE.nunique() == 1)

    # Save T0
    outdf["ccross_t0"] = df.groupby(["entry", "chunk"]).ccross_t0.first()
    # also save the cryostat number
    outdf["cryostat"] = df.groupby(["entry", "chunk"]).cryostat.first()
    # And run number
    outdf["run"] = df.groupby(["entry", "chunk"])[[('meta', 'run', '', ''),]].first()

    # Only save chunks that are all in one TPC
    outdf = outdf[outdf.tpcE | outdf.tpcW].drop(columns=["tpcW"])
   
    return outdf


def main(output, inputs):
    ntuples = NTupleGlob(inputs, branches.trkbranches + plane2branches)
    df = ntuples.dataframe(nproc="auto", f=reduce_df)
    df.to_hdf(output, key="df", mode="w")

if __name__ == "__main__":
    printhelp = len(sys.argv) < 3 or sys.argv[1] == "-h"
    if printhelp:
        print("Usage: python make_etau_df.py [output.df] [inputs.root,]")
    else:
        main(sys.argv[1], sys.argv[2:])
