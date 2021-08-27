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

# Get external information on electron lifetime
run_etaus = {}
with open("/icarus/app/users/gputnam/calib/plots2/etau_run_data.txt") as f:
    next(f) # Skip first (header) line
    for line in f:
        dat = line.split(" ")
        run_etaus[int(dat[0])] = [float(d) for d in dat[1:]]

plane2branches = [
    "h.time", "h.width", "h.tpc", "dqdx", "pitch", "rr", "dir.x",
]
plane2branches = ["hits2.%s" % s for s in plane2branches]

ray_branches = [
    "daughter_nsp",
    "daughter_pdg",
    "daughter_sp_toend_dist"
]

t0 = 0
def exp(t, *p):
    A,tau = p
    return A*np.exp(-(t - t0)/tau)

def isTPCE(df):
    return df.tpc <= 1

def reduce_df(df, raydf=None):
    # use the external input to build the t0
    ccross_t0_E = df.hit_max_time_p2_tpcE - tcathode_EE
    ccross_t0_E[df.cryostat==1] = df.hit_max_time_p2_tpcE - tcathode_WE

    ccross_t0_W = df.hit_max_time_p2_tpcW - tcathode_EW
    ccross_t0_W[df.cryostat==1] = df.hit_max_time_p2_tpcW - tcathode_WW

    # Select stopping tracks
    select_track = df.selected == 0

    df["ccross_t0"] = ((ccross_t0_E + ccross_t0_W) / 2.) * tick_period
    df["tpcE"] = isTPCE(df.hits2.h)
    
    # What to save
    outdf = df.loc[(df.hits2.dqdx > 0) & (df.hits2.rr < 200.) & select_track, 
                  [ ("hits2", "h", "time"),
                    ("tpcE", "", ""),
                    ("hits2", "dqdx", ""),
                    ("hits2", "pitch", ""),
                    ("hits2", "rr", ""),
                    ("hits2", "dir", "x"),
                    ("hits2", "h", "width"),
                    ("ccross_t0", "", ""),
                    ("meta", "run", ""),
                    ("cryostat", "", ""),
                    ("dir", "y", ""),
                    ("hit_min_time_p2_tpcE", "", ""),
                    ("hit_max_time_p2_tpcE", "", ""),
                    ("hit_min_time_p2_tpcW", "", ""),
                    ("hit_max_time_p2_tpcW", "", ""),
                   ]
                  ].copy()

    # Simplify column names
    outdf.columns = ["time", "tpcE", "dqdx_nocorr", "pitch", "rr", "dir_x", "width", "ccross_t0", "run", "cryostat", "tdir_y", "mint_tpcE", "maxt_tpcE", "mint_tpcW", "maxt_tpcW"]
    
    # Correct for electron lifetime
    outdf["thit"] = (outdf.time * tick_period - outdf.ccross_t0 - tanode*tick_period) / 1000.
    if len(outdf):
        thisrun = outdf.run.iloc[0]
        # Correct in each TPC
        outdf["dqdx_corr"] = outdf.dqdx_nocorr * exp(outdf.thit, 1., -run_etaus[thisrun][0]*1e3)
        outdf.loc[~outdf.tpcE & (outdf.cryostat==0), "dqdx_corr"] = (outdf.dqdx_nocorr * exp(outdf.thit, 1., -run_etaus[thisrun][1]*1e3))[~outdf.tpcE & (outdf.cryostat==0)]
        outdf.loc[outdf.tpcE &  (outdf.cryostat==1), "dqdx_corr"] = (outdf.dqdx_nocorr * exp(outdf.thit, 1., -run_etaus[thisrun][2]*1e3))[outdf.tpcE &  (outdf.cryostat==1)]
        outdf.loc[~outdf.tpcE & (outdf.cryostat==1), "dqdx_corr"] = (outdf.dqdx_nocorr * exp(outdf.thit, 1., -run_etaus[thisrun][3]*1e3))[~outdf.tpcE & (outdf.cryostat==1)]

    # Save information on PFP daughters
    if raydf is not None:
        closest_tdaughter = raydf.groupby(level=0).daughter_sp_toend_dist.min()
        outdf = outdf.join(closest_tdaughter.rename("closest_tdaughter"), on="entry", how="left")

    return outdf

def main(output, inputs):
    ntuples = NTupleGlob(inputs, branches.trkbranches + plane2branches + ray_branches)
    df = ntuples.dataframe(nproc="auto", f=reduce_df)
    df.to_hdf(output, key="df", mode="w")

if __name__ == "__main__":
    printhelp = len(sys.argv) < 3 or sys.argv[1] == "-h"
    if printhelp:
        print("Usage: python make_calib_df.py [output.df] [inputs.root,]")
    else:
        main(sys.argv[1], sys.argv[2:])
