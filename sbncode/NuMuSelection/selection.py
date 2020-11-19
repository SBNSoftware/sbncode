import numpy as np
def InBeam(t): # us
    return (t > 0.) & (t < 1800) # allow some wiggle

def InFV(x, y, z): # cm
    xmin = -199.15 + 15
    ymin = -200. + 15
    zmin = 0.0 + 15
    xmax = 199.15 - 15
    ymax =  200. - 15
    zmax =  500. - 50
    
    return (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax) & (z > zmin) & (z < zmax)

def InAV(x, y, z): # cm
    xmin = -199.15 + 5
    ymin = -200. + 5
    zmin = 0.0 + 5
    xmax = 199.15 - 5
    ymax =  200. - 5
    zmax =  500. - 5
    
    return (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax) & (z > zmin) & (z < zmax)

# Primary track calculation
def get_primary_tracks(data):
    # these are the tracks that we consider as coming from the neutrino vertex
    from_slice = data["slc.reco.trk.atslc"] & data["slc.reco.trk.parent_is_primary"]
    # If there is an exiting vertex track, pick the longest such track as the primary
    ret = data["slc.reco.trk.len"][from_slice].argmax().max()
    
    # In cases where all the tracks are contained, first apply cuts on the particle ID
    maybe_muon = (data["slc.reco.trk.bestplane.chi2_proton"] > 60) & \
            (data["slc.reco.trk.bestplane.chi2_muon"] < 30) & (data["slc.reco.trk.len"] > 50)
    
    # Cases where all the tracks are contained
    all_contained = data["slc.reco.trk.contained"][from_slice].all()
    ret[all_contained] = data["slc.reco.trk.len"][maybe_muon].argmax().max()[all_contained]
    # map very small to -1
    ret[ret < 0] = -1
    # setup for index
    #ret = group(ret, np.ones(ret.shape, dtype="int"))
    return ret

# Using truth to get the primary track
def get_true_primary_track(data):
    ret = ((np.abs(data["slc.reco.trk.truth.p.pdg"]) == 13) & \
    (data["slc.reco.trk.truth.bestmatch.energy"] / data["slc.reco.trk.truth.p.planeVisE"] > 0.5)).argmax().max()
    is_numu_cc = (data["slc.truth.index"] >= 0) & data["slc.truth.iscc"] & (np.abs(data["slc.truth.pdg"]) == 14)
    # ignore non cc numu cases
    ret[np.invert(is_numu_cc)] = -1
    ret[ret<0] = -1
    # setup for index
    #ret = group(ret, np.ones(ret.shape, dtype="int"))
    return ret

# Define the Cuts!!!!!
def fid(data):
    return InFV(data["slc.vertex.x"], data["slc.vertex.y"], data["slc.vertex.z"])
def nu_score(data):
    return (data["slc.nu_score"] > 0.3) & np.invert(data["slc.is_clear_cosmic"])
def f_time(data):
    return InBeam(data["slc.fmatch.time"]*1000.) # ms -> us
def f_score(data):
    return data["slc.fmatch.score"] < 7
def good_reco(data):
    return (data["slc.ptrk.contained"]) | (data["slc.ptrk.len"] > 100.)
def crttrack(data):
    return np.isnan(data["slc.ptrk.crttrack.angle"]) | (data["slc.ptrk.crttrack.angle"] > 0.2)
def crthit(data):
    return np.isnan(data["slc.ptrk.crthit.distance"]) | (data["slc.ptrk.crthit.distance"] > 5) | (np.invert(InBeam(data["slc.ptrk.crthit.hit.time"])))
def length(data):
    return (data["slc.ptrk.len"] > 50.)
