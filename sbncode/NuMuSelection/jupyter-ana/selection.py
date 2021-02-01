import numpy as np
from helpers import *

def InBeam(t): # us
    return (t > 0.) & (t < 1.800) # allow some wiggle

def InBeamVeto(t): # us
    return (t > 0.) & (t < 2.200) # allow some wiggle

def InBeamTrue(t): # us
    return (t > 0.) & (t < 1.596) # No wiggle!

def InFV(x, y, z): # cm
    xmin = -199.15 + 10
    ymin = -200. + 10
    zmin = 0.0 + 10
    xmax = 199.15 - 10
    ymax =  200. - 10
    zmax =  500. - 50
    
    return (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax) & (z > zmin) & (z < zmax)

# Primary track calculation
def get_primary_tracks(data):
    # these are the tracks that we consider as coming from the neutrino vertex
    from_slice = data["slc.reco.trk.atslc"] & data["slc.reco.trk.parent_is_primary"]
    
    maybe_muon_exiting = np.invert(data["slc.reco.trk.contained"]) & (data["slc.reco.trk.len"] > 100)
    maybe_muon_contained = data["slc.reco.trk.contained"] & (data["slc.reco.trk.bestplane.chi2_proton"] > 60) & (data["slc.reco.trk.bestplane.chi2_muon"] < 30) & (data["slc.reco.trk.len"] > 50)
    # In cases where all the tracks are contained, first apply cuts on the particle ID
    maybe_muon = from_slice & (maybe_muon_contained | maybe_muon_exiting)
    has_maybe_muon = maybe_muon.any()
    
    # Cases where all the tracks are contained
    ret = (data["slc.reco.trk.len"]*maybe_muon).argmax().max()
    # map very small to -1
    ret[ret < 0] = -1
    ret[np.invert(has_maybe_muon)] = -1
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
def crtveto(data):
    return broadcast(np.invert(InBeamVeto(data["crt_hits.time"]).any()), data["nslc"])
def crttrackveto_perevt(data):
    return np.invert(InBeamVeto(data["crt_tracks.time"]).any())
def crttrackveto(data):
    return broadcast(crttrackveto_perevt(data), data["nslc"])

def crttrackveto_nobottom_perevt(data):
    return np.invert((InBeamVeto(data["crt_tracks.time"]) & (data["crt_tracks.hita.position.y"] > -357.) & (data["crt_tracks.hitb.position.y"] > -357.)).any())
def crttrackveto_nobottom(data):
    return broadcast(crttrackveto_nobottom_perevt(data), data["nslc"])

def fid(data):
    return InFV(data["slc.vertex.x"], data["slc.vertex.y"], data["slc.vertex.z"])
def nu_score(data):
    return (data["slc.nu_score"] > 0.4) & np.invert(data["slc.is_clear_cosmic"])
def f_time(data):
    return InBeam(data["slc.fmatch.time"]) 
def f_score(data):
    return data["slc.fmatch.score"] < 7
def ptrk(data):
    return data["slc.has_ptrk"] & np.invert(np.isnan(data["slc.ptrk.recop"])) & (data["slc.ptrk.recop"] < 7.5) & (data["slc.ptrk.recop"] > 0.)
def crttrack(data):
    return np.isnan(data["slc.ptrk.crttrack.angle"]) | (data["slc.ptrk.crttrack.angle"] > 0.05)
def crthit(data):
    return np.isnan(data["slc.ptrk.crthit.distance"]) | (data["slc.ptrk.crthit.distance"] > 5) | InBeam(data["slc.ptrk.crthit.hit.time"])
