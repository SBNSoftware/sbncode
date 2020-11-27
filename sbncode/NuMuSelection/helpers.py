import numpy as np
import awkward as ak

# Broadcast a numpy array (var) over a sequence (nbroadcast)
def broadcast(var, nbroadcast):
    return np.repeat(var, nbroadcast)

# Broadcast an awkward array (var) over a sequence (nbroadcast)
def broadcast_ak(var, nbroadcast):
    flat_var = var.flatten()
    flat_nbroadcast = np.repeat(nbroadcast, var.counts)

    flat_var_repeat = np.repeat(flat_var, flat_nbroadcast)
    
    # TODO: VECTORIZE THIS STEP!
    reindex = np.hstack([np.add.outer(np.arange(0, var.counts[i]) * nbroadcast[i], np.arange(0, nbroadcast[i])).flatten("F") for i in range(len(nbroadcast))])

    reindex += np.repeat(np.hstack([[0], np.cumsum(nbroadcast * var.counts)[:-1]]), nbroadcast * var.counts)
    flat_var_repeat_ordered = flat_var_repeat[reindex]
    return group(flat_var_repeat_ordered, np.repeat(var.counts, nbroadcast))

# Group a numpy array (var) into an awkward array by groups (ngroup)
def group(var, ngroup):
    return ak.JaggedArray.fromcounts(ngroup, var)

# normalize
def NeutrinoPOT(data):
    _, ind = np.unique(data["hdr.subrun"] + data["hdr.run"]*100, return_index=True)
    return np.sum(data["hdr.pot"][ind])

def NGenEvt(data):
    _, ind = np.unique(data["hdr.subrun"] + data["hdr.run"]*100, return_index=True)
    return np.sum(data["hdr.ngenevt"][ind])

def NEvt(data):
    return len(data["hdr.evt"])

    
POT_PER_SPILL = 5e12
def CosmicPOT(cosmic, nu):
    neutrino_per_spill = (NGenEvt(nu) * POT_PER_SPILL) / NeutrinoPOT(nu)
    _, ind = np.unique(nu["hdr.subrun"] + nu["hdr.run"]*100, return_index=True)
    n_cosmic_evt = NGenEvt(cosmic)
    return n_cosmic_evt * POT_PER_SPILL / (1. - neutrino_per_spill)
