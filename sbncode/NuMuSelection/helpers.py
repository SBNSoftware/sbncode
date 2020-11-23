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
