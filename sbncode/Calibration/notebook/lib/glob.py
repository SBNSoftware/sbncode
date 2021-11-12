import glob
import numpy as np
import uproot
import pandas as pd
from tqdm.auto import tqdm
from multiprocessing import Pool
import multiprocessing
from . import names
import dill

class NTupleProc(object):
    def __init__(self, f=None, name="None"):
        self.name = name
        self.f = f

    def __call__(self, df):
        return self.f(df)

    def __bool__(self):
        return self.f is not None

    # make work with pickle for multiprocessing
    def __getstate__(self):
        return dill.dumps({"f": self.f, "name": self.name})

    def __setstate__(self, state):
        data = dill.loads(state)
        self.f = data["f"]
        self.name = data["name"]



def _makedf(dfs):
    if not isinstance(dfs, tuple):
        assert(isinstance(dfs, pd.DataFrame))
        dfs = [dfs]
    else:
        dfs = list(dfs)

    npad = [max([len(b.split(".")) for b in df.columns]) for df in dfs] 
    def pad(b, i):
        while len(b) < npad[i]:
            b.append("")
        return tuple(b)

    for i in range(len(dfs)):
        dfs[i].columns = pd.MultiIndex.from_tuples([pad(b.split("."), i) for b in dfs[i].columns])

    # set the index name if not present
    for i in range(len(dfs)):
        if len(dfs[i].index.names) == 1 and dfs[i].index.names[0] is None:
            dfs[i].index = dfs[i].index.set_names(["entry"])

    return dfs
        
def _loaddf(inp):
    fname, branches, index, applyf = inp
    with uproot.open(fname) as f:
        dfW = _makedf(f[names.folderW][names.tname].arrays(branches, library="pd"))
        dfE = _makedf(f[names.folderE][names.tname].arrays(branches, library="pd"))
        if applyf:
            dfW = applyf(*dfW)
            dfE = applyf(*dfE)
        else:
            dfW = dfW[0]
            dfE = dfE[0]
        # Set an index on the NTuple number to make sure we keep track of what is where
        dfW["__ntuple"] = index
        dfW.set_index("__ntuple", append=True, inplace=True)
        dfW = dfW.reorder_levels([dfW.index.nlevels-1] + list(range(0, dfW.index.nlevels-1)))

        dfE["__ntuple"] = index + 1
        dfE.set_index("__ntuple", append=True, inplace=True)
        dfE = dfE.reorder_levels([dfE.index.nlevels-1] + list(range(0, dfE.index.nlevels-1)))

        dfW = dfW.append(dfE)
    return dfW

def _process(inp):
    fname = inp[0]
    with uproot.open(fname) as f:
        return _do_process(f, *inp[1:])

def _do_process(rootf, branches, vars, whens, bins):
    hists = {}

    dfW = _makedf(rootf[names.folderW][names.tname].arrays(branches, library="pd"))
    dfE = _makedf(rootf[names.folderE][names.tname].arrays(branches, library="pd"))
    
    valW = [v(dfW) for v in vars]
    valE = [v(dfE) for v in vars]
    
    wW = [w(dfW) if w else None for w in whens]
    wE = [w(dfE) if w else None for w in whens]
    
    runs = dfW.meta.run.unique()
    
    hists["W"] = {}
    for r in runs:
        hists["W"][r] = {}
        for val, var in zip(valW, vars):
            hists["W"][r][var.name] = {}
            for w, when in zip(wW, whens):
                if w is None:
                    hist = np.histogram(val[dfW.meta.run == r], bins=bins)
                else:
                    hist = np.histogram(val[w & (dfW.meta.run == r)], bins=bins)

                hists["W"][r][var.name][when.name] = hist

    hists["E"] = {}
    for r in runs:
        hists["E"][r] = {}
        for val, var in zip(valE, vars):
            hists["E"][r][var.name] = {}
            for w, when in zip(wE, whens):
                if w is None:
                    hist = np.histogram(val[dfE.meta.run == r], bins=bins)
                else:
                    hist = np.histogram(val[w & (dfE.meta.run == r)], bins=bins)

                hists["E"][r][var.name][when.name] = hist

    return hists

class NTupleGlob(object):
    def __init__(self, g, branches):
        if isinstance(g, list):
            self.glob = g
        else:
            self.glob = glob.glob(g)
        self.branches = branches

    def dataframe(self, branches=None, maxfile=None, nproc=1, f=None):
        if nproc == "auto":
            nproc = multiprocessing.cpu_count()
        if branches is None:
            branches = self.branches

        thisglob = self.glob 
        if maxfile:
            thisglob = thisglob[:maxfile]

        ret = []
        with Pool(processes=nproc) as pool:
            thisglob = [(g, branches, i*2, f) for i,g in enumerate(thisglob)]
            for df in tqdm(pool.imap_unordered(_loaddf, thisglob), total=len(thisglob), unit="file", delay=5):
                ret.append(df)

        ret = pd.concat(ret, axis=0, ignore_index=False)

        # Fix the index So that we don't need __ntuple
        sub_index = ret.index.names[2:]
        ret = ret.reset_index()
        ret.entry = ret.groupby(["__ntuple", "entry"]).ngroup()
        ret.set_index(["entry"] + sub_index, inplace=True, verify_integrity=True)
        ret.sort_index(inplace=True)
        del ret["__ntuple"]

        return ret

    def histogram(self, var, bins, when=NTupleProc(), flatten_runs=False, flatten_cryo=False, maxfile=None, nproc=1):
        if nproc == "auto":
            nproc = multiprocessing.cpu_count()

        if not isinstance(var, list):
            var = [var]

        if not isinstance(when, list):
            when = [when]

        ret = {}

        thisglob = self.glob
        if maxfile:
            thisglob = thisglob[:maxfile]

        globdata = [(f, self.branches, var, when, bins) for f in thisglob]

        with Pool(processes=nproc) as pool:
            for hists in tqdm(pool.imap_unordered(_process, globdata), total=len(globdata), unit="file", delay=5):
                for cname in hists.keys():
                    for runname in hists[cname].keys():
                        for varname in hists[cname][runname].keys():
                            for whenname in hists[cname][runname][varname].keys():
                                hist = self._load_histogram(cname, runname, varname, whenname, ret)
                                if hist is None:
                                    ret[cname][runname][varname][whenname] = hists[cname][runname][varname][whenname]
                                else:
                                    ret[cname][runname][varname][whenname] = self._hadd(hist, hists[cname][runname][varname][whenname])

        # Do flattening
        flatret_cryo = {}
        if flatten_cryo:
            for runname in ret["E"].keys():
                flatret_cryo[runname] = {}
                for valname in ret["E"][runname].keys():
                    flatret_cryo[runname][valname] = {}
                    for whenname in ret["E"][runname][valname].keys():
                        flatret_cryo[runname][valname][whenname] = self._hadd(ret["E"][runname][valname][whenname], ret["W"][runname][valname][whenname])
            ret = flatret_cryo 

        flatret_run = {}
        if not flatten_cryo:
           flatret_run["E"] = {}
           flatret_run["W"] = {}

        if flatten_runs:
            histlist = [ret] if flatten_cryo else [ret["E"], ret["W"]]
            makeflatlist = [flatret_run] if flatten_cryo else [flatret_run["E"], flatret_run["W"]]

            for hists, makeflat in zip(histlist, makeflatlist):
                run0 = list(hists.keys())[0]

                for valname in hists[run0].keys():
                    makeflat[valname] = {}
                    for whenname in hists[run0][valname].keys():
                        makeflat[valname][whenname] = self._hadd(*[hists[runname][valname][whenname] for runname in hists.keys()])
            ret = flatret_run

        return ret

        if len(when) == 1 and not when[0]:
            if flatten_runs and flatten_cryo:
                for varname in ret.keys():
                    ret[varname] = ret[varname]["None"]
            elif flatten_runs:
                for cname in ret.keys():
                    for varname in ret[cname].keys():
                        ret[cname][varname] = ret[cname][varname]["None"]
            elif flatten_cryo:
                for runname in ret.keys():
                    for varname in ret[runname].keys():
                        ret[runname][varname] = ret[runname][varname]["None"]
            else:
                for cname in ret.keys():
                    for runname in ret[cname].keys():
                        for varname in ret[cname][runname].keys():
                            ret[cname][runname][varname] = ret[cname][runname][varname]["None"]

        return ret

    def _hadd(self, *hs):
        Ns = [N for N,_ in hs]
        return np.sum(Ns, axis=0), hs[0][1]
 
    def _load_histogram(self, cryo, run, valname, whenname, hist_dict):
        for level in [cryo, run, valname]:
            if level not in hist_dict:
                hist_dict[level] = {}
            hist_dict = hist_dict[level]

        if whenname not in hist_dict:
            return None
        else:
            return hist_dict[whenname]



