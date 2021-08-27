import sys
from lib.glob import NTupleGlob
from lib import branches

def main(output, inputs):
    ntuples = NTupleGlob(inputs, branches.trkbranches)
    df = ntuples.dataframe(nproc="auto")
    df.to_hdf(output, key="df", mode="w")

if __name__ == "__main__":
    printhelp = len(sys.argv) < 3 or sys.argv[1] == "-h"
    if printhelp:
        print("Usage: python make_driftV_df.py [output.df] [inputs.root,]")
    else:
        main(sys.argv[1], sys.argv[2:])
