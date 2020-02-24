import time 
import sys
import os.path

def run(filelists):
    for flist in filelists:
        with open(flist) as f:
            for line in f:
                checkfile(line.rstrip("\n"))

# make the file seem busy
def checkfile(fname):
    if os.path.isfile(fname):
        print fname

def hello(lst):
    print "Running! On files:"
    for name in lst: print name

if __name__ == "__main__":
    hello(sys.argv[1:])
    run(sys.argv[1:])
