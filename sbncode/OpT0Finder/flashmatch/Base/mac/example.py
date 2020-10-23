import sys
from ROOT import gSystem
gSystem.Load("libOpFlashAna_Base")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing Base..."

sys.exit(0)

