import sys
from ROOT import gSystem
gSystem.Load("libOpFlashAna_Algorithms")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing Algorithms..."

sys.exit(0)

