import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from seltool import ertool
from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("ertool_hist.root")

# Create ERTool algorithm (empty base class for this example)
my_algo = ertool.AlgoBase()

# Create ERTool analysis (empty base class for this example)
my_ana = ertool.AnaBase()

# Create larlite interfce analysis unit for ERTool
my_anaunit = fmwk.ExampleERSelection()

# Set Producers
# First Argument: True = MC, False = Reco
# Second Argument: producer module label
my_anaunit.SetShowerProducer(True,"mcreco");
my_anaunit.SetTrackProducer(True,"mcreco");
my_anaunit.SetVtxProducer(True,"generator");

# Implement manager
my_anaunit._mgr.AddAlgo(my_algo)
my_anaunit._mgr.AddAna(my_ana)
my_ana._mode =True # True = Select. False = Fill mode
my_proc.add_process(my_anaunit)

# run!
my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

#my_algo.StoreParams()
sys.exit(0)

