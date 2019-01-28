echo "Your Experiment is sbnd"

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

if [ "sbnd"  == "sbnd" ]; then
    source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
elif [ "sbnd"  == "uboone" ]; then
    source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
elif [ "sbnd"  == "icarus" ]; then
    source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
else
    echo "Sorry. Are you not on a gpvm?"
fi

