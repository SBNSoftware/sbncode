
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
echo "setting up uboonecode version: v08_30_00"
setup uboonecode v08_30_00 -q e17:prof
unsetup larbatch

source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh    
echo "setting up sbndcode version: v08_30_00"
setup sbndcode v08_30_00 -q e17:prof
unsetup larbatch

source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
echo "setting up icaruscode version: v08_30_00"
setup icaruscode  v08_30_00 -q e17:prof 

if [ "sbnd"  == "sbnd" ]; then
    source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
elif [ "sbnd"  == "uboone" ]; then
    source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
elif [ "sbnd"  == "icarus" ]; then
    source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
else
    echo "Sorry. Are you not on a gpvm?"
    exit
fi

