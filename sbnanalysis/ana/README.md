# Selection Functionality

## core

Added to the definition of Event and it's corresponding classes
- Metadata
- Neutrino
- FinalStateParticle
- Interaction
- RecoInteraction

Objects will be built to hold all interaction information given by LArSoft for the whole reconstruction chain

## util

Helper classes to hold
- PID functionality
- General analysis helper functions
- Other useful functionality for high-level analyses

## ana/SelectionAnalysis

Analysis macros which use selection functionality live here

---------------------------------------------------------------------------------------------------------------

This will eventually change to be cleaner, but is just a hub for now

---------------------------------------------------------------------------------

Location of the test reco file used to debug the merge:
/pnfs/sbnd/persistent/users/rsjones/reco-4960e937-cce4-4556-8591-09c05da633ae.root

Smaller testing reco file is here (with GENIE v3 implemented):
/pnfs/sbnd/persistent/users/rsjones/prodgenie_small_reco_test.root
