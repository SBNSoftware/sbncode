# Some details of the selection merge

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
- RecoUtils from LArSoft
- Analysis-based functions

## selection

This holds the bulk of the code for the topological selection as it exists currently (Sept. 2019)
- See selection/README.md for details

## ana/SelectionAnalysis

Analysis macros which use selection functionality live here

---------------------------------------------------------------------------------------------------------------

This will eventually change to be cleaner, but is just a hub for now

---------------------------------------------------------------------------------

Location of the test reco file used to debug the merge:
/pnfs/sbnd/persistent/users/rsjones/reco-4960e937-cce4-4556-8591-09c05da633ae.root

Smaller testing reco file is here (with GENIE v3 implemented):
/pnfs/sbnd/persistent/users/rsjones/prodgenie_small_reco_test.root
