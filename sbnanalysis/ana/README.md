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
/pnfs/sbnd/mc/reco/artroot/pre-production/MCP0.9/prodgenie_nu_singleinteraction_tpc_gsimple-configd-v1/sbndcode/v07_07_00_2_MCP0_9/run_number/00/00/00/06/subrun_number/00/00/00/00/reco-a095595c-af52-435a-9fb0-ea86e7f9f3d5.root
