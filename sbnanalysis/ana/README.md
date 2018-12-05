# SelectionTool

Directory to hold the tool which fills the selection objects and functions

Works at the back end of the selection procedure and should not be included as 
part of the analysis of the events

Prevously used to fill ROOT TTrees in LArSoft which were passed to selection
objects outside of LArSoft

# SelectionObjects

Directory to hold the objects used in the selection procedure

These objects are filled in the SelectionTool and will be used by the
SelectionAnalysis

# SelectionAnalysis

Directory to hold the ROOT macros currently used for analysis

This will eventually change to be cleaner, but is just a hub for now

---------------------------------------------------------------------------------

These directories may eventually be spread across 'core' 'ana' and 'utils' in
sbnanalysis. Keeping it all in one place while it is merged
