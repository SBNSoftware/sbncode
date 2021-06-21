 # PID Training Instructions
 To train the MVA run the following steps:

 - Run the PID module(s) with "MakeTree" set to true
 - Run the provided ROOT macros to produce the output weights files
 - Set the new weights files to as the input to the module, these will need to be in the FW_SEARCH_PATH (probably on experiment_data)
