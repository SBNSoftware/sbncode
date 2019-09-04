# Contents of this directory

Some of these may get moved around the sbnanalysis directory when I have finished, if somewhere else seems more appropriate (mostly separating out utilities from core selection features)

## Event

- Contains all objects associated with an 'Event' under the selection tool definition

## EventSelectionHelper

- Functions which contribute to the selection itself
- Includes PID for tracks & (basic, temporary procedure) for showers

## GeneralAnalysisHelper

- Functions which begin to apply the PID to analyses
- This will probably end up moving to util

## LoadEvents

- This will aid the filling of Particle and Event objects for use in the analysis functions

## Particle

- Contains all objects associated with an 'Particle' under the selection tool definition

## Plane 

- Functions to help determine whether a track or vertex is contained in the TPC or fiducial volume

## SelectionToolBase

- This will be the core mediator between the ProcessorBase output and the selection
