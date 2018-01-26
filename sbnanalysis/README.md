SBN Analysis
============
This package provides a lightweight event analysis and fitting framework for
measurements within the the Short-Baseline Neutrino program at Fermilab.

It is built on [gallery](http://art.fnal.gov/gallery/), a library designed
for reading [art](http://art.fnal.gov) ROOT files.

Analysis code is stored in subdirectories of `ana`, with one directory per
analysis. Common core code for event I/O is in `io`, and shared utilities
in `util`.

Usage
-----
A primary goal of this package is to centralize development of event selection
code and fitting frameworks, and establish a standard analysis tree format
that is compatible with all fitters. As such, the fundamental object in
`sbnanalysis` is an event processor class which writes a standard ROOT tree,
plus additional user-specified fields.

Example
-------
An example analysis package is provided in `ana/ExampleAnalysis`.

Authors
-------
* A. Mastbaum, UChicago
* G. Putnam, UChicago
* J. Zennamo, UChicago
* D. Schmitz, UChicago
* *Your name here!*

