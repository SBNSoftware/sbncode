#!/usr/bin/env python

import cafana

fname = '/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/test/icarus/gen-prodcorsika_genie_nooverburden__nuetest.caf.root'

loader = cafana.SpectrumLoader(fname)

kTruthEnergy = cafana.CSliceVar('if(sr.truth.index > 0) return sr.truth.E; else return -1.;')

binsEnergy = cafana.Binning.Simple(50, 0, 5)

axEnergy = cafana.HistAxis('True energy (GeV)', binsEnergy, kTruthEnergy)
sEnergy = cafana.Spectrum(loader, axEnergy, cafana.kNoCut)

loader.Go()

sEnergy.ToTH1(sEnergy.POT()).Draw('hist')
