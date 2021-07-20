# BNBFlux Uncertainty Evaluations

These codes are ported from MicroBooNE Flux Uncertainty Eventweight $v08_00_00_55$


## Structure
`FluxCalcPrep.*` implements the `FluxWeightCalc` class. `*WeightCalc()` are the weight calculators in different methods.

## FHiCL files
For each flux uncertainties, there are two parts of parameters: `Reweighting Environment` and `Calculator Settings`.

- Reweighting Environment
	- `type: Flux` request `WeightCalc` to use `FluxWeightCalc` class;
	- `random_seed` now still follows the one used in ubcode.
	- `mode: multisim` is used for flux uncertainty evaluations.
	- `number_of_multisims: 1000` sets the number of universes for the variation.
	- `parameter_list:[ "<flux uncertainty source>" ]` specifies a source of uncertainty.
	- `parameter_sigma: [n]` sets the variation with `n` sigma; it is only used in hadron production uncertainty evaluation
- Calculator Settings
	- `calc_type: "<methods>"` there are 5 methods to evaluate flux uncertainty weights: 
	Unisim, PrimaryHadronSWCentralSplineVariation, PrimaryHadronFeynmanScaling, 
	PrimaryHadronSanfordWang, and PrimaryHadronNormalization.
	- `PrimaryHadronGeantCode: [ptype]` is only used in hadron flux uncertainty;
	only particles with the specific parent type `ptype` would be evaluated.
	- `scale_factor_pos` is the scale factor to be applied when loading histograms from root files.
	- `scale_factor_neg` is another scale factor that is needed for the `Unisim` method.
	- `CentralValue_hist_file`, `PositiveSystematicVariation_hist_file`, and `NegativeSystematicVariation_hist_file` are MC files used in `Unisim`.
	- `ExternalData` and `ExternalFit` are data files used for hadron production uncertainty evaluations.

### External files needed
`Unisim` for non-hadron production uncertainty; it uses 
v
- CV simulation from BooBeamNT with parents redecayed 1000 times: *beamData/UnisimHists/may06_10kpot_ntrd1000_flux. root*; this simulation is reweighted to match up the may06 CV.
- Variation simulations are from MiniBooNE beam simulations by fine tuning one of the variables:
 - *beamData/UnisimHists/may06_horn175ka_rgen610.6_flux.root*	
 - *beamData/UnisimHists/may06_horn175ka_rgen610.6_flux.root*
 - *beamData/UnisimHists/may06_pioninexsec_up_rgen610.6_flux.root*	
 - *beamData/UnisimHists/may06_pioninexsec_down_rgen610.6_flux.root*
 - *beamData/UnisimHists/may06_pionqexsec_up_rgen610.6_flux.root*	
 - *beamData/UnisimHists/may06_pionqexsec_down_rgen610.6_flux.root*
 - *beamData/UnisimHists/may06_piontotxsec_up_rgen610.6_flux.root*	
 - *beamData/UnisimHists/may06_piontotxsec_down_rgen610.6_flux.root*
 - *beamData/UnisimHists/may06_nucleoninexsec_up_rgen610.6_flux.root*	
 - *beamData/UnisimHists/may06_nucleoninexsec_down_rgen610.6_flux.root*
 - *beamData/UnisimHists/may06_nucleonqexsec_up_rgen610.6_flux.root*	
 - *beamData/UnisimHists/ may06_nucleonqexsec_down_rgen610.6_flux.root*
 - *beamData/UnisimHists/may06_nucleontotxsec_up_rgen610.6_flux.root*	
 - *beamData/UnisimHists/may06_nucleontotxsec_down_rgen610.6_flux.root*
 - *beamData/UnisimHists/expskin_nrtd1000_flux.root*
 - *beamData/UnisimHists/expskin_nrtd1000_flux.root*

Hadron production uncertainty are from world data:
- MicroBooNE external data *beamData/ExternalData/BNBExternalData_uBooNE.root*
- HARP cross-section measurement *beamData/ExternalData/BNBExternalData_uBooNE_SplinesHARP.root*

#### Flux uncertainties breakdown
`Unisim` flux uncertainty reweights the following `parameter_list`:
- Horn current
- pion inelastic XS
- pion QE XS
- pion total XS
- Nucleon inelastic XS
- Nucleon QE XS
- Nucleon total XS
- Skin depth in horn

`PrimaryHadronSWCentralSplineVariation`
- piplus and piminus

`PrimaryHadronFeynmanScaling`
- kplus

PrimaryHadronSanfordWang  
- kzero

PrimaryHadronNormalization
- kzero
