# Ported from MicroBooNE Flux Uncertainty Eventweight $v08_00_00_55$

One header with five calculators.


## Structure
`Configure()` calls different `Configure*()` for specific calculators

`GetWeight()` calls different `*Calc()` for different types of uncertainties

### Functions
#### Virtual functions are defined in `/sbncode/App/` directiory
- `void Configure()` to prepare the calculator.
- `GetWeight()` to evaluate weights for 1000 universes (# of universes can be adjusted).
	- Param `inu` is the ith set of universes.

#### `*WeightCalc.cxx` 
Each of them contains only one `FluxWeightCalc::*WeightCalc()` function,
which carrys out the weight evaluation.


#### Non-hadron productions

Unisim flux uncertainty reweights the following `parameter_list`:

- Horn current
- pion inelastic XS
- pion QE XS
- pion total XS
- Nucleon inelastic XS
- Nucleon QE XS
- Nucleon total XS
- Skin depth in horn


