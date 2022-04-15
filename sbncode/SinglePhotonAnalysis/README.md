# SinglePhotonModule 
Ported from `ubana` in suite `uboonecode v08_00_00_43 e17:prof` in Mar. 2022

The original code can be found [here](https://cdcvs.fnal.gov/redmine/projects/ubana/repository?utf8=%E2%9C%93&rev=feature%2Fmarkross_Nov2021_merge)

NEED INTRODUCTION HERE.

## Overview
`NCRadiativeResonant_module.cc` provides module `NCDeltaRadiative` for NCDeltaRadiative events production

`Libraries` contains essential headers for this module.

`SEAview` is an additional module runs inside this module.

`job` contains jobs for running the module.

`HelperFunctions` contains helper functions for simple calculations.

Header structure

```mermaid
	flowchart TD;
	A--first formost-->C;
	A[SinglePhoton_module.cc]--also contains-->B[analyze*.h that look like *.cxx];
	C[SinglePhoton_module.h]-->D[helper_functions.h];
	C-->E[Atlas.h];
	C-->F[SBSCAN.h];
	C-->G[SEAviewer.h];
```

