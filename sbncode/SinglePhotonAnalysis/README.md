# SinglePhotonModule 
Ported from `ubana` in suite `uboonecode v08_00_00_43 e17:prof` in Mar. 2022

The original code can be found at [Fermilab Readmine](https://cdcvs.fnal.gov/redmine/projects/ubana/repository?utf8=%E2%9C%93&rev=feature%2Fmarkross_Nov2021_merge)

NEED INTRODUCTION HERE.

* [Quick Start](#quick-start)
* [Update Log](#update-log)
* [Overview](#overview)
* [Glossary](#glossary)

---
## Quick Start
```
mrb i -j4
lar -c run_singlephoton_sbnd.fcl -s <input_reco2_sbnd_ROOT_FILE>
```

## Update Log
The module does not work out of the box, so there are some modifications to make it fit into the SBND.

- June 2022
	- New electric gains for SBND MC at `singlephoton_sbnd.fcl`, see SBN-doc-19505-v1;
	- Update the neutrino slice definition in the `Single Photon` module.
	- Update on geometry variables such as, `reco_*_dist_to_CathodePlane`
	- Disable Kalman dEdx variables
	- Disable Second Shower Search 3D (`sss3d`) variables
	- Adapate SBND FHiCLs for filters: `NCRadiativeResonant` and `NCDeltaRadiative`

These updates are to accommodate the change of Pandora features listed below:

### Pandora features updates
|Items|Pandora @ MicroBooNE|Pandora @ SBND|
|---|---|---|
|Neutrino Slice|One nu slice in each event | Multiple nu slices in each event|
|Kalman Fitter| In-use | Not in-use|
|3dShowers Objects| Available | Unavailable|
|MVA Method (for track/neutrino scores)|Support Vector Machines|Boosted Decision Trees|


---
## Overview

Three sub-modules are included in the Single Photon Analysis Module:
- `SinglePhoton_module.cc` impements the `SinglePhoton` module to read reco2 files and produce n-tuples with varaibles targeting photon reconstructions.
- `NCRadiativeResonant_module.cc`  implements the `NCRadiativeResonant` filter to select events with photons coming out from the nucleus.
- `NCDeltaRadiative_module.cc` implements the `NCDeltaRadiative` filter for NCDeltaRadiative events

### Code Structure
- `Libraries/` contains essential headers for the `SinglePhoton` module.

- `SEAview/` is an additional module runs inside the `SinglePhoton` module.

- `jobs/` contains FHiCL files for running these modules

- `HelperFunctions/` contains some useful functions to simplify the code

### Headers structure

```mermaid
flowchart TD;
A--first formost-->C;
A[SinglePhoton_module.cc]--contains-->B[analyze*.h that look like *.cxx];
C[SinglePhoton_module.h]-->D[helper_*.h];
C-->F[DBSCAN.h];
C-->G[SEAviewer.h];
```

### Flows of the `SinglePhoton` module

```mermaid
flowchart TB
A[ClearVertex & prepare PFParticles objects]-->AA;
subgraph AA[1. Gather Objects]
direction TB
B[AnalyzeSlices]-->C;
C[AnalyzeFlashes]-->D;
D[AnalyzeTracks]-->E;
E[AnalyzeShowers]-->F;
F[AnalyzeGeant4]-->G;
G[BuildParticleHitMaps];
end  

subgraph BB[2. Reco MC Matching]
direction TB
H[Match Reco. Showers to MCParticles]-->I;
I[Match Reco. Tracks to MCParticles]-->J;
J[CollectMCTruth Information];
end

subgraph CC[3. Gather Other Info. ]
direction TB
K[AnalyzeEventWeight]-->L;
L[AnalyzeRecoSlices];
end

subgraph DD[4. 2nd Shower earch]
direction TB
M[Isolation Study for 2nd Shower Veto]-->N;
N[Match 2nd Shower]-->O[Search 2nd Shower];
end

subgraph THIS[flow chart]
direction LR
AA-->BB
CC-->DD;
end
DD-->End[Output ROOT];
```


### Pandora Dependency
Objects are obtained from Pandora reconstruction, and they are connected  via the following:
```mermaid
graph TD
PF[PFParticle];T[Track];Sh[Shower];
H[Hit]; C[Cluster]; S[Slice];MCT[MCTruth];MCP[MCParticle];ID[ID#'s];M[Metadata];SP[SpacePoint];
PF-->ID;
PF-->M;
PF-->SP;
PF-->C;
C-->H;
PF-->Sh;
PF-->T;
PF-->S;
S-->H;
Sh-->MCP;
T-->MCP;
MCP-->MCT;
```

These objects are connected via labels shown as the following:
|Alias|Objects|Label|
|---|---|---|
|PandoraLabel|`std::vector<recob::PFParticle>`<br>`std::vector<recob::Cluster>`<br>`std::vector<recob::Slice>`<br>`std::vector<recob::Vertex>`|pandora|
|TrackLabel|`std::vector<recob::Track>`|pandoraTrack|
|ShowerLabel|`std::vector<recob::Shower>`|pandoraShower|
|ParticleIDLabel|`anab::ParticleID`(`art::FindOneP`)|pandoraSCEPid|
|CaloLabel|`anab::Calorimetry`(`art::FindManyP`)|pandoraSCECalo|
|FlashLabel|`std::vector<recob::OpFlash>`|opflashtpc0|
|POTLabel|`sumdata::POTSummary`|generator|

Below alias auto-configured by default
|Alias|Objects|Label|
|---|---|---|
|HitFinderModule|`std::vector<recob::Hit>`|gaushit|
|BadChannelLabel|`std::vector<int>`|badmasks(REMOVED)|
|ShowerTrackFitter|`art::Assns<recob::PFParticle,recob::Track,void>`|pandoraTrack|
|ShowerTrackFitterCalo|`anab::Calorimetry`(`art::FindManyP`)|pandoraCalo|
|GeneratorLabel|`std::vector<simb::GTruth>`<br>`std::vector<simb::MCTruth>`|generator|
|GeantModule|`simb::MCParticle>`(`art::FindManyP`)|largeant|
|HitMCParticleAssnLabel|`simb::MCParticle,anab::BackTrackerHitMatchingData`<br>(`art::FindManyP`)|gaushitTruthMatch|
|Shower3DLabel|`recob::Shower`(`art::FindOneP`)|pandoraShower|


---
## Glossary

### Parameters in FHiCL

### Variables

In `TTree vertex_tree`, variables prefix have the following meaning:
- `sss_*` 
- `trackstub_*` 
- `reco_*` 
- `sim_*` 
- `mctruth_*` 
