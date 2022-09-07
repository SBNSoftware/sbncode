Use `singlephoton_sbnd.fcl` to configure the module for reading sample files;

Use `run_singlePhoton_sbnd.fcl` to configure services;


### Labels for objects

*xLabel* is the alias for the actual label to locate some `<objects>`: 

*PandoraLabel*:
- `std::vector<recob::PFParticle>`
- `std::vector<recob::Cluster>`
- `std::vector<recob::Slice>`
- `std::vector<recob::Vertex>`

*TrackLabel*:
- ...
    TrackLabel:     "pandora"
    ShowerLabel:    "pandora"
    ParticleIDLabel:"pandoracalipidSCE"
    CaloLabel:      "pandoracaliSCE"
    FlashLabel:     "simpleFlashBeam"
    POTLabel:       "generator"
    input_param:    "optional" 

