add_library(sbncode_CosmicIdAlgs SHARED
    ApaCrossCosmicIdAlg.cc
    # CosmicIdAlg.cc
    # CpaCrossCosmicIdAlg.cc
    CrtHitCosmicIdAlg.cc
    CrtTrackCosmicIdAlg.cc
    FiducialVolumeCosmicIdAlg.cc
    GeometryCosmicIdAlg.cc
    PandoraT0CosmicIdAlg.cc
    StoppingParticleCosmicIdAlg.cc
    Util.cc
)

target_link_libraries(sbncode_CosmicIdAlgs
	larcorealg_Geometry
	larsim_Simulation
	lardataobj_Simulation
	lardata_Utilities
	lardataobj_RecoBase
	sbndcode_CRT
	sbndcode_CRTUtils
	sbnobj_Common_CRT
)

install(TARGETS sbncode_CosmicIdAlgs DESTINATION lib)

