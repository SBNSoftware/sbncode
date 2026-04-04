set(	MODULE_LIBRARIES
		lardata::DetectorClocksService
		larcore::headers
		sbncode::Utilities
		larcorealg::Geometry
		lardataobj::RawData
		lardataobj::Simulation
		sbnobj::ICARUS_PMT_Data
		)

cet_build_plugin(AdjustSimForTrigger art::module LIBRARIES PRIVATE ${MODULE_LIBRARIES})
cet_build_plugin(FilterSimEnergyDeposits art::module LIBRARIES PRIVATE ${MODULE_LIBRARIES})

add_subdirectory(fcl)

#install_headers()
install_fhicl()
install_source()
