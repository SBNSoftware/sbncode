set(	MODULE_LIBRARIES
		larcorealg::Geometry
		larcore::Geometry_Geometry_service
		larsim::Simulation
		nug4::ParticleNavigation
		lardataobj::Simulation
		lardata::Utilities
		lardataalg::DetectorInfo
		larevt::Filters
		lardataobj::RawData
		larevt::CalibrationDBI_Providers
		nurandom::RandomUtils_NuRandomService_service
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art_root_io::tfile_support ROOT::Core
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art_root_io::TFileService_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		messagefacility::headers
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		CLHEP::Random
		ROOT::Geom
		ROOT::XMLIO
		ROOT::Gdml
		ROOT::Core
		FFTW3::FFTW3
		)

cet_build_plugin(AdjustSimForTrigger art::module LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(FilterSimEnergyDeposits art::module LIBRARIES ${MODULE_LIBRARIES})

add_subdirectory(fcl)

#install_headers()
install_fhicl()
install_source()
