cet_build_plugin(BlipReco3D art::module
        LIBRARIES
	art::Framework_Core
	art::Framework_Principal
	art::Framework_Services_Registry
	art::Persistency_Common
	art::Utilities canvas::canvas
	cetlib::cetlib cetlib_except::cetlib_except
	messagefacility::MF_MessageLogger
	lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
	sbnobj::Common_Reco
	lardataobj::AnalysisBase
	larcorealg::Geometry
	larreco::Calorimetry
	lardataalg::DetectorInfo
	lardataobj::RecoBase
	larsim::MCCheater_BackTrackerService_service
	larsim::MCCheater_ParticleInventoryService_service
	nusimdata::SimulationBase
)
install_headers()
install_fhicl()
install_source()
