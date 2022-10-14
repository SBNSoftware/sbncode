art_make_library(
        SOURCE StubBuilder.cxx StubMergeAlgorithms.cxx
        LIBRARIES
	art::Framework_Core
	art::Framework_Principal
	art::Persistency_Common
	art::Persistency_Provenance
	canvas::canvas art::Utilities
	cetlib::cetlib cetlib_except::cetlib_except
	messagefacility::MF_MessageLogger
        lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
	sbnobj::Common_Reco
	lardataobj::AnalysisBase
        larcorealg::Geometry
	larreco::Calorimetry
	lardataobj::RecoBase
	lardataalg::DetectorInfo
        larevt::SpaceCharge
)
simple_plugin(VertexChargeVacuum module
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
	lardataobj::RecoBase
	lardataalg::DetectorInfo
	sbncode::TPCReco_VertexStub
        larevt::SpaceCharge
)

simple_plugin(VertexStubTracker module
	art::Framework_Core
	art::Framework_Services_Registry
	art_root_io::TFileService_service
	art::Framework_Principal
	art::Persistency_Common
	art::Persistency_Provenance
	art::Utilities canvas::canvas
	cetlib::cetlib cetlib_except::cetlib_except
	messagefacility::MF_MessageLogger
        lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
	sbnobj::Common_Reco
	lardataobj::AnalysisBase
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
	larreco::Calorimetry
	lardataalg::DetectorInfo
	sbncode::TPCReco_VertexStub
        larevt::SpaceCharge
)

simple_plugin(PlaneStubMerge tool
        art::Framework_Core
        art::Framework_Services_Registry
        art_root_io::TFileService_service
        art::Framework_Principal
        art::Persistency_Common
	art::Persistency_Provenance
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
	sbnobj::Common_Reco
	lardataobj::AnalysisBase
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
	larreco::Calorimetry
	lardataalg::DetectorInfo
	sbncode::TPCReco_VertexStub
        larevt::SpaceCharge
)

simple_plugin(TwoPlaneStubMerge tool
        art::Framework_Core
        art::Framework_Services_Registry
        art_root_io::TFileService_service
        art::Framework_Principal
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
	sbnobj::Common_Reco
	lardataobj::AnalysisBase
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
	larreco::Calorimetry
	lardataalg::DetectorInfo
	sbncode::TPCReco_VertexStub
        larevt::SpaceCharge
)

add_subdirectory(config)

install_headers()
install_fhicl()
install_source()
