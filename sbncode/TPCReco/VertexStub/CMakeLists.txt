cet_make_library(
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
cet_build_plugin(VertexChargeVacuum art::module
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
        lardataobj::RecoBase
        lardataalg::DetectorInfo
        sbncode_TPCReco_VertexStub
        larevt::SpaceCharge
)

cet_build_plugin(VertexStubTracker art::module
        LIBRARIES
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
        sbncode_TPCReco_VertexStub
        larevt::SpaceCharge
)

cet_build_plugin(PlaneStubMerge art::tool
        LIBRARIES
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
        sbncode_TPCReco_VertexStub
        larevt::SpaceCharge
)

cet_build_plugin(TwoPlaneStubMerge art::tool
        LIBRARIES
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
        sbncode_TPCReco_VertexStub
        larevt::SpaceCharge
)

add_subdirectory(config)

install_headers()
install_fhicl()
install_source()
