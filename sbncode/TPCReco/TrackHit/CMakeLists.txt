simple_plugin(TrackAreaHit module
        art::Framework_Core
        art::Framework_Services_Registry
        art_root_io::TFileService_service
        art::Framework_Principal
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        lardataobj::RecoBase_AssnsDicts_dict
        lardataobj::RecoBase     
        lardataobj::RecoBase_dict        
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
        lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
        lardataalg::DetectorInfo
        larsim::Simulation
)

simple_plugin(TrackHitFilter module
        art::Framework_Core
        art::Framework_Services_Registry
        art_root_io::TFileService_service
        art::Framework_Principal
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        lardataobj::RecoBase_AssnsDicts_dict
        lardataobj::RecoBase     
        lardataobj::RecoBase_dict        
        lardataalg::DetectorInfo
        lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
)

simple_plugin(TrackHitDumper module
        art::Framework_Core
        art::Framework_Services_Registry
        art_root_io::TFileService_service
        art::Framework_Principal
        art::Persistency_Common
        art::Utilities canvas::canvas
        cetlib::cetlib cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        lardataobj::RecoBase_AssnsDicts_dict
        lardataobj::RecoBase
        lardataobj::RecoBase_dict
        lardataalg::DetectorInfo
        lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
        larcorealg::Geometry
        larcore::Geometry_Geometry_service
)

install_headers()
install_fhicl()
install_source()

