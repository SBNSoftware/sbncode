
art_make_library(
		LIBRARIES
			lardataobj::RawData
			lardataobj::RecoBase
			lardata::Utilities
			fhiclcpp::fhiclcpp
			cetlib::cetlib
)
set(		MODULE_LIBRARIES
			larcorealg::Geometry
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardata::Utilities
			larevt::Filters
			lardataobj::RawData
			larevt::CalibrationDBI_IOVData
			larevt::CalibrationDBI_Providers
			lardataobj::RecoBase
			lardata::ArtDataHelper
			larreco::RecoAlg
			sbnobj::ICARUS_TPC
			sbnobj::Common_Utilities
      	    art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			messagefacility::headers
			fhiclcpp::fhiclcpp
			cetlib::cetlib
			ROOT::Geom
			ROOT::XMLIO
			ROOT::Gdml
			ROOT::FFTW

	)
cet_build_plugin(GaussHitFinderSBN art::module 
				 LIBRARIES ${MODULE_LIBRARIES}
				 larreco::HitFinder
				 larreco::CandidateHitFinderTool
				 larreco::PeakFitterTool
				)


install_headers()
install_fhicl()
install_source()

