# Add the subfolder containing the deconvolution tools
art_make( 
          MODULE_LIBRARIES larcorealg_Geometry
            larcore_Geometry_Geometry_service
            lardata_Utilities
            larevt_Filters
            lardataobj_RawData
            larevt_CalibrationDBI_IOVData
            larevt_CalibrationDBI_Providers
            lardataobj_RecoBase
            lardata_ArtDataHelper
            larsim_Simulation nug4::ParticleNavigation 
            lardataobj_Simulation
            nusimdata::SimulationBase
            lardata_Utilities
            lardataobj_MCBase
            larcoreobj_SummaryData
            nusimdata::SimulationBase
            nurandom::RandomUtils_NuRandomService_service
            lardataobj_MCBase
            CLHEP::CLHEP
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
            
            fhiclcpp::fhiclcpp
            cetlib::cetlib
            ROOT::FFTW
            ${ROOT_BASIC_LIB_LIST}
	    ROOT::EG
        )

install_headers()
install_fhicl()
install_source()

