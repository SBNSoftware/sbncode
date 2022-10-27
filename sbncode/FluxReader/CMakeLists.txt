# for classes_def.xml
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )

cet_enable_asserts()

file(GLOB INTERFACEFILES *Interface* )

cet_build_plugin( FluxReader art::source
          LIBRARIES
               sbncode_FluxReader
               larcore::Geometry_Geometry_service
               nutools::EventGeneratorBase_CRY
               lardata::Utilities
)

cet_build_plugin ( FluxGeoFilter art::module
            LIBRARIES
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                art::Framework_IO_Sources
                art_root_io::TFileService_service
                art_root_io::tfile_support
                art::Persistency_Common canvas::canvas
                art::Persistency_Provenance canvas::canvas
                art::Utilities canvas::canvas
                nusimdata::SimulationBase
                nugen::EventGeneratorBase_GENIE
                larcore::Geometry_Geometry_service
                lardata::Utilities
                larcorealg::Geometry
                ROOT::Geom
)

cet_make_library( LIBRARY_NAME sbncode_FluxReader_FluxInterface
                  SOURCE        ${INTERFACEFILES}
                  LIBRARIES
                  art::Framework_IO_Sources
                  nusimdata::SimulationBase
                  larcore::Geometry_Geometry_service
                  # GENIE STUFF #####################################
                  LibXml2::LibXml2
                  log4cpp::log4cpp
                  GENIE::GFwEG
                  GENIE::GFwAlg
                  GENIE::GFwMsg
                  GENIE::GFwGHEP
                  GENIE::GTlFlx
                  dk2nu::Tree
                  GENIE::GFwParDat
                  dk2nu::Genie
                  ####################################################
                  cetlib::cetlib cetlib_except::cetlib_except
                  lardata::Utilities
                  CLHEP::CLHEP
                  ROOT::Core ROOT::Tree
)

cet_make_library(
          SOURCE FluxReader.cxx
          LIBRARIES     nusimdata::SimulationBase
                        sbncode_FluxReader_FluxInterface
                        larcoreobj::SummaryData
                        art_root_io::TFileService_service
                        art_root_io::tfile_support
                        art::Framework_IO_Sources
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        ROOT::Core ROOT::Tree
                        cetlib::cetlib cetlib_except::cetlib_except
                        nutools::EventGeneratorBase_CRY
                        lardata::Utilities
)

cet_build_plugin(FluxReaderAna art::module LIBRARIES
                        sbncode_FluxReader
                        sbncode_FluxReader_FluxInterface
                        nusimdata::SimulationBase
                        larcoreobj::SummaryData
                        larcore::Geometry_Geometry_service
                        lardata::Utilities
                        art::Framework_Core
                        art::Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art::Framework_IO_Sources
                        art_root_io::TFileService_service
                        art_root_io::TFileService_service
                        art_root_io::tfile_support
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except::cetlib_except
                        # GENIE STUFF #####################################
                        LibXml2::LibXml2
                        log4cpp::log4cpp
                        GENIE::GFwEG
                        GENIE::GFwAlg
                        GENIE::GFwMsg
                        GENIE::GFwGHEP
                        GENIE::GTlFlx
                        GENIE::GFwParDat
                        ####################################################
                        ROOT::Core ROOT::Tree ROOT::EG
)

art_dictionary(DICTIONARY_LIBRARIES PRIVATE nusimdata::SimulationBase dk2nu::Tree)

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)
