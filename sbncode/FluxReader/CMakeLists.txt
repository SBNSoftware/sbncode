cet_enable_asserts()
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{DK2NUDATA_INC} )
include_directories ( $ENV{DK2NUGENIE_INC} )
link_directories( $ENV{DK2NUDATA_LIB})

file(GLOB INTERFACEFILES *Interface* )

simple_plugin( FluxReader
               "source"
               sbncode_FluxReader
               ${ART_FRAMEWORK_PERSISTENCY_COMMON}
)

simple_plugin ( FluxGeoFilter "module"
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                art_root_io::TFileService_service
                art_root_io::tfile_support
                art::Persistency_Common canvas::canvas
                art::Persistency_Provenance canvas::canvas
                art::Utilities canvas::canvas
                messagefacility::MF_MessageLogger
                nusimdata::SimulationBase
                larcore_Geometry_Geometry_service
                larcorealg_Geometry
                ROOT::Geom
                ${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME sbncode_FluxReader_FluxInterface
                  SOURCE	${INTERFACEFILES}
                  LIBRARIES
                  nusimdata::SimulationBase
                  # GENIE STUFF #####################################
                  xml2
                  log4cpp::log4cpp
                  GENIE::GFwEG
                  GENIE::GFwAlg
                  GENIE::GFwMsg
                  GENIE::GFwGHEP
                  GENIE::GTlFlx
                  GENIE::GFwParDat
                  dk2nuTree
                  #${DK2NUDATA}
                  ${DK2NUGENIE}
                  ####################################################
                  cetlib::cetlib cetlib_except

                  fhiclcpp::fhiclcpp
                  messagefacility::MF_MessageLogger
                  CLHEP::CLHEP
                  ${ROOT_BASIC_LIB_LIST}

)

art_make( EXCLUDE FluxReader_source.cc ${INTERFACEFILES} FluxGeoFilter_module.cc
          LIB_LIBRARIES	nusimdata::SimulationBase
                        sbncode_FluxReader_FluxInterface
                        larcoreobj_SummaryData
                        art::Framework_Core
                        art::Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}

                        art_root_io::TFileService_service
                        art_root_io::tfile_support
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp

                        ${ROOT_BASIC_LIB_LIST}
                        cetlib::cetlib cetlib_except

          MODULE_LIBRARIES sbncode_FluxReader
                        sbncode_FluxReader_FluxInterface
                        nusimdata::SimulationBase
                        larcoreobj_SummaryData
                        art::Framework_Core
                        art::Framework_IO_Sources
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        art_root_io::TFileService_service
                        art_root_io::tfile_support
                        fhiclcpp::fhiclcpp
                        # GENIE STUFF #####################################
                        xml2
                        log4cpp::log4cpp
                        GENIE::GFwEG
                        GENIE::GFwAlg
                        GENIE::GFwMsg
                        GENIE::GFwGHEP
                        GENIE::GTlFlx
                        GENIE::GFwParDat
                        ${DK2NUDATA}
                        ${DK2NUGENIE}
                        ####################################################
                        ${ROOT_EVE_LIB_LIST}
                        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)
