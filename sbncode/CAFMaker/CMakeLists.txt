# for classes_def.xml
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
#set( PACKAGE CAFMaker )
FILE( GLOB src_files *.cxx )

add_subdirectory(RecoUtils)
add_subdirectory(bin)

art_make_library( LIBRARY_NAME sbncode_CAFMaker
                  SOURCE ${src_files}
                  LIBRARIES
                  sbnanaobj_StandardRecord
                  art::Framework_Core
                  art::Framework_Services_Registry
                  art_root_io::TFileService_service
                  ${ART_FRAMEWORK_SERVICES}
                  art::Framework_Principal
                  art::Persistency_Common
                  art::Utilities canvas::canvas
                  messagefacility::MF_MessageLogger
                  fhiclcpp::fhiclcpp
                  cetlib::cetlib cetlib_except
                  ${ROOT_EVE_LIB_LIST}
                  ${ROOT_X3d}
                  ${ROOT_BASIC_LIB_LIST}
                  caf_RecoUtils
                  lardataobj_AnalysisBase
                  lardataobj_RecoBase
                  larcorealg_Geometry
                  larcore_Geometry_Geometry_service
                  larcorealg_GeoAlgo
                  larsim_MCCheater_BackTrackerService_service
                  larsim_Utils
                  nusimdata::SimulationBase
                  larsim_MCCheater_ParticleInventoryService_service
                  sbnobj_Common_CRT
                  sbnobj_Common_Reco
                  sbnobj_Common_Analysis
		  lardataalg_DetectorInfo
                  art_Framework_Services_System_TriggerNamesService_service
                  sbncode_Metadata_MetadataSBN_service
                  )

cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
##include_directories( $ENV{IFDHC_FQ_DIR}/inc )
##include_directories( $ENV{IFDH_ART_FQ_DIR}/inc )

cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB )

simple_plugin ( CAFMaker module
               sbncode_CAFMaker
               sbnanaobj_StandardRecord
               ifdh_art::IFDH_service
               art::Framework_Core
               art::Framework_Principal
               art::Framework_Services_Registry
               ${ART_FRAMEWORK_SERVICES_BASIC}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               art::Utilities canvas::canvas
               messagefacility::MF_MessageLogger
               fhiclcpp::fhiclcpp
               cetlib::cetlib cetlib_except
               ${ROOT_BASIC_LIB_LIST}
               art_root_io_RootDB
               hep_concurrency
               nurandom::RandomUtils_NuRandomService_service
               BASENAME_ONLY
            )

build_dictionary( CAFMaker
                  DICTIONARY_LIBRARIES
                  sbnanaobj_StandardRecord
                  sbnanaobj_StandardRecord_dict
                  art::Framework_Core
                  art::Framework_Principal
                  art::Persistency_Common
                  art::Persistency_Provenance
                  art::Utilities canvas::canvas
                  cetlib::cetlib cetlib_except
                  messagefacility::MF_MessageLogger
                  
                  )


install_headers()
install_fhicl()
install_source()
