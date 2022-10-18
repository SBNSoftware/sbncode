# for classes_def.xml
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
#set( PACKAGE CAFMaker )
FILE( GLOB src_files *.cxx )

add_subdirectory(RecoUtils)
add_subdirectory(bin)

include_directories( $ENV{SRPROXY_DIR}/include )
art_make_library( LIBRARY_NAME sbncode_CAFMaker
                  SOURCE ${src_files}
                  LIBRARIES
                  ifdh_art::IFDH_service
                  sbnanaobj::StandardRecord
                  sbnanaobj::StandardRecordFlat
                  art::Framework_Core
                  art::Framework_Services_Registry
                  art_root_io::TFileService_service
                  art::Framework_Principal
                  art::Persistency_Common
                  art::Utilities canvas::canvas
                  cetlib::cetlib cetlib_except::cetlib_except
                  ROOT::X3d
                  ROOT::EG
                  caf_RecoUtils
                  lardataobj::AnalysisBase
                  lardataobj::RecoBase
                  larcorealg::Geometry
                  larcore::Geometry_Geometry_service
                  larcorealg::GeoAlgo
                  larsim::MCCheater_BackTrackerService_service
                  nusimdata::SimulationBase
                  larsim::MCCheater_ParticleInventoryService_service
                  sbnobj::Common_CRT
                  sbnobj::Common_Reco
                  sbnobj::Common_Analysis
                  lardataalg::DetectorInfo
                  art::Framework_Services_System_TriggerNamesService_service
                  sbncode::Metadata_MetadataSBN_service
                  larsim::Utils
                  )


cet_build_plugin ( CAFMaker art::module
               LIBRARIES
               sbncode::CAFMaker
               sbnanaobj::StandardRecord
               larcorealg::Geometry
               ifdh_art::IFDH_service
               art::Framework_Core
               art::Framework_Principal
               art::Framework_Services_Registry
               art::Utilities canvas::canvas
               cetlib::cetlib cetlib_except::cetlib_except
               ROOT::Core ROOT::Tree
               art_root_io::RootDB
               hep_concurrency::hep_concurrency
               lardataobj::RecoBase
               nurandom::RandomUtils_NuRandomService_service
               lardata::DetectorClocksService
               BASENAME_ONLY
            )

build_dictionary( CAFMaker
                  DICTIONARY_LIBRARIES
                  sbnanaobj::StandardRecord
                  sbnanaobj::StandardRecord_dict
                  lardataobj::RecoBase
                  larcorealg::Geometry
                  art::Framework_Core
                  art::Framework_Principal
                  art::Persistency_Common
                  art::Persistency_Provenance
                  art::Utilities canvas::canvas
                  cetlib::cetlib cetlib_except::cetlib_except
                  )


install_headers()
install_fhicl()
install_source()
