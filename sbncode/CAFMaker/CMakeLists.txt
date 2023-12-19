# for classes_def.xml
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
#set( PACKAGE CAFMaker )

include_directories( $ENV{SRPROXY_DIR}/include )

add_subdirectory(RecoUtils)
add_subdirectory(bin)

art_make_library( LIBRARY_NAME sbncode_CAFMaker
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
                  CLHEP::Random
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
                  sbnobj::SBND_CRT
                  lardataalg::DetectorInfo
                  art::Framework_Services_System_TriggerNamesService_service
                  sbncode_Metadata_MetadataSBN_service
                  larsim::Utils
		  larevt::SpaceCharge
                  systematicstools::interface
                  systematicstools::interpreters
                  systematicstools::utility
                  ${GENIE_LIB_LIST}
                  nugen::EventGeneratorBase_GENIE
                  )

# this can be removed once icaruscode moves to cetmodules
cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB )

cet_build_plugin ( CAFMaker art::module
               LIBRARIES
               sbncode_CAFMaker
               sbnanaobj::StandardRecord
               larcorealg::Geometry
               ifdh_art::IFDH_service
               art::Framework_Core
               art::Framework_Principal
               art::Framework_Services_Registry
               art::Utilities canvas::canvas
               cetlib::cetlib cetlib_except::cetlib_except
               CLHEP::Random
               ROOT::Core ROOT::Tree
               art_root_io::RootDB
               hep_concurrency::hep_concurrency
               lardataobj::RecoBase
               nurandom::RandomUtils_NuRandomService_service
               lardata::DetectorClocksService
               BASENAME_ONLY
            )

art_dictionary( CAFMaker
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
