cet_enable_asserts()
include_directories ( $ENV{GENIE}/include/GENIE )
include_directories ( $ENV{GENIE}/include )
find_package(dk2nudata REQUIRED)

file(GLOB INTERFACEFILES *Interface* )

simple_plugin( FluxReader
               "source"
               sbncode::FluxReader
               larcore::Geometry_Geometry_service
               nutools::EventGeneratorBase_CRY
               lardata::Utilities
)

simple_plugin ( FluxGeoFilter "module"
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

art_make_library( LIBRARY_NAME sbncode::FluxReader_FluxInterface
                  SOURCE	${INTERFACEFILES}
                  LIBRARIES
                art::Framework_IO_Sources
                  nusimdata::SimulationBase
               larcore::Geometry_Geometry_service
                  # GENIE STUFF #####################################
                  xml2

                  log4cpp::log4cpp
                  GENIE::GFwEG
                  GENIE::GFwAlg
                  GENIE::GFwMsg
                  GENIE::GFwGHEP
                  GENIE::GTlFlx
                  dk2nu::Tree
                  GENIE::GFwParDat
                  #${DK2NUDATA}
                  ${DK2NUGENIE}
                  ####################################################
                  cetlib::cetlib cetlib_except::cetlib_except
               lardata::Utilities

                  CLHEP::CLHEP
                  ROOT::Core ROOT::Tree
)

art_make( EXCLUDE FluxReader_source.cc ${INTERFACEFILES} FluxGeoFilter_module.cc
          LIB_LIBRARIES	nusimdata_SimulationBase
                        sbncode_FluxReader_FluxInterface
                        larcoreobj_SummaryData
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

          MODULE_LIBRARIES sbncode::FluxReader
                        sbncode::FluxReader_FluxInterface
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
                        cetlib cetlib_except
                        # GENIE STUFF #####################################
                        xml2
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

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)
