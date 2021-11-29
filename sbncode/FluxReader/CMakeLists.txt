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
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                ${ART_ROOT_IO_TFILE_SUPPORT}
                art_Persistency_Common canvas
                art_Persistency_Provenance canvas
                art_Utilities canvas
                ${MF_MESSAGELOGGER}
                nusimdata_SimulationBase
                larcore_Geometry_Geometry_service
                larcorealg_Geometry
                ${ROOT_GEOM}
                ${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME sbncode_FluxReader_FluxInterface
                  SOURCE	${INTERFACEFILES}
                  LIBRARIES
                  nusimdata_SimulationBase
                  # GENIE STUFF #####################################
                  xml2
                  LOG4CPP
                  GFwEG
                  GFwAlg
                  GFwMsg
                  GFwGHEP
                  GTlFlx
                  GFwParDat
                  dk2nuTree
                  #${DK2NUDATA}
                  ${DK2NUGENIE}
                  ####################################################
                  cetlib cetlib_except

                  ${FHICLCPP}
                  ${MF_MESSAGELOGGER}
                  ${CLHEP}
                  ${ROOT_BASIC_LIB_LIST}
                  ${Boost_SYSTEM_LIBRARY}
)

art_make( EXCLUDE FluxReader_source.cc ${INTERFACEFILES} FluxGeoFilter_module.cc
          LIB_LIBRARIES	nusimdata_SimulationBase
                        sbncode_FluxReader_FluxInterface
                        larcoreobj_SummaryData
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
                        cetlib cetlib_except

          MODULE_LIBRARIES sbncode_FluxReader
                        sbncode_FluxReader_FluxInterface
                        nusimdata_SimulationBase
                        larcoreobj_SummaryData
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        # GENIE STUFF #####################################
                        xml2
                        LOG4CPP
                        ${GALGORITHM}
                        ${GBASE}
                        ${GEVGCORE}
                        ${GEVGDRIVERS}
                        ${GFLUXDRIVERS}
                        ${GINTERACTION}
                        ${GMESSENGER}
                        ${GHEP}
                        ${DK2NUDATA}
                        ${DK2NUGENIE}
                        ####################################################
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_source()
install_fhicl()

add_subdirectory(job)
