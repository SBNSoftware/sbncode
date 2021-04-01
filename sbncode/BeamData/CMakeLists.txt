if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
cet_add_compiler_flags(CXX -Wno-dynamic-exception-spec -Wno-return-type-c-linkage )
endif()

include_directories ( $ENV{IFDHC_FQ_DIR}/inc )

find_library( UB_DATA_TYPES NAMES ubdata_types PATHS $ENV{UBOONEDAQ_DATATYPES_LIB} )
cet_make(LIBRARY_NAME BeamUtil LIBRARIES ${UB_DATA_TYPES} ${ROOT_SPECTRUM} ${ROOT_BASIC_LIB_LIST} EXCLUDE BeamData_module.cc BeamDataQualityFilter_module.cc)

art_make( EXCLUDE BeamDataQualityFilter_module.cc
          BASENAME_ONLY
          MODULE_LIBRARIES
          lardata_Utilities
          larevt_Filters
          lardataobj_RawData
			    larcoreobj_SummaryData
			    BeamUtil
			    BeamDAQ
			    ${UB_DATA_TYPES}
			    ${Boost_SERIALIZATION_LIBRARY}
          nurandom_RandomUtils_NuRandomService_service
          ${ART_FRAMEWORK_CORE}
          ${ART_FRAMEWORK_PRINCIPAL}
          ${ART_FRAMEWORK_SERVICES_REGISTRY}
          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
          ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
          art_Persistency_Common canvas
          art_Persistency_Provenance canvas
          art_Utilities canvas
			    ifdh
          ${MF_MESSAGELOGGER}
          ${FHICLCPP}
          cetlib cetlib_except
          ${CLHEP}
          ${ROOT_BASIC_LIB_LIST}
		      
          LIB_LIBRARIES
			    ${UB_DATA_TYPES}
			    ${ROOT_SPECTRUM}
          ${ROOT_BASIC_LIB_LIST}
        )


# install_headers()
install_fhicl()
install_source()
