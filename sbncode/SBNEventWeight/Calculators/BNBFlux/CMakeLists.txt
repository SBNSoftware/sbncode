include_directories($ENV{GENIE_INC}/GENIE)

art_make(
		LIBRARY_NAME sbncode_SBNEventWeight_Calculators_BNBFlux
		LIB_LIBRARIES
		sbncode_SBNEventWeight_Base
		nugen_EventGeneratorBase_GENIE
		nugen_NuReweight_art
		nugen_NuReweight
		nurandom_RandomUtils_NuRandomService_service
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_PERSISTENCY_PROVENANCE}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${FHICLCPP}
		cetlib_except
		)

install_headers()
install_fhicl()
install_source()

add_subdirectory(BNBPrimaryHadron)
