simple_plugin(G4InfoReducer "module"
							larsim_Simulation
							${ART_FRAMEWORK_CORE}
							${ART_FRAMEWORK_PRINCIPAL}
							${ART_FRAMEWORK_SERVICES_REGISTRY}
							${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
							${ART_ROOT_IO_TFILESERVICE}
							art_Persistency_Common
							art_Persistency_Provenance
							art_Utilities canvas
							${MF_MESSAGELOGGER}
							${MF_UTILITIES}
							${CETLIB}
							${CLHEP}
							${ROOT_BASIC_LIB_LIST}
							BASENAME_ONLY)

install_fhicl()
install_source()
