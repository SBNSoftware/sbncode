# fix required for c7 and e19 with python 2
# also fine with c2 and e17
cet_add_compiler_flags(CXX -Wno-register )

# modernize finding python
find_package(Python COMPONENTS Interpreter Development)

include_directories(${Python_INCLUDE_DIRS})

art_make( LIBRARY_NAME fcl
	  LIB_LIBRARIES fhiclcpp::fhiclcpp
	  		cetlib::cetlib
			cetlib_except::cetlib_except
			Python::Python )
set_target_properties(fcl PROPERTIES PREFIX "" SUFFIX ".so")

install_source()
